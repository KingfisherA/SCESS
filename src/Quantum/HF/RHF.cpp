#include<Quantum/Q_PCH.h>
//
#include <chrono>
#include <algorithm>
#include <Eigen/Core>
//
#include "HF.h"
#include "DIIS.h"
#include "Quantum/Basis/BasisSet.h"
#include "Quantum/Basis/CGF.h"
#include "Quantum/Integrals/Integrals.h"
#include "Quantum/Molecule.h"
#include "Quantum/Integrals/IntegralTransforms.h"
#include "Dahlia/Tools/Threading.h"


//
//
namespace SCESS::HF {
    using namespace Eigen;
    //
    //
    MatrixXd ConstructGRHF(const MatrixXd &D,
                           const TEITensor &eri);
    MatrixXd ProjectOrbitalDensity(const BasisSet &set1,
                                   const BasisSet &set2,
                                   const Integrals::MolInts &mints,
                                   const HF_Wfn &rhfWfn);
    MatrixXd ConstructDensityRHF(const MatrixXd &Coeffs,
                                 Index nElecPairs);
    double CalculateDelta(const MatrixXd &oldDensity,
                          const MatrixXd &density);
    //
    HF_Wfn RHF(const Molecule &molecule,
              const BasisSet &basis,
              double convergeLimit,
              size_t scfMaxIterations,
              Integrals::MolInts *mintsOut,
              const DH::Thread *thread,
              std::stringstream *console) {
        using namespace std::chrono;

        if (molecule.GetNoElectrons() % 2 != 0)
            throw std::runtime_error("Closed-shell RHF requires an even number of electrons");

        std::cout << "\nBeginning integral evaluation...\n";
        auto integralStartTP = high_resolution_clock::now();
        Integrals::MolInts integrals = Integrals::EvaluateMolecularIntegrals(basis, molecule.GetAtoms(), thread);
        auto integralEndTP = high_resolution_clock::now();


        if (console)
            *console << "Integral evaluation completed in " <<
                    duration_cast<milliseconds>(integralEndTP - integralStartTP).count() << " milliseconds\n";


        HF_Wfn wfn = RHF(integrals,
                        molecule.GetNoElectrons() / 2,
                        molecule.CalculateVnn(),
                        convergeLimit,
                        scfMaxIterations,
                        thread,
                        console);

        if (mintsOut)
            (*mintsOut) = std::move(integrals);

        return wfn;
    }
    //
    HF_Wfn RHF(const Integrals::MolInts &mints,
              uint32_t nElectronPairs,
              double Vnn,
              double convergeLimit,
              size_t scfMaxIterations,
              const DH::Thread *thread,
              std::stringstream *console) {
        if (console)
            *console << "---Restricted Hartree Fock---\nUsing " << mints.S.rows() << " basis set functions\n";

        if (thread)
            if (thread->ShouldStop())
                return {};

        const MatrixXd &S = mints.S;
        const MatrixXd &T = mints.T;
        const MatrixXd &V = mints.V;
        const TEITensor &eri = mints.eri; //Two-electron repulsion integral tensor

        Index dimensions = S.rows();
        double delta = 1.0;
        double E_RHF = 0.0;
        double lastEnergy = DBL_MAX;
        bool converged = false;
        //The RMS of the orbital gradient (FDS - SDF) required for convergence. The energy error is quadratic in the
        //orbital gradient so pair it with the square root of the energy limit
        const double commutatorLimit = std::sqrt(convergeLimit);

        //Core one-electron elements of the Hamiltonian
        MatrixXd h = T + V;

        //Canonical orthogonalisation transformation matrix
        SelfAdjointEigenSolver<MatrixXd> S_eig(S);
        MatrixXd s = S_eig.eigenvalues().real().cwiseInverse().cwiseSqrt().asDiagonal();
        MatrixXd U = S_eig.eigenvectors().real();
        MatrixXd X = U * s;

        MatrixXd F, C;
        //Density matrix
        MatrixXd D = MatrixXd::Zero(dimensions, dimensions);
        VectorXd Ep;

        //Diagonalise a Fock matrix in the orthonormal basis and rebuild the orbitals and density from it
        auto DiagonaliseFock = [&](const MatrixXd &fock) {
            //Transform the Fock matrix in the orthonormal basis set
            MatrixXd F_Prime = X.transpose() * fock * X;

            //Diagonalise the transformed Fock matrix
            SelfAdjointEigenSolver<MatrixXd> eigensolver2(F_Prime);
            Ep = eigensolver2.eigenvalues().real();

            //Transform coefficients back into normal basis
            C = X * eigensolver2.eigenvectors().real();

            //Recompute the density matrix
            D = ConstructDensityRHF(C, nElectronPairs);
        };

        if (console)
            *console << "Beginning SCF iterations...\n";

        //DIIS is turned on at the fourth iteration, discarding the first two (poorest) Fock matrices
        FockDIIS diis(6);

        Index scfIter = 1;
        for (; scfIter < scfMaxIterations + 1; scfIter++) {
            if (thread)
                if (thread->ShouldStop())
                    return {};

            //Generate the Fock matrix
            F = h + ConstructGRHF(D, eri);
            F = (F + F.transpose()) * 0.5;

            //The energy of the current density, evaluated with the Fock matrix built from that same density
            lastEnergy = E_RHF;
            E_RHF = CalculateEnergyRHF(D, h, F) + Vnn;

            //Build the DIIS error vector (FDS-SDF)
            MatrixXd DIIS_Error = F * D * S - S * D * F;
            DIIS_Error = X.transpose() * DIIS_Error * X;
            delta = sqrt(DIIS_Error.cwiseSquare().mean());

            if (console)
                *console << "   RHF(" << scfIter << "): E = " << E_RHF << "      dE = " << delta << '\n';

            if (scfIter > 1 && delta < commutatorLimit && abs(E_RHF - lastEnergy) < convergeLimit) {
                converged = true;
                break;
            }

            diis.Add({F}, {DIIS_Error});

            MatrixXd F_Extrapolated = F;
            if (scfIter >= 4) {
                if (scfIter == 4) {
                    diis.DiscardOldest(2);
                    if (console)
                        *console << "Turning on DIIS \n";
                }
                F_Extrapolated = diis.Extrapolate()[0];
            }

            DiagonaliseFock(F_Extrapolated);
        }

        //Canonical orbitals of the (non-extrapolated) Fock matrix of the final density so that the returned orbitals,
        //orbital energies and Fock matrix are consistent with each other
        DiagonaliseFock(F);

        if (console) {
            if (converged)
                *console << "---- RHF procedure converged! ----\nRHF converged after " <<
                        scfIter << " iterations.\nE(RHF) = " << E_RHF << " hartrees\n\n\n";
            else
                *console << "!!! RHF procedure failed to converge !!!\n\n\n";
        }

        VectorXi OccAlpha = VectorXi::Zero(h.rows());
        VectorXi OccBeta = VectorXi::Zero(h.rows());
        for (uint32_t i = 0; i < nElectronPairs; i++) {
            OccAlpha(i) = 1;
            OccBeta(i) = 1;
        }


        return HF_Wfn{
            E_RHF,
            F,
            {},
            D,
            {},
            C,
            {},
            Ep,
            {},
            OccAlpha,
            OccBeta,
            false
        };
    }
    //
    MatrixXd ConstructGRHF(const MatrixXd &D,
                           const TEITensor &eri) {
        //G = J - K / 2 for the closed-shell density D
        return Integrals::Transforms::CoulombMatrix(eri, D) - 0.5 * Integrals::Transforms::ExchangeMatrix(eri, D);
    }
    //
    MatrixXd ProjectOrbitalDensity(const BasisSet &set1,
                                   const BasisSet &set2,
                                   const Integrals::MolInts &mints,
                                   const HF_Wfn &rhfWfn) {
        MatrixXd S21 = Integrals::CalculateMixedOverlapMatrix(set1, set2);
        MatrixXd S22Inv = mints.S.completeOrthogonalDecomposition().pseudoInverse();

        MatrixXd Density = S22Inv * (S21 * rhfWfn.D_Alpha * S21.transpose()) * S22Inv;
        Density.conservativeResize(mints.S.rows(), mints.S.cols());
        return Density;
    }
    //
    MatrixXd ConstructDensityRHF(const MatrixXd &Coeffs,
                                 Index nElecPairs) {
        const auto occupied = Coeffs.leftCols(nElecPairs);
        return 2.0 * occupied * occupied.transpose();
    }
    //
    double CalculateDelta(const MatrixXd &oldDensity,
                          const MatrixXd &density) {
        double delta = 0.0;
        for (int i = 0; i < density.rows(); i++) {
            for (int j = 0; j < density.rows(); j++) {
                delta += std::pow(density(i, j) - oldDensity(i, j), 2);
            }
        }

        return std::sqrt(delta);
    }
    //
    double CalculateEnergyRHF(const MatrixXd &Density,
                              const MatrixXd &Hcore, const MatrixXd &Fock) {
        double result = 0.0;

        for (int i = 0; i < Fock.rows(); i++) {
            for (int j = 0; j < Fock.rows(); j++) {
                result += Density(i, j) * (Hcore(i, j) + Fock(i, j));
            }
        }

        return 0.5 * result;

        // MatrixXd HPlusF = Hcore + Fock;
        // Map<const VectorXd> vHPlusF(HPlusF.data(), HPlusF.size());
        // Map<const VectorXd> vDensity(Density.data(), Density.size());
        //
        // return 0.5 * vDensity.dot(vHPlusF);
    }
    //
    double MP2_RHF_ClosedShellAO(const HF_Wfn &wfn,
                                 const TEITensor &aoERI,
                                 const DH::Thread *thread,
                                 std::stringstream *console) {
        const VectorXi &occ = wfn.Occ_Alpha;
        const VectorXd &ep = wfn.Ep_Alpha;
        const Index nOcc = std::distance(occ.begin(), std::find(occ.begin(), occ.end(), 0));
        const Index nVirt = occ.size() - nOcc;

        if (console) {
            *console << std::setprecision(12);
            *console << "--- Moller Plesset Order 2 Correction ---\n";
        }

        //Only (a i | b j) is needed. The last index is transformed first, so starting with an occupied index the
        //largest intermediate is N^3 x nOcc rather than N^4
        const MatrixXd C_Occ = wfn.C_Alpha.leftCols(nOcc), C_Virt = wfn.C_Alpha.rightCols(nVirt);
        const TEITensor aibj = Integrals::Transforms::TransformTEI(aoERI, C_Virt, C_Occ, C_Virt, C_Occ);
        if (thread && thread->ShouldStop())
            return 0.0;

        double E_MP2 = 0.0;
        for (Index j = 0; j < nOcc; j++)
            for (Index b = 0; b < nVirt; b++)
                for (Index i = 0; i < nOcc; i++)
                    for (Index a = 0; a < nVirt; a++) {
                        const double iajb = aibj(a, i, b, j), ibja = aibj(b, i, a, j);
                        E_MP2 -= iajb * (2.0 * iajb - ibja) /
                                (ep[nOcc + a] + ep[nOcc + b] - ep[i] - ep[j]);
                    }

        if (console) {
            *console << "MP2 Correction Energy: " << E_MP2 << '\n';
            *console << "E(Hartree Fock) + E(MP2) = " << wfn.E_HF + E_MP2 << " hartrees\n\n\n";
        }
        return E_MP2;
    }
    //
    double MP2_RHF_ClosedShell(const HF_Wfn &wfn,
                               const TEITensor &moERI,
                               const DH::Thread *thread,
                               std::stringstream *console) {
        const VectorXi &occ = wfn.Occ_Alpha;
        const VectorXd &ep = wfn.Ep_Alpha;

        double E_MP2 = 0.0;

        //Number of doubly-occupied molecular orbitals
        Index nOcc = std::distance(occ.begin(), std::find(occ.begin(), occ.end(), 0));

        if (console) {
            *console << std::setprecision(12);
            *console << "--- Moller Plesset Order 2 Correction ---\n";
        }

        //i and j run over the occupied orbitals
        for (Index i = 0; i < nOcc; i++)
            for (Index j = 0; j < nOcc; j++) {
                if (thread)
                    if (thread->ShouldStop())
                        return 0.0;

                //a and b run over the virtual orbitals
                for (Index a = nOcc; a < occ.size(); a++)
                    for (Index b = nOcc; b < occ.size(); b++) {
                        //Note that integrals are in chemist's notation
                        double intIJAB = moERI(i, a, j, b);
                        double intIJBA = moERI(i, b, j, a);
                        double epsilonIJAB = ep[a] + ep[b] - ep[i] - ep[j];

                        //MP2 formula below already integrates spin out for a closed shell
                        //configuration
                        E_MP2 -= intIJAB * (2.0 * intIJAB - intIJBA) / epsilonIJAB;
                    }
            }

        if (console) {
            *console << "MP2 Correction Energy: " << E_MP2 << '\n';
            *console << "E(Hartree Fock) + E(MP2) = " << wfn.E_HF + E_MP2 << " hartrees\n\n\n";
        }

        return E_MP2;
    }
}
