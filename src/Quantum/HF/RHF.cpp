#include<Quantum/Q_PCH.h>
//
#include <chrono>
#include <algorithm>
#include <Eigen/Core>
//
#include "HF.h"
#include "Quantum/Basis/BasisSet.h"
#include "Quantum/Basis/CGF.h"
#include "Quantum/Integrals/Integrals.h"
#include "Quantum/Molecule.h"
#include "Quantum/Integrals/IntegralTransforms.h"
#include "GraphicsEngine/Tools/Threading.h"


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
              const GEng::Thread *thread,
              std::stringstream *console) {
        using namespace std::chrono;

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
              const GEng::Thread *thread,
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

        if (console)
            *console << "Beginning SCF iterations...\n";

        std::vector<MatrixXd> DIIS_Errors;
        DIIS_Errors.reserve(7);

        //Fock DIIS matrices
        std::vector<MatrixXd> F_DIIS;
        F_DIIS.reserve(7);

        Index scfIter = 1;
        for (; scfIter < scfMaxIterations + 1; scfIter++) {
            if (thread)
                if (thread->ShouldStop())
                    return {};

            //Generate the Fock matrix
            F = h + ConstructGRHF(D, eri);
            F = (F + F.transpose()) * 0.5;

            //Build the DIIS error vector (FDS-SDF)
            MatrixXd DIIS_Error = F * D * S - S * D * F;
            DIIS_Error = X.transpose() * DIIS_Error * X;
            DIIS_Errors.push_back(DIIS_Error);
            F_DIIS.push_back(F);
            delta = sqrt(DIIS_Error.cwiseSquare().mean());

            if (scfIter >= 4) {
                //Turn on DIIS

                if (scfIter == 4) {
                    F_DIIS.erase(F_DIIS.begin(), F_DIIS.begin() + 2);
                    DIIS_Errors.erase(DIIS_Errors.begin(), DIIS_Errors.begin() + 2);
                    if (console)
                        *console << "Turning on DIIS \n";
                }

                Index DIIS_Size = F_DIIS.size();

                if (DIIS_Size > 6) {
                    F_DIIS.erase(F_DIIS.begin());
                    DIIS_Errors.erase(DIIS_Errors.begin());
                    DIIS_Size -= 1;
                }

                MatrixXd ErrorMat = MatrixXd::Constant(DIIS_Size + 1, DIIS_Size + 1, -1.0);
                ErrorMat(ErrorMat.rows() - 1, ErrorMat.cols() - 1) = 0.0;

                double maxVal = 0.0;
                for (Index i = 0; i < DIIS_Errors.size(); i++)
                    for (Index j = i; j < DIIS_Errors.size(); j++) {
                        ErrorMat(i, j) = ErrorMat(j, i) = DIIS_Errors[i].cwiseProduct(DIIS_Errors[j]).sum();
                        maxVal = std::max(maxVal, abs(ErrorMat(i, j)));
                    }

                ErrorMat.block(0, 0, ErrorMat.rows() - 1, ErrorMat.cols() - 1) /= maxVal;

                VectorXd Residual = VectorXd::Zero(DIIS_Size + 1);
                Residual(Residual.rows() - 1) = -1.0;

                VectorXd DIISCoeffs = ErrorMat.lu().solve(Residual);

                F = MatrixXd::Zero(F.rows(), F.cols());
                for (Index i = 0; i < F_DIIS.size(); i++)
                    F += DIISCoeffs(i) * F_DIIS[i];
            }

            //Transform the Fock matrix in the orthonormal basis set
            MatrixXd F_Prime = X.transpose() * F * X;

            //Diagonalise the transformed Fock matrix
            SelfAdjointEigenSolver<MatrixXd> eigensolver2(F_Prime);
            Ep = eigensolver2.eigenvalues().real();

            //Transform coefficients back into normal basis
            const MatrixXd &C_Prime = eigensolver2.eigenvectors().real();
            C = X * C_Prime;

            //Recompute the density matrix
            D = ConstructDensityRHF(C, nElectronPairs);

            //Calculate delta change to the density matrix
            // const MatrixXd DensityMatriciesDelta = Density - OldDensity;
            // delta = std::sqrt(DensityMatriciesDelta.cwiseProduct(DensityMatriciesDelta).sum());

            lastEnergy = E_RHF;
            E_RHF = CalculateEnergyRHF(D, h, F) + Vnn;

            if (delta < 1E-3 && abs(E_RHF - lastEnergy) < convergeLimit)
                break;

            //Calculate the current energy
            if (console)
                *console << "   RHF(" << scfIter << "): E = " << E_RHF << "      dE = " << delta << '\n';
        }

        E_RHF = CalculateEnergyRHF(D, h, F) + Vnn;
        if (console) {
            if (scfIter < scfMaxIterations)
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
        MatrixXd G = MatrixXd::Zero(D.rows(), D.rows());

        for (Index i = 0; i < G.rows(); i++)
            for (Index j = 0; j < G.rows(); j++)
                for (Index k = 0; k < G.rows(); k++)
                    for (Index l = 0; l < G.rows(); l++) {
                        //Coulomb integral; <ik | jl> = (ij | kl)
                        double J = eri(i, j, k, l);
                        //Exchange integral; <ij | kl> = (ik | jl)
                        double K = eri(i, k, j, l);
                        G(i, j) += D(k, l) * (J - 0.5 * K);
                    }

        return G;
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
        MatrixXd Density = MatrixXd::Zero(Coeffs.rows(), Coeffs.rows());
        for (Index i = 0; i < Density.rows(); i++) {
            for (Index j = 0; j < Density.rows(); j++) {
                for (Index k = 0; k < nElecPairs; k++) {
                    Density(i, j) += 2.0 * Coeffs(i, k) * Coeffs(j, k);
                }
            }
        }

        return Density;
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
                                 const GEng::Thread *thread,
                                 std::stringstream *console) {
        //<ab | ij > [MO Basis] = <mu kappa | nu lambda> [AO Basis]
        const TEITensor moTEI = Integrals::Transforms::TransformTEI(aoERI, wfn.C_Alpha);
        return MP2_RHF_ClosedShell(wfn, moTEI, thread, console);
    }
    //
    double MP2_RHF_ClosedShell(const HF_Wfn &wfn,
                               const TEITensor &moERI,
                               const GEng::Thread *thread,
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
