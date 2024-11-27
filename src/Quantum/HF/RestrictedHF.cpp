#include<pch.h>

#include <chrono>
#include <algorithm>
#include <Eigen/Core>

#include "GraphicsEngine/Utilities/Threading.h"

#include "HartreeFock.h"
#include "Quantum/Basis/Gaussian.h"
#include "Quantum/Integrals/Integrals.h"
#include "Quantum/Molecule.h"
#include "Quantum/Integrals/IntegralTransforms.h"

//
//
namespace SCESS {
    using namespace Eigen;

    uRef<RHFWavefunction> HartreeFock::RHF(const Molecule &molecule,
                                           const Molecule *lastMolecule,
                                           const RHFWavefunction *lastResult,
                                           const double convergeLimit,
                                           const size_t scfMaxIterations,
                                           const GEng::Thread *thread,
                                           std::stringstream *console,
                                           const std::string_view basisSet) {
        using namespace std::chrono;


        if (console)
            *console << "---Restricted Hartree Fock---\nUsing " << basisSet << " basis set with " << molecule.
                    GetOrbitals().size() << " basis set functions \nBeginning integral evaluation...\n";

        auto integralStartTP = high_resolution_clock::now();
        ConfigIntegrals integrals = Integrals::CalculateConfiguration(molecule, thread);
        auto integralEndTP = high_resolution_clock::now();

        if (thread)
            if (thread->ShouldStop())
                return nullptr;

        const MatrixXd &S = integrals.S;
        const MatrixXd &T = integrals.T;
        const MatrixXd &V = integrals.V;
        const TEITensor &eri = integrals.eri; //Two-electron repulsion integral tensor

        if (console)
            *console << "Integral evaluation completed in " <<
                    duration_cast<milliseconds>(integralEndTP - integralStartTP).count() << " milliseconds\n";

        const Eigen::Index dimensions = S.rows();
        double delta = 1.0;
        double Vnn = molecule.CalculateNuclearRepulsion();
        double energy = 0.0, lastEnergy = DBL_MAX;

        MatrixXd Hcore = T + V;

        SelfAdjointEigenSolver<MatrixXd> eigensolver;
        eigensolver.compute(S);
        const VectorXd &SEigenValues = eigensolver.eigenvalues().real();
        const MatrixXd &SEigenVects = eigensolver.eigenvectors().real();
        MatrixXd SValsInvSqrt = SEigenValues.cwiseInverse().cwiseSqrt().asDiagonal();
        MatrixXd TransformMat = SEigenVects * SValsInvSqrt * SEigenVects.transpose();
        TransformMat = (TransformMat + TransformMat.transpose()) * 0.5;

        MatrixXd Fock, Coeffs, Density;
        VectorXd Energies;

        //Project the last result from the smaller into the larger basis set
        if (lastResult && lastMolecule) {
            MatrixXd S21 = Integrals::CalculateMixedOverlapMatrix(molecule, *lastMolecule);
            MatrixXd S22Inv = S.completeOrthogonalDecomposition().pseudoInverse();

            Density = S22Inv * (S21 * lastResult->Density * S21.transpose()) * S22Inv;
            Density.conservativeResize(dimensions, dimensions);
        } else
            Density = MatrixXd::Zero(dimensions, dimensions);

        if (console)
            *console << "Beginning SCF iterations...\n";


        std::vector<MatrixXd> DIIS_Errors;
        DIIS_Errors.reserve(7);
        std::vector<MatrixXd> FockMatrices;
        FockMatrices.reserve(7);

        Eigen::Index scfIter = 1;
        for (; scfIter < scfMaxIterations + 1; scfIter++) {
            if (thread)
                if (thread->ShouldStop())
                    return nullptr;

            //Generate the Fock matrix
            Fock = Hcore + ConstructGRHF(Density, eri);
            Fock = (Fock + Fock.transpose()) * 0.5;

            //Build the DIIS error vector (FDS-SDF)
            MatrixXd DIIS_Error = Fock * Density * S - S * Density * Fock;
            DIIS_Error = TransformMat.transpose() * DIIS_Error * TransformMat;
            DIIS_Errors.push_back(DIIS_Error);
            FockMatrices.push_back(Fock);
            delta = sqrt(DIIS_Error.cwiseSquare().mean());

            if (scfIter >= 4) {
                //Turn on DIIS

                if (scfIter == 4) {
                    FockMatrices.erase(FockMatrices.begin(), FockMatrices.begin() + 2);
                    DIIS_Errors.erase(DIIS_Errors.begin(), DIIS_Errors.begin() + 2);
                    if (console)
                        *console << "Turning on DIIS \n";
                }

                Eigen::Index DIIS_Size = FockMatrices.size();

                if (DIIS_Size > 6) {
                    FockMatrices.erase(FockMatrices.begin());
                    DIIS_Errors.erase(DIIS_Errors.begin());
                    DIIS_Size -= 1;
                }

                MatrixXd ErrorMat = MatrixXd::Constant(DIIS_Size + 1, DIIS_Size + 1, -1.0);
                ErrorMat(ErrorMat.rows() - 1, ErrorMat.cols() - 1) = 0.0;

                double maxVal = 0.0;
                for (Eigen::Index i = 0; i < DIIS_Errors.size(); i++)
                    for (Eigen::Index j = i; j < DIIS_Errors.size(); j++) {
                        ErrorMat(i, j) = ErrorMat(j, i) = DIIS_Errors[i].cwiseProduct(DIIS_Errors[j]).sum();
                        maxVal = std::max(maxVal, abs(ErrorMat(i, j)));
                    }

                ErrorMat.block(0, 0, ErrorMat.rows() - 1, ErrorMat.cols() - 1) /= maxVal;

                VectorXd Residual = VectorXd::Zero(DIIS_Size + 1);
                Residual(Residual.rows() - 1) = -1.0;

                VectorXd DIISCoeffs = ErrorMat.lu().solve(Residual);

                Fock = MatrixXd::Zero(Fock.rows(), Fock.cols());
                for (Eigen::Index i = 0; i < FockMatrices.size(); i++)
                    Fock += DIISCoeffs(i) * FockMatrices[i];
            }

            //Transform the Fock matrix in the orthonormal basis set
            MatrixXd FockPrime = TransformMat.transpose() * Fock * TransformMat;

            //Diagonalise the transformed Fock matrix
            SelfAdjointEigenSolver<MatrixXd> eigensolver2(FockPrime);
            //eigensolver2.compute(FockPrime);
            Energies = eigensolver2.eigenvalues().real();

            //Transform coefficients back into normal basis
            const MatrixXd &CPrime = eigensolver2.eigenvectors().real();
            Coeffs = TransformMat * CPrime;

            //Recompute the density matrix
            Density = ConstructDensityRHF(Coeffs, molecule.GetNoElectrons());

            //Calculate delta change to the density matrix
            // const  MatrixXd DensityMatriciesDelta = Density - OldDensity;
            // delta = std::sqrt(DensityMatriciesDelta.cwiseProduct(DensityMatriciesDelta).sum());

            lastEnergy = energy;
            energy = CalculateEnergyRHF(Density, Hcore, Fock) + Vnn;

            if (delta < 1E-3 && abs(energy - lastEnergy) < convergeLimit)
                break;

            //Calculate the current energy
            if (console)
                *console << "   RHF(" << scfIter << "): E = " << energy << "      dE = " << delta << '\n';
        }

        energy = CalculateEnergyRHF(Density, Hcore, Fock) + Vnn;
        if (console) {
            if (scfIter < scfMaxIterations)
                *console << "---- RHF procedure converged! ----\nRHF converged after " <<
                        scfIter << " iterations.\nE(RHF) = " << energy << " hartrees\n\n\n";
            else
                *console << "!!! RHF procedure failed to converge !!!\n\n\n";
        }

        VectorXi Occupation = VectorXi::Zero(Hcore.rows());

        for (int64_t i = ((int64_t) molecule.GetNoElectrons() / 2) - 1; i >= 0; --i) {
            Occupation(i) = 2;
        }


        return std::make_unique<RHFWavefunction>(Fock, Hcore, Density, Coeffs, Energies, Occupation, integrals, energy);
    }

    MatrixXd HartreeFock::ConstructGRHF(const MatrixXd &Density,
                                        const TEITensor &teiTensor) {
        MatrixXd G = MatrixXd::Zero(Density.rows(), Density.rows());

        for (Eigen::Index i = 0; i < G.rows(); i++)
            for (Eigen::Index j = 0; j < G.rows(); j++)
                for (Eigen::Index k = 0; k < G.rows(); k++)
                    for (Eigen::Index l = 0; l < G.rows(); l++) {
                        //Coulomb integral; <ik | jl> = (ij | kl)
                        double J = teiTensor(i, j, k, l);
                        //Exchange integral; <ij | kl> = (ik | jl)
                        double K = teiTensor(i, k, j, l);
                        G(i, j) += Density(k, l) * (J - 0.5 * K);
                    }

        return G;
    }

    MatrixXd HartreeFock::ConstructDensityRHF(const MatrixXd &Coeffs,
                                              const Eigen::Index Nelecs) {
        MatrixXd Density = MatrixXd::Zero(Coeffs.rows(), Coeffs.rows());
        for (int i = 0; i < Density.rows(); i++) {
            for (int j = 0; j < Density.rows(); j++) {
                for (int k = 0; k < (Nelecs / 2); k++) {
                    Density(i, j) += 2.0 * Coeffs(i, k) * Coeffs(j, k);
                }
            }
        }

        return Density;
    }

    double HartreeFock::CalculateDelta(const MatrixXd &oldDensity,
                                       const MatrixXd &density) {
        double delta = 0.0;
        for (int i = 0; i < density.rows(); i++) {
            for (int j = 0; j < density.rows(); j++) {
                delta += std::pow(density(i, j) - oldDensity(i, j), 2);
            }
        }

        return std::sqrt(delta);
    }

    double HartreeFock::CalculateEnergyRHF(const MatrixXd &Density,
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

    double HartreeFock::MP2ClosedShell(const RHFWavefunction &result, const GEng::Thread *thread,
                                       std::stringstream *console) {
        //<ab | ij > [MO Basis] = <mu kappa | nu lambda> [AO Basis]
        const TEITensor moTEI = Integrals::Transforms::TransformTEI(result.Integrals.eri, result.C);
        return MP2ClosedShell(result, moTEI, thread, console);
    }

    double HartreeFock::MP2ClosedShell(const RHFWavefunction &result, const TEITensor &moTEI,
                                       const GEng::Thread *thread,
                                       std::stringstream *console) {
        const VectorXi &occ = result.Occupation;
        const VectorXd &energies = result.Energies;

        double MP2_E = 0.0;

        //Number of doubly-occupied molecular orbitals
        Eigen::Index nOcc = std::distance(occ.begin(), std::find(occ.begin(), occ.end(), 0));

        if (console) {
            *console << std::setprecision(12);
            *console << "--- Moller Plesset Order 2 Correction ---\n";
        }

        //i and j run over the occupied orbitals
        for (Eigen::Index i = 0; i < nOcc; i++)
            for (Eigen::Index j = 0; j < nOcc; j++) {
                if (thread)
                    if (thread->ShouldStop())
                        return 0.0;

                //a and b run over the virtual orbitals
                for (Eigen::Index a = nOcc; a < occ.size(); a++)
                    for (Eigen::Index b = nOcc; b < occ.size(); b++) {
                        //Note that integrals are in chemist's notation
                        double intIJAB = moTEI(i, a, j, b);
                        double intIJBA = moTEI(i, b, j, a);
                        double epsilonIJAB = energies[a] + energies[b] - energies[i] - energies[j];

                        //MP2 formula below already integrates spin out for a closed shell
                        //configuration
                        MP2_E += intIJAB * (2.0 * intIJAB - intIJBA) / epsilonIJAB;
                    }
            }

        MP2_E = -MP2_E;

        if (console) {
            *console << "MP2 Correction Energy: " << MP2_E << '\n';
            *console << "E(Hartree Fock) + E(MP2) = " << result.Energy + MP2_E << " hartrees\n\n\n";
        }

        return MP2_E;
    }
}
