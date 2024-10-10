#include<pch.h>

#include "Chemistry/Quantum/HF/HartreeFock.h"
#include"Chemistry/Quantum/Basis/Gaussian.h"
#include"Chemistry/Quantum/Integrals/Integrals.h"
#include"Chemistry/Quantum/Molecule.h"

#include<chrono>
#include<algorithm>
#include <Eigen/QR>

#include "Utilities/Utilities.h"


//
//
namespace SCESS {
    uRef<RHFResult> HartreeFock::Restricted(const Molecule &molecule,
                                            const Molecule *lastMolecule,
                                            const RHFResult *lastResult,
                                            const double convergeLimit,
                                            const size_t scfMaxIterations,
                                            std::stringstream *console,
                                            const std::string_view basisSet) {
        using namespace std::chrono;

        if (console)
            *console << "---Restricted Hartree Fock---\nUsing " << basisSet << " basis set with " << molecule.
                    GetOrbitals().size() << " basis set functions \nBeginning integral evaluation...\n";

        auto integralStartTP = high_resolution_clock::now();
        ConfigIntegrals integrals = Integrals::CalculateConfiguration(molecule);
        auto integralEndTP = high_resolution_clock::now();

        const Eigen::MatrixXd &S = integrals.S;
        const Eigen::MatrixXd &T = integrals.T;
        const Eigen::MatrixXd &V = integrals.V;
        const TEITensor &tei = integrals.tei; //Two-electron tensor

        if (console)
            *console << "Integral evaluation completed in " <<
                    duration_cast<milliseconds>(integralEndTP - integralStartTP).count() << " milliseconds\n";

        const Eigen::Index dimensions = S.rows();
        double delta = 1.0;
        double Vnn = molecule.CalculateNuclearRepulsion();
        double energy = 0.0, lastEnergy = DBL_MAX;

        Eigen::MatrixXd Hcore = T + V;

        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver;
        eigensolver.compute(S);
        const Eigen::VectorXd &SEigenValues = eigensolver.eigenvalues().real();
        const Eigen::MatrixXd &SEigenVects = eigensolver.eigenvectors().real();
        Eigen::MatrixXd SValsInvSqrt = SEigenValues.cwiseInverse().cwiseSqrt().asDiagonal();
        Eigen::MatrixXd TransformMat = SEigenVects * SValsInvSqrt * SEigenVects.transpose();
        TransformMat = (TransformMat + TransformMat.transpose()) * 0.5;

        Eigen::MatrixXd Fock, Coeffs, Density, OldDensity;
        Eigen::VectorXd Energies;

        if (lastResult && lastMolecule) {
            Eigen::MatrixXd S21 = Integrals::CalculateMixedOverlapMatrix(molecule, *lastMolecule);
            Eigen::MatrixXd S22Inv = S.completeOrthogonalDecomposition().pseudoInverse();

            Density = S22Inv * (S21 * lastResult->Density * S21.transpose()) * S22Inv;
            Density.conservativeResize(dimensions, dimensions);
        } else
            Density = Eigen::MatrixXd::Zero(dimensions, dimensions);

        if (console)
            *console << "Beginning SCF iterations...\n";


        std::vector<Eigen::MatrixXd> DIIS_Errors;
        DIIS_Errors.reserve(7);
        std::vector<Eigen::MatrixXd> FockMatrices;
        FockMatrices.reserve(7);

        Eigen::Index scfIter = 1;
        for (; scfIter < scfMaxIterations + 1; scfIter++) {
            //Generate the Fock matrix
            Fock = Hcore + ConstructGRHF(Density, tei);
            Fock = (Fock + Fock.transpose()) * 0.5;

            //Build the DIIS error vector (FDS-SDF)
            Eigen::MatrixXd DIIS_Error = Fock * Density * S - S * Density * Fock;
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

                Eigen::MatrixXd ErrorMat = Eigen::MatrixXd::Constant(DIIS_Size + 1, DIIS_Size + 1, -1.0);
                ErrorMat(ErrorMat.rows() - 1, ErrorMat.cols() - 1) = 0.0;

                double maxVal = 0.0;
                for (Eigen::Index i = 0; i < DIIS_Errors.size(); i++)
                    for (Eigen::Index j = i; j < DIIS_Errors.size(); j++) {
                        ErrorMat(i, j) = ErrorMat(j, i) = DIIS_Errors[i].cwiseProduct(DIIS_Errors[j]).sum();
                        maxVal = std::max(maxVal, abs(ErrorMat(i, j)));
                    }

                ErrorMat.block(0, 0, ErrorMat.rows() - 1, ErrorMat.cols() - 1) /= maxVal;

                Eigen::VectorXd Residual = Eigen::VectorXd::Zero(DIIS_Size + 1);
                Residual(Residual.rows() - 1) = -1.0;

                Eigen::VectorXd DIISCoeffs = ErrorMat.lu().solve(Residual);

                Fock = Eigen::MatrixXd::Zero(Fock.rows(), Fock.cols());
                for (Eigen::Index i = 0; i < FockMatrices.size(); i++)
                    Fock += DIISCoeffs(i) * FockMatrices[i];
            }

            //Transform the Fock matrix in the orthonormal basis set
            Eigen::MatrixXd FockPrime = TransformMat.transpose() * Fock * TransformMat;

            //Diagonalise the transformed Fock matrix
            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver2(FockPrime);
            //eigensolver2.compute(FockPrime);
            Energies = eigensolver2.eigenvalues().real();

            //Transform coefficients back into normal basis
            const Eigen::MatrixXd &CPrime = eigensolver2.eigenvectors().real();
            Coeffs = TransformMat * CPrime;

            //Recompute the density matrix
            OldDensity = Density;
            Density = ConstructDensityRHF(Coeffs, molecule.GetNoElectrons());

            //Calculate delta change to the density matrix
            // const Eigen::MatrixXd DensityMatriciesDelta = Density - OldDensity;
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

        Eigen::VectorXi Occupation = Eigen::VectorXi::Zero(Hcore.rows());

        for (long i = (molecule.GetNoElectrons() / 2) - 1; i >= 0; --i) {
            Occupation(i) = 2;
        }


        return std::make_unique<RHFResult>(Fock, Hcore, Density, Coeffs, Energies, Occupation, integrals, energy);
    }

    uRef<UHFResult> HartreeFock::Unrestricted(const Molecule &molecule, const double convergeLimit,
                                              const size_t scfMaxIterations,
                                              std::stringstream *console, const std::string_view basisSet) {
        using namespace std::chrono;

        if (console != nullptr)
            *console << "---Unrestricted Hartree Fock---\nUsing " << basisSet << " basis set with " << molecule.
                    GetOrbitals().size() << " basis set functions\nBeginning integral evaluation...\n";

        //Spin-up electrons
        Eigen::Index alphaElectrons = molecule.GetNoElectrons() / 2;
        //Spin-down electrons
        Eigen::Index betaElectrons = molecule.GetNoElectrons() - alphaElectrons;

        auto integralStartTP = high_resolution_clock::now();
        ConfigIntegrals integrals = Integrals::CalculateConfiguration(molecule);
        auto integralEndTP = high_resolution_clock::now();

        const Eigen::MatrixXd &S = integrals.S;
        const Eigen::MatrixXd &T = integrals.T;
        const Eigen::MatrixXd &V = integrals.V;
        const TEITensor &tei = integrals.tei; //Two-electron tensor

        if (console)
            *console << "Integral evaluation completed in " <<
                    duration_cast<milliseconds>(integralEndTP - integralStartTP).count() << " milliseconds\n";

        const Eigen::Index dimensions = S.rows();


        double delta = 1.0;
        Eigen::Index scfCount = 0;
        double Vnn = molecule.CalculateNuclearRepulsion();
        double currentE = 0.0;

        const Eigen::MatrixXd &Hcore = T + V;

        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver;
        eigensolver.compute(S);
        const Eigen::VectorXd &SEigenValues = eigensolver.eigenvalues().real();
        const Eigen::MatrixXd &SEigenVects = eigensolver.eigenvectors().real();
        Eigen::MatrixXd SValsInvSqrt = SEigenValues.cwiseInverse().cwiseSqrt().asDiagonal();
        Eigen::MatrixXd TransformMat = SEigenVects * SValsInvSqrt * SEigenVects.transpose();
        TransformMat = (TransformMat + TransformMat.transpose()) * 0.5;

        Eigen::MatrixXd DensityAlpha = Eigen::MatrixXd::Zero(dimensions, dimensions);
        Eigen::MatrixXd DensityBeta = Eigen::MatrixXd::Zero(dimensions, dimensions);
        Eigen::MatrixXd OldDensityAlpha, OldDensityBeta;
        Eigen::MatrixXd FockAlpha, FockBeta;
        Eigen::MatrixXd CoeffsAlpha, CoeffsBeta;
        Eigen::VectorXd EnergiesAlpha, EnergiesBeta;

        while (delta > convergeLimit && scfCount < scfMaxIterations) {
            const auto [GAlpha, GBeta] = ConstructGUHF(DensityAlpha, DensityBeta, tei);
            FockAlpha = Hcore + GAlpha;
            FockAlpha = (FockAlpha + FockAlpha.transpose()) * 0.5;
            FockBeta = Hcore + GBeta;
            FockBeta = (FockBeta + FockBeta.transpose()) * 0.5;

            Eigen::MatrixXd FockAlphaPrime = TransformMat.transpose() * FockAlpha * TransformMat;
            Eigen::MatrixXd FockBetaPrime = TransformMat.transpose() * FockBeta * TransformMat;

            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> alphaEigen(FockAlphaPrime);
            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> betaEigen(FockBetaPrime);

            EnergiesAlpha = alphaEigen.eigenvalues().real();
            EnergiesBeta = betaEigen.eigenvalues().real();

            CoeffsAlpha = TransformMat * alphaEigen.eigenvectors().real();
            CoeffsBeta = TransformMat * betaEigen.eigenvectors().real();

            OldDensityAlpha = DensityAlpha;
            OldDensityBeta = DensityBeta;

            DensityAlpha = Eigen::MatrixXd::Zero(Hcore.rows(), Hcore.rows());
            DensityBeta = Eigen::MatrixXd::Zero(Hcore.rows(), Hcore.rows());
            for (Eigen::Index i = 0; i < Hcore.rows(); ++i) {
                for (Eigen::Index j = 0; j < Hcore.cols(); ++j) {
                    Eigen::Index k = 0;
                    for (k = 0; k < alphaElectrons; ++k) {
                        DensityAlpha(i, j) += CoeffsAlpha(i, k) * CoeffsAlpha(j, k);
                    }

                    for (k = 0; k < betaElectrons; ++k) {
                        DensityBeta(i, j) += CoeffsBeta(i, k) * CoeffsBeta(j, k);
                    }
                }
            }

            Eigen::MatrixXd DensityMatriciesDelta = DensityAlpha - OldDensityAlpha;
            delta = DensityMatriciesDelta.cwiseProduct(DensityMatriciesDelta).sum();
            DensityMatriciesDelta = DensityBeta - OldDensityBeta;
            delta += DensityMatriciesDelta.cwiseProduct(DensityMatriciesDelta).sum();
            delta = sqrt(delta);

            currentE = CalculateEnergyUHF(DensityAlpha,
                                          DensityBeta,
                                          Hcore,
                                          FockAlpha,
                                          FockBeta) + Vnn;

            if (console)
                *console << "   UHF(" << scfCount << "): E = " << currentE << "      dE = " << delta << '\n';

            scfCount++;
        }

        currentE = CalculateEnergyUHF(DensityAlpha,
                                      DensityBeta,
                                      Hcore,
                                      FockAlpha,
                                      FockBeta) + Vnn;

        if (console) {
            if (scfCount < scfMaxIterations)
                *console << "---- UHF procedure converged! ----\nRHF converged after " <<
                        scfCount << " iterations.\nE(UHF) = " << currentE << " hartrees\n\n\n";
            else
                *console << "!!! UHF procedure failed to converge !!!\n\n\n";
        }

        Eigen::VectorXi OccupationAlpha = Eigen::VectorXi::Zero(Hcore.rows());
        Eigen::VectorXi OccupationBeta = Eigen::VectorXi::Zero(Hcore.rows());
        for (long i = alphaElectrons - 1; i >= 0; i--) {
            OccupationAlpha(i) = 1;
        }
        for (long i = betaElectrons - 1; i >= 0; i--) {
            OccupationBeta(i) = 1;
        }


        return std::make_unique<UHFResult>(
            FockAlpha, FockBeta, Hcore, DensityAlpha, DensityBeta,
            CoeffsAlpha, CoeffsBeta, EnergiesAlpha, EnergiesBeta, OccupationAlpha, OccupationBeta,
            integrals, currentE
        );
    }

    Eigen::MatrixXd HartreeFock::ConstructGRHF(const Eigen::MatrixXd &Density,
                                               const TEITensor &teiTensor) {
        Eigen::MatrixXd G = Eigen::MatrixXd::Zero(Density.rows(), Density.rows());

        for (Eigen::DenseIndex i = 0; i < G.rows(); i++) {
            for (Eigen::DenseIndex j = 0; j < G.rows(); j++) {
                for (Eigen::DenseIndex k = 0; k < G.rows(); k++) {
                    for (Eigen::DenseIndex l = 0; l < G.rows(); l++) {
                        //Coulomb integral
                        double J = teiTensor(i, j, k, l);
                        //Exchange integral
                        double K = teiTensor(i, l, k, j);
                        G(i, j) += Density(k, l) * (J - 0.5 * K);
                    }
                }
            }
        }

        return G;
    }

    std::pair<Eigen::MatrixXd, Eigen::MatrixXd> HartreeFock::ConstructGUHF(const Eigen::MatrixXd &DensityAlpha,
                                                                           const Eigen::MatrixXd &DensityBeta,
                                                                           const TEITensor &teiTensor) {
        Eigen::MatrixXd GAlpha = Eigen::MatrixXd::Zero(DensityAlpha.rows(), DensityAlpha.rows());
        Eigen::MatrixXd GBeta = Eigen::MatrixXd::Zero(DensityAlpha.rows(), DensityAlpha.rows());

        for (Eigen::Index i = 0; i < DensityAlpha.rows(); i++) {
            for (Eigen::Index j = 0; j < DensityAlpha.rows(); j++) {
                for (Eigen::Index k = 0; k < DensityAlpha.rows(); k++) {
                    for (Eigen::Index l = 0; l < DensityAlpha.rows(); l++) {
                        const double J = teiTensor(i, j, k, l);
                        const double K = teiTensor(i, l, k, j);
                        GAlpha(i, j) += DensityAlpha(k, l) * (J - K) + DensityBeta(k, l) * J;
                        GBeta(i, j) += DensityBeta(k, l) * (J - K) + DensityAlpha(k, l) * J;
                    }
                }
            }
        }
        return {GAlpha, GBeta};
    }

    Eigen::MatrixXd HartreeFock::ConstructDensityRHF(const Eigen::MatrixXd &Coeffs,
                                                     const Eigen::Index Nelecs) {
        Eigen::MatrixXd Density = Eigen::MatrixXd::Zero(Coeffs.rows(), Coeffs.rows());
        for (int i = 0; i < Density.rows(); i++) {
            for (int j = 0; j < Density.rows(); j++) {
                for (int k = 0; k < (Nelecs / 2); k++) {
                    Density(i, j) += 2.0 * Coeffs(i, k) * Coeffs(j, k);
                }
            }
        }
        return Density;
    }

    double HartreeFock::CalculateDelta(const Eigen::MatrixXd &oldDensity,
                                       const Eigen::MatrixXd &density) {
        double delta = 0.0;
        for (int i = 0; i < density.rows(); i++) {
            for (int j = 0; j < density.rows(); j++) {
                delta += std::pow(density(i, j) - oldDensity(i, j), 2);
            }
        }

        return std::sqrt(delta);
    }

    double HartreeFock::CalculateEnergyRHF(const Eigen::MatrixXd &Density,
                                           const Eigen::MatrixXd &Hcore, const Eigen::MatrixXd &Fock) {
        double result = 0.0;

        for (int i = 0; i < Fock.rows(); i++) {
            for (int j = 0; j < Fock.rows(); j++) {
                result += Density(i, j) * (Hcore(i, j) + Fock(i, j));
            }
        }

        return 0.5 * result;
    }

    double HartreeFock::CalculateEnergyUHF(
        const Eigen::MatrixXd &DensityAlpha,
        const Eigen::MatrixXd &DensityBeta,
        const Eigen::MatrixXd &Hcore,
        const Eigen::MatrixXd &FockAlpha,
        const Eigen::MatrixXd &FockBeta) {
        double energy = 0.0;
        for (Eigen::Index i = 0; i < Hcore.rows(); i++) {
            for (Eigen::Index j = 0; j < Hcore.cols(); j++) {
                energy += (DensityAlpha(i, j) + DensityBeta(i, j)) * Hcore(i, j);
                energy += DensityAlpha(i, j) * FockAlpha(i, j) + DensityBeta(i, j) * FockBeta(i, j);
            }
        }
        energy *= 0.5;
        return energy;
    }

    double HartreeFock::MP2ClosedShell(const RHFResult &result, std::stringstream *console) {
        //<ab | ij > [MO Basis] = <mu kappa | nu lambda> [AO Basis]
        const TEITensor moTEI = Integrals::Transforms::ToMolecularTEI(result.Integrals.tei, result.Coeffs);
        return MP2ClosedShell(result, moTEI, console);
    }

    double HartreeFock::MP2ClosedShell(const RHFResult &result, const TEITensor &moTEI,
                                       std::stringstream *console) {
        const Eigen::VectorXi &occ = result.Occupation;
        const Eigen::VectorXd &energies = result.Energies;

        double MP2_E = 0.0;

        //Number of doubly-occupied molecular orbitals
        Eigen::Index nOcc = std::distance(occ.begin(), std::find(occ.begin(), occ.end(), 0));

        if (console) {
            *console << std::setprecision(12);
            *console << "--- Moller Plesset Order 2 Correction ---\n";
        }

        //i and j run over the occupied orbitals
        for (Eigen::Index i = 0; i < nOcc; i++)
            for (Eigen::Index j = 0; j < nOcc; j++)
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

        MP2_E = -MP2_E;

        if (console) {
            *console << "MP2 Correction Energy: " << MP2_E << '\n';
            *console << "E(Hartree Fock) + E(MP2) = " << result.Energy + MP2_E << " hartrees\n\n\n";
        }

        return MP2_E;
    }
}
