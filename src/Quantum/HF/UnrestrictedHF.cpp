#include<pch.h>

#include <chrono>
#include <Eigen/Core>
#include <GraphicsEngine/Utilities/Threading.h>

#include "HartreeFock.h"
#include "Quantum/Integrals/Integrals.h"
#include "Quantum/Basis/BasisSet.h"
#include "Quantum/Molecule.h"
#include "Utilities/Utilities.h"

namespace SCESS {
    using namespace Eigen;

    uRef<UHFWavefunction> HartreeFock::UHF(const Molecule &molecule,
                                           const BasisSet &basis,
                                           double convergeLimit,
                                           size_t scfMaxIterations,
                                           const GEng::Thread *thread,
                                           std::stringstream *console,
                                           std::string_view basisSetName) {
        using namespace std::chrono;

        if (console != nullptr)
            *console << "Starting integral evaluation for unrestricted Hartree-Fock\n";


        auto integralStartTP = high_resolution_clock::now();
        Integrals::ConfigInts integrals = Integrals::CalculateConfiguration(basis, molecule.GetAtoms(), thread);
        auto integralEndTP = high_resolution_clock::now();

        if (thread)
            if (thread->ShouldStop())
                return nullptr;

        if (console)
            *console << "Integral evaluation completed in " <<
                    duration_cast<milliseconds>(integralEndTP - integralStartTP).count() << " milliseconds\n";

        auto alphaElectrons = molecule.GetNoElectrons() / 2;

        return UHF(integrals, alphaElectrons, molecule.GetNoElectrons() - alphaElectrons,
                   molecule.CalculateNuclearRepulsion(), convergeLimit, scfMaxIterations, thread, console,
                   basisSetName);
    }

    uRef<UHFWavefunction> HartreeFock::UHF(const Integrals::ConfigInts &mints,
                                           uint32_t alphaElectrons,
                                           uint32_t betaElectrons,
                                           double Vnn,
                                           double convergeLimit,
                                           size_t scfMaxIterations,
                                           const GEng::Thread *thread,
                                           std::stringstream *console,
                                           std::string_view basisSet) {
        if (console != nullptr)
            *console << "---Unrestricted Hartree Fock---\nUsing " << basisSet << " basis set with " << mints.S.rows() <<
                    " basis set functions\nBeginning integral evaluation...\n";

        const MatrixXd &S = mints.S;
        const MatrixXd &T = mints.T;
        const MatrixXd &V = mints.V;
        const TEITensor &tei = mints.eri; //Two-electron tensor

        const Index dimensions = S.rows();

        double delta = 10.0;
        Index scfCount = 0;
        double currentE = 0.0;

        const MatrixXd &Hcore = T + V;

        SelfAdjointEigenSolver<MatrixXd> eigensolver;
        eigensolver.compute(S);
        const VectorXd &SEigenValues = eigensolver.eigenvalues().real();
        const MatrixXd &SEigenVects = eigensolver.eigenvectors().real();
        MatrixXd SValsInvSqrt = SEigenValues.cwiseInverse().cwiseSqrt().asDiagonal();
        MatrixXd TransformMat = SEigenVects * SValsInvSqrt * SEigenVects.transpose();
        TransformMat = (TransformMat + TransformMat.transpose()) * 0.5;

        MatrixXd DensityAlpha = MatrixXd::Zero(dimensions, dimensions);
        MatrixXd DensityBeta = MatrixXd::Zero(dimensions, dimensions);
        MatrixXd OldDensityAlpha, OldDensityBeta;
        MatrixXd FockAlpha, FockBeta;
        MatrixXd CoeffsAlpha, CoeffsBeta;
        VectorXd EnergiesAlpha, EnergiesBeta;

        while (delta > convergeLimit && scfCount < scfMaxIterations) {
            if (thread)
                if (thread->ShouldStop())
                    return nullptr;

            const auto [GAlpha, GBeta] = ConstructGUHF(DensityAlpha, DensityBeta, tei);
            FockAlpha = Hcore + GAlpha;
            FockAlpha = (FockAlpha + FockAlpha.transpose()) * 0.5;
            FockBeta = Hcore + GBeta;
            FockBeta = (FockBeta + FockBeta.transpose()) * 0.5;

            MatrixXd FockAlphaPrime = TransformMat.transpose() * FockAlpha * TransformMat;
            MatrixXd FockBetaPrime = TransformMat.transpose() * FockBeta * TransformMat;

            SelfAdjointEigenSolver<MatrixXd> alphaEigen(FockAlphaPrime);
            SelfAdjointEigenSolver<MatrixXd> betaEigen(FockBetaPrime);

            EnergiesAlpha = alphaEigen.eigenvalues().real();
            EnergiesBeta = betaEigen.eigenvalues().real();

            CoeffsAlpha = TransformMat * alphaEigen.eigenvectors().real();
            CoeffsBeta = TransformMat * betaEigen.eigenvectors().real();

            OldDensityAlpha = DensityAlpha;
            OldDensityBeta = DensityBeta;

            DensityAlpha = MatrixXd::Zero(Hcore.rows(), Hcore.rows());
            DensityBeta = MatrixXd::Zero(Hcore.rows(), Hcore.rows());
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

            MatrixXd DensityMatriciesDelta = DensityAlpha - OldDensityAlpha;
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

        VectorXi OccupationAlpha = VectorXi::Zero(Hcore.rows());
        VectorXi OccupationBeta = VectorXi::Zero(Hcore.rows());
        for (long i = long(alphaElectrons) - 1; i >= 0; i--) {
            OccupationAlpha(i) = 1;
        }
        for (long i = long(betaElectrons) - 1; i >= 0; i--) {
            OccupationBeta(i) = 1;
        }


        return std::make_unique<UHFWavefunction>(
            FockAlpha, FockBeta, Hcore, DensityAlpha, DensityBeta,
            CoeffsAlpha, CoeffsBeta, EnergiesAlpha, EnergiesBeta, OccupationAlpha, OccupationBeta,
            mints, currentE
        );
    }

    std::pair<MatrixXd, MatrixXd> HartreeFock::ConstructGUHF(const MatrixXd &DensityAlpha,
                                                             const MatrixXd &DensityBeta,
                                                             const TEITensor &teiTensor) {
        MatrixXd GAlpha = MatrixXd::Zero(DensityAlpha.rows(), DensityAlpha.rows());
        MatrixXd GBeta = MatrixXd::Zero(DensityAlpha.rows(), DensityAlpha.rows());

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

    double HartreeFock::CalculateEnergyUHF(
        const MatrixXd &DensityAlpha,
        const MatrixXd &DensityBeta,
        const MatrixXd &Hcore,
        const MatrixXd &FockAlpha,
        const MatrixXd &FockBeta) {
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
}
