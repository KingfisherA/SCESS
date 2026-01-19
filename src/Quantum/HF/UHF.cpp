#include<Quantum/Q_PCH.h>
//
#include <chrono>
#include <Eigen/Core>
#include <Dahlia/Tools/Threading.h>
//
#include "HF.h"
#include "Quantum/Integrals/Integrals.h"
#include "Quantum/Basis/BasisSet.h"
#include "Quantum/Molecule.h"
#include "Utilities/Utilities.h"

//
//
namespace SCESS::HF {
    using namespace Eigen;
    //
    //
    void ConstructGUHF(const MatrixXd &D_Alpha,
                       const MatrixXd &D_Beta,
                       const TEITensor &eri,
                       MatrixXd &G_Alpha,
                       MatrixXd &G_Beta);

    //
    //
    HF_Wfn UHF(const Molecule &molecule,
              const BasisSet &basis,
              double convergeLimit,
              size_t scfMaxIterations,
              Integrals::MolInts *mintsOut,
              const DH::Thread *thread,
              std::stringstream *console) {
        using namespace std::chrono;

        if (console != nullptr)
            *console << "Starting integral evaluation for unrestricted Hartree-Fock\n";


        auto integralStartTP = high_resolution_clock::now();
        Integrals::MolInts mints = Integrals::EvaluateMolecularIntegrals(basis, molecule.GetAtoms(), thread);
        auto integralEndTP = high_resolution_clock::now();

        if (thread)
            if (thread->ShouldStop())
                return {};

        if (console)
            *console << "Integral evaluation completed in " <<
                    duration_cast<milliseconds>(integralEndTP - integralStartTP).count() << " milliseconds\n";

        auto alphaElectrons = molecule.GetNoElectrons() / 2;

        HF_Wfn wfn = UHF(mints,
                        alphaElectrons,
                        molecule.GetNoElectrons() - alphaElectrons,
                        molecule.CalculateVnn(),
                        convergeLimit,
                        scfMaxIterations,
                        thread,
                        console);

        if (mintsOut)
            (*mintsOut) = std::move(mints);

        return wfn;
    }
    //
    HF_Wfn UHF(const Integrals::MolInts &mints,
              uint32_t alphaElectrons,
              uint32_t betaElectrons,
              double Vnn,
              double convergeLimit,
              size_t scfMaxIterations,
              const DH::Thread *thread,
              std::stringstream *console) {
        if (console != nullptr)
            *console << "---Unrestricted Hartree Fock---\nUsing " << mints.S.rows() <<
                    " basis set functions\nBeginning integral evaluation...\n";

        const MatrixXd &S = mints.S;
        const MatrixXd &T = mints.T;
        const MatrixXd &V = mints.V;
        const TEITensor &tei = mints.eri; //Two-electron tensor

        const Index dimensions = S.rows();

        double delta = 10.0;
        Index scfCount = 0;
        double E_UHF = 0.0;

        MatrixXd h = T + V;

        //Canonical orthogonalisation transformation matrix
        SelfAdjointEigenSolver<MatrixXd> S_eig(S);
        MatrixXd s = S_eig.eigenvalues().real().cwiseInverse().cwiseSqrt().asDiagonal();
        MatrixXd U = S_eig.eigenvectors().real();
        MatrixXd X = U * s;

        MatrixXd D_Alpha = MatrixXd::Zero(dimensions, dimensions);
        MatrixXd D_Beta = MatrixXd::Zero(dimensions, dimensions);
        MatrixXd D_AlphaOld, D_BetaOld;
        MatrixXd F_Alpha, F_Beta;
        MatrixXd G_Alpha = MatrixXd::Zero(dimensions, dimensions);
        MatrixXd G_Beta = MatrixXd::Zero(dimensions, dimensions);
        MatrixXd C_Alpha, C_Beta;
        VectorXd Ep_Alpha, Ep_Beta;

        while (delta > convergeLimit && scfCount < scfMaxIterations) {
            if (thread)
                if (thread->ShouldStop())
                    return {};

            ConstructGUHF(D_Alpha, D_Beta, tei, G_Alpha, G_Beta);
            F_Alpha = h + G_Alpha;
            F_Alpha = (F_Alpha + F_Alpha.transpose()) * 0.5;
            F_Beta = h + G_Beta;
            F_Beta = (F_Beta + F_Beta.transpose()) * 0.5;

            MatrixXd F_AlphaPrime = X.transpose() * F_Alpha * X;
            MatrixXd F_BetaPrime = X.transpose() * F_Beta * X;

            SelfAdjointEigenSolver<MatrixXd> alphaEigen(F_AlphaPrime);
            SelfAdjointEigenSolver<MatrixXd> betaEigen(F_BetaPrime);

            Ep_Alpha = alphaEigen.eigenvalues().real();
            Ep_Beta = betaEigen.eigenvalues().real();

            C_Alpha = X * alphaEigen.eigenvectors().real();
            C_Beta = X * betaEigen.eigenvectors().real();

            D_AlphaOld = D_Alpha;
            D_BetaOld = D_Beta;

            D_Alpha = MatrixXd::Zero(h.rows(), h.rows());
            D_Beta = MatrixXd::Zero(h.rows(), h.rows());
            for (Index i = 0; i < h.rows(); ++i) {
                for (Index j = 0; j < h.cols(); ++j) {
                    Index k = 0;
                    for (k = 0; k < alphaElectrons; ++k) {
                        D_Alpha(i, j) += C_Alpha(i, k) * C_Alpha(j, k);
                    }

                    for (k = 0; k < betaElectrons; ++k) {
                        D_Beta(i, j) += C_Beta(i, k) * C_Beta(j, k);
                    }
                }
            }

            MatrixXd D_Delta = D_Alpha - D_AlphaOld;
            delta = D_Delta.cwiseProduct(D_Delta).sum();
            D_Delta = D_Beta - D_BetaOld;
            delta += D_Delta.cwiseProduct(D_Delta).sum();
            delta = sqrt(delta);

            E_UHF = CalculateEnergyUHF(D_Alpha,
                                       D_Beta,
                                       h,
                                       F_Alpha,
                                       F_Beta) + Vnn;

            if (console)
                *console << "   UHF(" << scfCount << "): E = " << E_UHF << "      dE = " << delta << '\n';

            scfCount++;
        }

        E_UHF = CalculateEnergyUHF(D_Alpha,
                                   D_Beta,
                                   h,
                                   F_Alpha,
                                   F_Beta) + Vnn;

        if (console) {
            if (scfCount < scfMaxIterations)
                *console << "---- UHF procedure converged! ----\nRHF converged after " <<
                        scfCount << " iterations.\nE(UHF) = " << E_UHF << " hartrees\n\n\n";
            else
                *console << "!!! UHF procedure failed to converge !!!\n\n\n";
        }

        VectorXi Occ_Alpha = VectorXi::Zero(h.rows());
        VectorXi Occ_Beta = VectorXi::Zero(h.rows());
        for (long i = long(alphaElectrons) - 1; i >= 0; i--) {
            Occ_Alpha(i) = 1;
        }
        for (long i = long(betaElectrons) - 1; i >= 0; i--) {
            Occ_Beta(i) = 1;
        }


        return HF_Wfn{
            E_UHF,
            F_Alpha,
            F_Beta,
            D_Alpha,
            D_Beta,
            C_Alpha,
            C_Beta,
            Ep_Alpha,
            Ep_Beta,
            Occ_Alpha,
            Occ_Beta,
            true
        };
    }
    //
    void ConstructGUHF(const MatrixXd &D_Alpha,
                       const MatrixXd &D_Beta,
                       const TEITensor &eri,
                       MatrixXd &G_Alpha,
                       MatrixXd &G_Beta) {
        G_Alpha.setZero();
        G_Beta.setZero();

        for (Index i = 0; i < D_Alpha.rows(); i++)
            for (Index j = 0; j < D_Alpha.rows(); j++)
                for (Index k = 0; k < D_Alpha.rows(); k++)
                    for (Index l = 0; l < D_Alpha.rows(); l++) {
                        double J = eri(i, j, k, l);
                        double K = eri(i, l, k, j);
                        G_Alpha(i, j) += D_Alpha(k, l) * (J - K) + D_Beta(k, l) * J;
                        G_Beta(i, j) += D_Beta(k, l) * (J - K) + D_Alpha(k, l) * J;
                    }
    }
    //
    double CalculateEnergyUHF(
        const MatrixXd &DensityAlpha,
        const MatrixXd &DensityBeta,
        const MatrixXd &Hcore,
        const MatrixXd &FockAlpha,
        const MatrixXd &FockBeta) {
        double energy = 0.0;
        for (Index i = 0; i < Hcore.rows(); i++) {
            for (Index j = 0; j < Hcore.cols(); j++) {
                energy += (DensityAlpha(i, j) + DensityBeta(i, j)) * Hcore(i, j);
                energy += DensityAlpha(i, j) * FockAlpha(i, j) + DensityBeta(i, j) * FockBeta(i, j);
            }
        }
        energy *= 0.5;
        return energy;
    }
}
