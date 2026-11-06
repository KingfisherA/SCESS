#include<Quantum/Q_PCH.h>
//
#include <chrono>
#include <Eigen/Core>
#include <Dahlia/Tools/Threading.h>
//
#include "HF.h"
#include "DIIS.h"
#include "Quantum/Integrals/Integrals.h"
#include "Quantum/Basis/BasisSet.h"
#include "Quantum/Molecule.h"
#include "Utilities/Utilities.h"
#include "Quantum/Integrals/IntegralTransforms.h"

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

        auto BuildFock = [&]() {
            ConstructGUHF(D_Alpha, D_Beta, tei, G_Alpha, G_Beta);
            F_Alpha = h + G_Alpha;
            F_Alpha = (F_Alpha + F_Alpha.transpose()) * 0.5;
            F_Beta = h + G_Beta;
            F_Beta = (F_Beta + F_Beta.transpose()) * 0.5;
        };

        auto DiagonaliseFock = [&](const MatrixXd &fockAlpha, const MatrixXd &fockBeta) {
            MatrixXd F_AlphaPrime = X.transpose() * fockAlpha * X;
            MatrixXd F_BetaPrime = X.transpose() * fockBeta * X;

            SelfAdjointEigenSolver<MatrixXd> alphaEigen(F_AlphaPrime);
            SelfAdjointEigenSolver<MatrixXd> betaEigen(F_BetaPrime);

            Ep_Alpha = alphaEigen.eigenvalues().real();
            Ep_Beta = betaEigen.eigenvalues().real();

            C_Alpha = X * alphaEigen.eigenvectors().real();
            C_Beta = X * betaEigen.eigenvectors().real();

            D_Alpha = C_Alpha.leftCols(alphaElectrons) * C_Alpha.leftCols(alphaElectrons).transpose();
            D_Beta = C_Beta.leftCols(betaElectrons) * C_Beta.leftCols(betaElectrons).transpose();
        };

        FockDIIS diis(8);

        while (delta > convergeLimit && scfCount < scfMaxIterations) {
            if (thread)
                if (thread->ShouldStop())
                    return {};

            BuildFock();
            E_UHF = CalculateEnergyUHF(D_Alpha, D_Beta, h, F_Alpha, F_Beta) + Vnn;

            D_AlphaOld = D_Alpha;
            D_BetaOld = D_Beta;

            //Pulay DIIS on both spins from the second iteration (the first density is the zero guess)
            if (scfCount > 0) {
                diis.Add({F_Alpha, F_Beta},
                         {
                             X.transpose() * (F_Alpha * D_Alpha * S - S * D_Alpha * F_Alpha) * X,
                             X.transpose() * (F_Beta * D_Beta * S - S * D_Beta * F_Beta) * X
                         });
            }
            if (diis.Size() >= 2) {
                const std::vector<MatrixXd> extrapolated = diis.Extrapolate();
                DiagonaliseFock(extrapolated[0], extrapolated[1]);
            } else
                DiagonaliseFock(F_Alpha, F_Beta);

            MatrixXd D_Delta = D_Alpha - D_AlphaOld;
            delta = D_Delta.cwiseProduct(D_Delta).sum();
            D_Delta = D_Beta - D_BetaOld;
            delta += D_Delta.cwiseProduct(D_Delta).sum();
            delta = sqrt(delta);

            if (console)
                *console << "   UHF(" << scfCount << "): E = " << E_UHF << "      dD = " << delta << '\n';

            scfCount++;
        }
        const bool converged = delta <= convergeLimit;

        //Final Fock build so that the energy pairs each density with the Fock matrix built from it, then take the
        //canonical orbitals of that Fock matrix
        BuildFock();
        E_UHF = CalculateEnergyUHF(D_Alpha,
                                   D_Beta,
                                   h,
                                   F_Alpha,
                                   F_Beta) + Vnn;
        DiagonaliseFock(F_Alpha, F_Beta);

        if (console) {
            if (converged)
                *console << "---- UHF procedure converged! ----\nUHF converged after " <<
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
        using namespace Integrals::Transforms;
        const MatrixXd J = CoulombMatrix(eri, D_Alpha + D_Beta);
        G_Alpha = J - ExchangeMatrix(eri, D_Alpha);
        G_Beta = J - ExchangeMatrix(eri, D_Beta);
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
