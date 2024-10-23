#include <pch.h>
#include<iostream>

#include "Chemistry/Quantum/DFT/DFT.h"
#include "Chemistry/Quantum/Integrals/Integrals.h"
#include "Chemistry/Quantum/Integrals/Grid/MolecularGrid.h"

#include "xc.h"
#include "XCFunctionals.h"
#include "Chemistry/Quantum/Integrals/Grid/Quadrature.h"

//
//
namespace SCESS {
    using namespace Eigen;

    double DensityFunctional::DFT(const Molecule &mol) {
        std::cout << "-- Density Functional Theory with B3LYP functional --\nEvaluating Molecular Integrals ....\n";

        constexpr double ALPHA = 0.9;
        const ConfigIntegrals &mints = Integrals::CalculateConfiguration(mol);

        const MatrixXd &S = mints.S;
        const MatrixXd &V = mints.V;
        const MatrixXd &T = mints.T;
        const TEITensor &tei = mints.tei;

        auto nElectrons = static_cast<uint32_t>(mol.GetNoElectrons());

        SelfAdjointEigenSolver<MatrixXd> eigensolver;
        eigensolver.compute(S);
        const VectorXd &SEigenValues = eigensolver.eigenvalues().real();
        const MatrixXd &SEigenVects = eigensolver.eigenvectors().real();
        MatrixXd SValsInvSqrt = SEigenValues.cwiseInverse().cwiseSqrt().asDiagonal();
        //Matrix to transform from AO to orthogonal basis set
        MatrixXd U = SEigenVects * SValsInvSqrt * SEigenVects.transpose();

        MolecularGrid molecularGrid(&mol);
        //Coarse
        //molecularGrid.Create(10, Quadrature::LEBEDEV_50);
        molecularGrid.Create(15, Quadrature::LEBEDEV_110);
        //Ultrafine
        //molecularGrid.Create(30, Quadrature::LEBEDEV_194);

        //Hamiltonian matrix
        MatrixXd H = T + V;
        // The Kohn-Sham Fock matrix
        MatrixXd F;
        // The one particle density matrices
        MatrixXd Density = MatrixXd::Zero(S.rows(), S.cols());
        //The J coulomb matrix
        MatrixXd J = MatrixXd::Zero(S.rows(), S.cols());
        //The exchange-correlation matrix
        MatrixXd ExchangeCorr = MatrixXd::Zero(S.rows(), S.cols());
        // The wavefunction coefficient matrix
        MatrixXd C;
        // The Kohn-Sham orbital energies
        VectorXd Energies;

        //The total electronic energy of the system
        double E = 0.0, lastE = 0.0;
        //Nuclear-nuclear repulsion energy
        double Vnn = mol.CalculateNuclearRepulsion();

        constexpr double E_CONVERGENCE = 1e-8;
        constexpr uint32_t MAX_ITERATIONS = 500;

        for (uint32_t iter = 0; iter < MAX_ITERATIONS; ++iter) {
            F = H + 2.0 * J + ExchangeCorr;
            MatrixXd Fprime = U.transpose() * F * U;

            SelfAdjointEigenSolver<MatrixXd> eigensolver2(Fprime);
            Energies = eigensolver2.eigenvalues().real();
            C = U * eigensolver2.eigenvectors().real();

            MatrixXd NewDensity = CreateDensityMatrix(C, nElectrons);

            if (iter == 0)
                Density = NewDensity;
            else
                Density = ALPHA * Density + (1 - ALPHA) * NewDensity;

            molecularGrid.SetDensity(Density);
            molecularGrid.CorrectDensities();

            J = CreateColoumbMatrixTEI(tei, Density);
            double E_XC = 0.0;
            CalculateXCMatrix(mol, molecularGrid, ExchangeCorr, E_XC);

            double kineticPotentialE = 2.0 * (Density * H).trace();
            double electronRepulsion = 2.0 * (Density * J).trace();

            E = kineticPotentialE + electronRepulsion + Vnn + E_XC;

            double dE = abs(E - lastE);
            std::cout << "DFT(" << iter << ") E(DFT) = " << E << " dE = " << dE << '\n';

            if (dE < E_CONVERGENCE) {
                std::cout << "--- DFT procedure has converged! ---\n";
                break;
            }
            lastE = E;
        }

        std::cout << "--- Final DFT Energy ---\n";
        std::cout << "E(DFT) = " << E << " hartrees\n";

        return E;
    }

    MatrixXd DensityFunctional::CreateColoumbMatrixTEI(const TEITensor &tei,
                                                       const MatrixXd &Density) {
        MatrixXd J = MatrixXd::Zero(Density.rows(), Density.cols());

        for (Eigen::Index i = 0; i < J.rows(); i++)
            for (Eigen::Index j = 0; j < J.rows(); j++)
                for (Eigen::Index k = 0; k < J.rows(); k++)
                    for (Eigen::Index l = 0; l < J.rows(); l++) {
                        //Coulomb integral; <ik | jl> = (ij | kl)
                        J(i, j) += Density(k, l) * tei(i, j, k, l);
                    }

        return J;
    }

    MatrixXd DensityFunctional::CreateDensityMatrix(const MatrixXd &C, uint32_t Nelecs) {
        MatrixXd Density = MatrixXd::Zero(C.rows(), C.rows());
        for (int i = 0; i < Density.rows(); i++) {
            for (int j = 0; j < Density.rows(); j++) {
                for (int k = 0; k < (Nelecs / 2); k++) {
                    Density(i, j) += C(i, k) * C(j, k);
                }
            }
        }

        return Density;
    }

    void DensityFunctional::CalculateXCMatrix(const Molecule &mol, const MolecularGrid &grid,
                                              MatrixXd &ExchangeCorr, double &E_XC) {
        for (Eigen::Index i = 0; i < ExchangeCorr.size(); i++)
            ExchangeCorr.data()[i] = 0.0;

        VectorXd densities = grid.GetDensities();
        VectorXd weights = grid.GetWeights();
        MatrixXd amplitudes = grid.GetAmplitudes();

        VectorXd ex; // exchange energy
        VectorXd vxa; // exchange potential (spin-up)
        VectorXd vxb; // exchange potential (spin-down)

        VectorXd ec; // correlation energy
        VectorXd vca; // correlation potential (spin-up)
        VectorXd vcb; // correlation potential (spin-down)

        VectorXd densitiesa = densities * 0.5;
        VectorXd densitiesb = densities * 0.5;

        // calculate exchange energies and potential
        XCFunctionals::xalpha_x_functional(densitiesa, densitiesb, ex, vxa, vxb);

        // calculate correlation energy and potential
        XCFunctionals::vwm_c_functional(densitiesa, densitiesb, ec, vca, vcb);

        // calculate total exchange-correlation energy
        E_XC = weights.dot(ex + ec);

        // calculate gridpoint-wise xc potential

        // I still don't understand why I need to put a 0.5 here... (see also the function above)
        // Does it have anything to do how I construct the density matrix P?
        // Can someone who knows the answer send me an e-mail?
        VectorXd wva = weights.cwiseProduct((vxa + vxb + vca + vcb) * 0.5);

        for (uint32_t i = 0; i < mol.GetNoOrbitals(); i++) {
            VectorXd row = amplitudes.row(i);
            VectorXd wva_i = wva.cwiseProduct(row);
            for (uint32_t j = 0; j < mol.GetNoOrbitals(); j++) {
                ExchangeCorr(i, j) = wva_i.dot(amplitudes.row(j));
            }
        }
    }
}
