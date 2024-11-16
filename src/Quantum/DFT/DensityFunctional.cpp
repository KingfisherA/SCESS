#include <pch.h>
#include<iostream>

#include  "Quantum/DFT/DensityFunctional.h"

#include <GraphicsEngine/Utilities/Threading.h>

#include  "Quantum/Integrals/Integrals.h"
#include  "Quantum/Integrals/Grid/MolecularGrid.h"

#include "xc.h"
#include "XCFunctionals.h"
#include  "Quantum/Integrals/Grid/Quadrature.h"

//
//
namespace SCESS {
    using namespace Eigen;

    uRef<KSWavefunction> DensityFunctional::KohnShamDFT(const Molecule &mol, MolecularGrid::Resolution gridRes,
                                                        const GEng::Thread *thread,
                                                        std::stringstream *console, std::string_view basisSetName) {
        if (console) {
            *console << "--- Kohn-Sham Density Functional Theory ---\n";
            *console << "--- Selected Functionals ---\n";
            *console << "   - Exchange = Slater-X-LDA\n";
            *console << "   - Correlation = VWM-C-LDA\n\n";
            *console << "Evaluating Molecular Integrals ....\n";
        }

        constexpr double ALPHA = 0.9;
        const ConfigIntegrals &mints = Integrals::CalculateConfiguration(mol, thread);

        if(thread)
            if(thread->ShouldStop())
                return nullptr;

        if (console)
            *console << "Molecular Integrals Evaluated!\n";

        const MatrixXd &S = mints.S;
        const MatrixXd &V = mints.V;
        const MatrixXd &T = mints.T;
        const TEITensor &tei = mints.eri;

        auto nElectrons = static_cast<uint32_t>(mol.GetNoElectrons());

        SelfAdjointEigenSolver<MatrixXd> eigensolver;
        eigensolver.compute(S);
        const VectorXd &SEigenValues = eigensolver.eigenvalues().real();
        const MatrixXd &SEigenVects = eigensolver.eigenvectors().real();
        MatrixXd SValsInvSqrt = SEigenValues.cwiseInverse().cwiseSqrt().asDiagonal();
        //Matrix to transform from AO to orthogonal basis set
        MatrixXd U = SEigenVects * SValsInvSqrt * SEigenVects.transpose();

        if (console)
            *console << "Creating Molecular Gauss-Chebyshev-Lebedev grid...\n";
        MolecularGrid molecularGrid(&mol);
        molecularGrid.Create(gridRes);
        if (console) {
            *console << "Molecular Gauss-Chebyshev-Lebedev grid created with:\n";
            *console << "   - Radial Points: " << gridRes.radialPoints << '\n';
            *console << "   - Lebedev Points: " << Quadrature::nLebedevPoints[(int) gridRes.lebedevOrder] << '\n';
            *console << "   - Atomic Grids: " << molecularGrid.GetAtomicGrids().size() << '\n';
            *console << "   - Total Number of Points: " << molecularGrid.GetGridSize() << '\n';
        }


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

        constexpr double E_CONVERGENCE = 1e-5;
        constexpr uint32_t MAX_ITERATIONS = 500;

        if (console)
            *console << "--- Starting SCF cycles ---\n";

        for (uint32_t iter = 0; iter < MAX_ITERATIONS; ++iter) {

            //Poll to check whether we should cancel the calculation
            if(thread)
                if(thread->ShouldStop())
                    return nullptr;

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
            //Exchange-correlation energy
            double E_XC = 0.0;
            CalculateXCMatrix(mol, molecularGrid, ExchangeCorr, E_XC);

            //Single electron energy
            double kineticPotentialE = 2.0 * (Density * H).trace();
            //Coulomb repulsion energy
            double electronRepulsion = 2.0 * (Density * J).trace();

            E = kineticPotentialE + electronRepulsion + Vnn + E_XC;

            double dE = abs(E - lastE);
            if (console) *console << "   DFT(" << iter << "): E(DFT) = " << E << "          dE = " << dE << '\n';

            if (dE < E_CONVERGENCE) {
                if (console) *console << "--- DFT procedure has converged! ---\n";
                break;
            }
            lastE = E;
        }

        if (console) {
            *console << "--- Final DFT Energy ---\n";
            *console << "  E(DFT) = " << E << " hartrees\n";
        }

        VectorXi Occupation = VectorXi::Zero(H.rows());

        for (int64_t i = ((int64_t) mol.GetNoElectrons() / 2) - 1; i >= 0; --i)
            Occupation(i) = 2;

        return std::make_unique<KSWavefunction>(F, C, Density, E, Energies, Occupation, mints);
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
