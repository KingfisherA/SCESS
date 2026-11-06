#include <Quantum/Q_PCH.h>
#include  <Dahlia/Tools/Threading.h>
//
#include "Quantum/DFT/DFT.h"
#include "Quantum/Basis/BasisSet.h"
#include "Quantum/Integrals/Integrals.h"
#include "Quantum/Integrals/IntegralTransforms.h"
#include "Quantum/HF/DIIS.h"
#include "Quantum/Integrals/Grid/MolGrid.h"
#include "XC.h"
#include "Quantum/Integrals/Grid/Quadrature.h"

//
//q
namespace SCESS::DFT {
    using namespace Eigen;
    using namespace Quadrature;
    //
    //
    MatrixXd CreateDensityMatrix(const MatrixXd &C, uint32_t Nelecs);
    void CalculateXCMatrix(const MolGrid &grid,
                           MatrixXd &ExchangeCorr,
                           double &E_XC);
    //
    //
    DFT::KS_Wfn KS_DFT(const Molecule &mol,
                       const BasisSet &basis,
                       MolGrid::Resolution gridRes,
                       const DH::Thread *thread,
                       std::stringstream *console) {
        if (console) {
            *console << "--- Kohn-Sham Density Functional Theory ---\n";
            *console << "--- Selected Functionals ---\n";
            *console << "   - Exchange = Slater-X-LDA\n";
            *console << "   - Correlation = VWM-C-LDA\n\n";
            *console << "Evaluating Molecular Integrals ....\n";
        }

        const Integrals::MolInts &mints = Integrals::EvaluateMolecularIntegrals(basis, mol.GetAtoms(), thread);

        if (thread)
            if (thread->ShouldStop())
                return {};

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
        MolGrid MolGrid(&mol, &basis);
        MolGrid.Create(gridRes);
        if (console) {
            *console << "Molecular Gauss-Chebyshev-Lebedev grid created with:\n";
            *console << "   - Radial Points: " << gridRes.radialPoints << '\n';
            *console << "   - Lebedev Points: " << Quadrature::nLebedevPoints[(int) gridRes.lebedevOrder] << '\n';
            *console << "   - Atomic Grids: " << MolGrid.GetAtomicGrids().size() << '\n';
            *console << "   - Total Number of Points: " << MolGrid.GetGridSize() << '\n';
        }


        //Hamiltonian matrix
        MatrixXd H = T + V;
        // The Kohn-Sham Fock matrix
        MatrixXd F;
        // The one particle density matrices
        MatrixXd D = MatrixXd::Zero(S.rows(), S.cols());
        //The J coulomb matrix
        MatrixXd J = MatrixXd::Zero(S.rows(), S.cols());
        //The exchange-correlation matrix
        MatrixXd ExchangeCorr = MatrixXd::Zero(S.rows(), S.cols());
        // The wavefunction coefficient matrix
        MatrixXd C;
        // The Kohn-Sham orbital energies
        VectorXd Ep;

        //The total electronic energy of the system
        double E_KS = 0.0, lastE = 0.0;
        //Nuclear-nuclear repulsion energy
        double Vnn = mol.CalculateVnn();

        constexpr double E_CONVERGENCE = 1e-8;
        //RMS of the orbital gradient FDS - SDF required for convergence
        constexpr double COMMUTATOR_CONVERGENCE = 1e-5;
        constexpr uint32_t MAX_ITERATIONS = 500;
        constexpr size_t MAX_DIIS_SIZE = 6;

        //Pulay DIIS extrapolation of the Fock matrices, from the second Fock matrix
        HF::FockDIIS diis(MAX_DIIS_SIZE);
        double commutatorRMS = 1.0;
        bool converged = false;

        if (console)
            *console << "--- Starting SCF cycles ---\n";

        for (uint32_t iter = 0; iter < MAX_ITERATIONS; ++iter) {
            //Poll to check whether we should cancel the calculation
            if (thread)
                if (thread->ShouldStop())
                    return {};

            F = H + 2.0 * J + ExchangeCorr;

            MatrixXd F_Extrapolated = F;
            if (iter > 0) {
                //DIIS error vector (FDS - SDF) in the orthonormal basis
                MatrixXd error = U.transpose() * (F * D * S - S * D * F) * U;
                commutatorRMS = std::sqrt(error.cwiseAbs2().mean());

                diis.Add({F}, {error});
                if (diis.Size() >= 2)
                    F_Extrapolated = diis.Extrapolate()[0];
            }

            MatrixXd Fprime = U.transpose() * F_Extrapolated * U;

            SelfAdjointEigenSolver<MatrixXd> eigensolver2(Fprime);
            Ep = eigensolver2.eigenvalues().real();
            C = U * eigensolver2.eigenvectors().real();

            D = CreateDensityMatrix(C, nElectrons);

            MolGrid.SetDensity(D);

            J = Integrals::Transforms::CoulombMatrix(tei, D);
            //Exchange-correlation energy
            double E_XC = 0.0;
            CalculateXCMatrix(MolGrid, ExchangeCorr, E_XC);

            //Single electron energy
            double kineticPotentialE = 2.0 * (D * H).trace();
            //Coulomb repulsion energy
            double electronRepulsion = 2.0 * (D * J).trace();

            E_KS = kineticPotentialE + electronRepulsion + Vnn + E_XC;

            double dE = abs(E_KS - lastE);
            if (console)
                *console << "   DFT(" << iter << "): E(DFT) = " << E_KS << "          dE = " << dE
                        << "          |FDS-SDF| = " << commutatorRMS << '\n';

            if (iter > 0 && dE < E_CONVERGENCE && commutatorRMS < COMMUTATOR_CONVERGENCE) {
                converged = true;
                if (console) *console << "--- DFT procedure has converged! ---\n";
                break;
            }
            lastE = E_KS;
        }

        if (!converged && console)
            *console << "!!! DFT procedure failed to converge !!!\n";

        if (console) {
            *console << "--- Final DFT Energy ---\n";
            *console << "  E(DFT) = " << E_KS << " hartrees\n";
        }

        VectorXi Occ = VectorXi::Zero(H.rows());

        for (int64_t i = ((int64_t) mol.GetNoElectrons() / 2) - 1; i >= 0; --i)
            Occ(i) = 2;

        return KS_Wfn{E_KS, F, D, C, Ep, Occ};
    }
    //
    MatrixXd CreateDensityMatrix(const MatrixXd &C, uint32_t Nelecs) {
        const auto occupied = C.leftCols(Nelecs / 2);
        return occupied * occupied.transpose();
    }
    //
    void CalculateXCMatrix(const MolGrid &grid,
                           MatrixXd &ExchangeCorr,
                           double &E_XC) {
        VectorXd densities = grid.GetDensities();
        const VectorXd &weights = grid.GetWeights();
        //Basis function values, one row per grid point
        const MatrixXd &amplitudes = grid.GetAmplitudes();

        VectorXd ex; // exchange energy
        VectorXd vxa; // exchange potential (spin-up)
        VectorXd vxb; // exchange potential (spin-down)

        VectorXd ec; // correlation energy
        VectorXd vca; // correlation potential (spin-up)
        VectorXd vcb; // correlation potential (spin-down)

        VectorXd densitiesa = densities * 0.5;
        VectorXd densitiesb = densities * 0.5;

        // calculate exchange energies and potential
        XC::xalpha_x_functional(densitiesa, densitiesb, ex, vxa, vxb);

        // calculate correlation energy and potential
        XC::vwm_c_functional(densitiesa, densitiesb, ec, vca, vcb);

        // calculate total exchange-correlation energy
        E_XC = weights.dot(ex + ec);

        // calculate gridpoint-wise xc potential

        VectorXd wva = weights.cwiseProduct((vxa + vxb + vca + vcb) * 0.5);

        //V_xc(mu, nu) = sum_p w_p v_p phi_mu(p) phi_nu(p)
        ExchangeCorr = amplitudes.transpose() * (amplitudes.array().colwise() * wva.array()).matrix();
    }
    //
    //
}
