#pragma once
#include <Eigen/Dense>
#include "Quantum/Integrals/Integrals.h"

namespace SCESS {
    class BasisSet;
    class CGF;
    class Molecule;
}

namespace SCESS::Integrals {
    struct MolInts;
}

namespace DH {
    class Thread;
}

//Performs Hartree-Fock Calculations
namespace SCESS::HF {
    using TEITensor = Eigen::Tensor<double, 4>;
    //
    //Represents a Hartree-Fock wavefunction.
    //If RHF then only the alpha parameters are set except for occupation where the alpha and beta are always both set.
    //UHF HF wavefunctions have both the alpha and beta parameters set
    struct HF_Wfn {
        //The total electronic energy of the Hartree Fock calculation including the nuclear-nuclear repulsion
        //term Vnn
        double E_HF = 0.0;

        //The alpha Fock matrices
        Eigen::MatrixXd F_Alpha;
        Eigen::MatrixXd F_Beta;

        //The alpha density matrix
        Eigen::MatrixXd D_Alpha;
        Eigen::MatrixXd D_Beta;

        //The HF molecular orbital coefficients
        Eigen::MatrixXd C_Alpha;
        Eigen::MatrixXd C_Beta;

        //The Hartree-Fock orbital energy eigenvalues
        Eigen::VectorXd Ep_Alpha;
        Eigen::VectorXd Ep_Beta;

        //The alpha and beta electron occupation numbers
        Eigen::VectorXi Occ_Alpha;
        Eigen::VectorXi Occ_Beta;

        bool isUHF = false;

        //
        int GetTotalNoElectrons() const {
            int elec = 0;

            for (const auto &ele: Occ_Alpha)
                elec += ele;

            for (const auto &ele: Occ_Beta)
                elec += ele;

            return elec;
        }
        //
        //
    };
    //
    //
    /// <summary>
    /// Computes the restricted closed-shell hartree fock energy and orbitals in the given basis set defined by the molecule.
    /// </summary>
    /// <param name="molecule">The molecule to compute RHF on. This molecule must be finalised with a set of orbital basis functions</param>
    /// <param name="convergeLimit">Convergence limit of the SCF procedure</param>
    /// <param name="scfMaxIterations">Maximum amount of scf iterations to perform. The iterations are automatically terminated when the convergence limit has been achieved.</param>
    /// <param name="thread">The parent thread running this calculation</param>
    /// <param name="console">An optional pointer to a string stream to output progress information to</param>
    /// <param name="basisSetName">The basis set name for this calculation</param>
    /// <returns>A RHFResult structure with MO coefficients and energies</returns>
    HF_Wfn RHF(const Molecule &molecule,
               const BasisSet &basis,
               double convergeLimit = 1E-5,
               size_t scfMaxIterations = 100,
               Integrals::MolInts *mintsOut = nullptr,
               const DH::Thread *thread = nullptr,
               std::stringstream *console = nullptr);
    //
    HF_Wfn RHF(const Integrals::MolInts &mints,
               uint32_t nElectronPairs,
               double Vnn,
               double convergeLimit = 1E-5,
               size_t scfMaxIterations = 100,
               const DH::Thread *thread = nullptr,
               std::stringstream *console = nullptr);
    //
    /// <summary>
    /// Computes the unrestricted open-shell hartree fock energy and orbitals in the given basis set defined by the molecule.
    /// </summary>
    /// <param name="molecule">The molecule to compute RHF on. This molecule must be finalised with a set of orbital basis functions</param>
    /// <param name="basis">The basis set for this calculation</param>
    /// <param name="convergeLimit">Convergence limit of the SCF procedure</param>
    /// <param name="scfMaxIterations">Maximum amount of scf iterations to perform. The iterations are automatically terminated when the convergence limit has been achieved.</param>
    /// <param name="thread">The parent thread running this calculation</param>
    /// <param name="console">An optional pointer to a string stream to output progress information to</param>
    /// <returns>A RHFResult structure with MO coefficients and energies</returns>
    HF_Wfn UHF(const Molecule &molecule,
               const BasisSet &basis,
               double convergeLimit = 1E-5,
               size_t scfMaxIterations = 100,
               Integrals::MolInts *mintsOut = nullptr,
               const DH::Thread *thread = nullptr,
               std::stringstream *console = nullptr);
    //
    HF_Wfn UHF(const Integrals::MolInts &mints,
               uint32_t alphaElectrons,
               uint32_t betaElectrons,
               double Vnn,
               double convergeLimit = 1E-5,
               size_t scfMaxIterations = 100,
               const DH::Thread *thread = nullptr,
               std::stringstream *console = nullptr);
    //
    double CalculateEnergyRHF(const Eigen::MatrixXd &Density,
                              const Eigen::MatrixXd &Hcore,
                              const Eigen::MatrixXd &Fock);

    /// <summary>
    /// Evaluate the energy for a given density, hamiltonian and fock matrix for an open-shell UHF procedure.
    /// </summary>
    /// <param name="DensityAlpha">Density matrix for spin-up electrons</param>
    /// <param name="DensityBeta">Density matrix for spin-down electrons</param>
    /// <param name="Hcore">Hamiltonian matrix</param>
    /// <param name="FockAlpha">Fock matrix for spin-up electrons</param>
    /// <param name="FockBeta">Fock matrix for spin-down electrons</param>
    /// <returns>Computed energy value</returns>
    double CalculateEnergyUHF(
        const Eigen::MatrixXd &DensityAlpha,
        const Eigen::MatrixXd &DensityBeta,
        const Eigen::MatrixXd &Hcore,
        const Eigen::MatrixXd &FockAlpha,
        const Eigen::MatrixXd &FockBeta);
    //
    /// <summary>
    /// Compute the closed-shell restricted HF Moller-Plesset Order 2 correction energy
    /// </summary>
    /// <param name="result">Result from a RHF computation</param>
    /// <returns>MP2 correction energy</returns>
    double MP2_RHF_ClosedShellAO(const HF_Wfn &wfn,
                                 const TEITensor &aoERI,
                                 const DH::Thread *thread = nullptr,
                                 std::stringstream *console = nullptr);
    //
    /**
     * Compute the closed-shell restricted HF Moller-Plesset Order 2 correction energy
     * @param result Result from a RHF computation
     * @param moTEI The tei transformed into a molecular orbital basis in chemist's notation
     * @return MP2 correction energy
     */
    double MP2_RHF_ClosedShell(const HF_Wfn &wfn,
                               const TEITensor &moERI,
                               const DH::Thread *thread = nullptr,
                               std::stringstream *console = nullptr);
    //
    //
}
