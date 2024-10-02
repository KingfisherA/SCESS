#pragma once

#include "Eigen/Dense"
#include "Chemistry/Quantum/Integrals/Integrals.h"

namespace SCESS {
    class GaussianCollection;
    struct RHFResult;
    struct UHFResult;
    class Molecule;

    class HartreeFock {
    public:
        /// <summary>
        /// Computes the restricted closed-shell hartree fock energy and orbitals in the given basis set defined by the molecule.
        /// </summary>
        /// <param name="molecule">The molecule to compute RHF on. This molecule must be finalised with a set of orbital basis functions</param>
        /// <param name="convergeLimit">Convergence limit of the SCF procedure</param>
        /// <param name="scfMaxIterations">Maximum amount of scf iterations to perform. The iterations are automatically terminated when the convergence limit has been achieved.</param>
        /// <param name="console">An optional pointer to a string stream to output progress information to</param>
        /// <returns>A RHFResult structure with MO coefficients and energies</returns>
        static RHFResult Restricted(const Molecule &molecule, const double convergeLimit = 1E-5,
                                    const size_t scfMaxIterations = 100, std::stringstream *console = nullptr,
                                    const std::string_view basisSet = "");

        /// <summary>
        /// Computes the unrestricted open-shell hartree fock energy and orbitals in the given basis set defined by the molecule.
        /// </summary>
        /// <param name="molecule">The molecule to compute RHF on. This molecule must be finalised with a set of orbital basis functions</param>
        /// <param name="convergeLimit">Convergence limit of the SCF procedure</param>
        /// <param name="scfMaxIterations">Maximum amount of scf iterations to perform. The iterations are automatically terminated when the convergence limit has been achieved.</param>
        /// <param name="console">An optional pointer to a string stream to output progress information to</param>
        /// <returns>A RHFResult structure with MO coefficients and energies</returns>
        static UHFResult Unrestricted(const Molecule &molecule, const double convergeLimit = 1E-5,
                                      const size_t scfMaxIterations = 100, std::stringstream *console = nullptr,
                                      const std::string_view basisSet = "");

        /// <summary>
        /// Evaluate the energy for a given density, hamiltonian and fock matrix for a closed-shell RHF procedure.
        /// </summary>
        /// <param name="Density">RHF density matrix</param>
        /// <param name="Hcore">Hamiltonian matrix</param>
        /// <param name="Fock">RHF Fock matrix</param>
        /// <returns>Computed energy value</returns>
        static double CalculateEnergyRHF(const Eigen::MatrixXd &Density,
                                         const Eigen::MatrixXd &Hcore, const Eigen::MatrixXd &Fock);

        /// <summary>
        /// Evaluate the energy for a given density, hamiltonian and fock matrix for an open-shell UHF procedure.
        /// </summary>
        /// <param name="DensityAlpha">Density matrix for spin-up electrons</param>
        /// <param name="DensityBeta">Density matrix for spin-down electrons</param>
        /// <param name="Hcore">Hamiltonian matrix</param>
        /// <param name="FockAlpha">Fock matrix for spin-up electrons</param>
        /// <param name="FockBeta">Fock matrix for spin-down electrons</param>
        /// <returns>Computed energy value</returns>
        static double CalculateEnergyUHF(
            const Eigen::MatrixXd &DensityAlpha,
            const Eigen::MatrixXd &DensityBeta,
            const Eigen::MatrixXd &Hcore,
            const Eigen::MatrixXd &FockAlpha,
            const Eigen::MatrixXd &FockBeta);

        /// <summary>
        /// Compute the closed-shell restricted HF Moller-Plesset Order 2 correction energy
        /// </summary>
        /// <param name="result">Result from a RHF computation</param>
        /// <returns>MP2 correction energy</returns>
        static double MP2ClosedShell(const RHFResult &result, std::stringstream *console = nullptr);
        /**
         * Compute the closed-shell restricted HF Moller-Plesset Order 2 correction energy
         * @param result Result from a RHF computation
         * @param moTEI The tei transformed into a molecular orbital basis in chemist's notation
         * @return MP2 correction energy
         */
        static double MP2ClosedShell(const RHFResult &result, const TEITensor &moTEI,std::stringstream *console = nullptr);

    private:
        static Eigen::MatrixXd ConstructGRHF(const Eigen::MatrixXd &Density, const TEITensor &teiTensor);
        static std::pair<Eigen::MatrixXd, Eigen::MatrixXd> ConstructGUHF(const Eigen::MatrixXd &DensityAlpha,
                                                                         const Eigen::MatrixXd &DensityBeta,
                                                                         const TEITensor &teiTensor);


        static Eigen::MatrixXd ConstructDensityRHF(const Eigen::MatrixXd &Coeffs,
                                                   const Eigen::Index Nelecs);

        static double CalculateDelta(const Eigen::MatrixXd &oldDensity,
                                     const Eigen::MatrixXd &density);
    };


    struct HFResult {
        HFResult() : Energy(0), Integrals() {
        }
        HFResult(const double energyIn, const ConfigIntegrals &integrals) : Energy(energyIn), Integrals(integrals) {
        }
        HFResult(const HFResult &res) = default;


        const double Energy;
        const ConfigIntegrals Integrals;
    };

    struct RHFResult : public HFResult {
        RHFResult() : HFResult() {
        }
        RHFResult(
            const Eigen::MatrixXd &Fock,
            const Eigen::MatrixXd &Hcore,
            const Eigen::MatrixXd &Density,
            const Eigen::MatrixXd &Coeffs,
            const Eigen::VectorXd &Energies,
            const Eigen::VectorXi &Occupation,
            const ConfigIntegrals &integrals,
            const double Energy) : Fock(Fock), Hcore(Hcore),
                                   Density(Density), Coeffs(Coeffs), Energies(Energies), Occupation(Occupation),
                                   HFResult(Energy, integrals) {
        }
        RHFResult(const RHFResult &) = default;
        Eigen::MatrixXd Fock;
        Eigen::MatrixXd Hcore;
        Eigen::MatrixXd Density;
        Eigen::MatrixXd Coeffs;
        Eigen::VectorXd Energies;
        Eigen::VectorXi Occupation;
    };

    struct UHFResult : public HFResult {
        UHFResult() : HFResult() {
        }
        UHFResult(
            const Eigen::MatrixXd &FockAlpha,
            const Eigen::MatrixXd &FockBeta,
            const Eigen::MatrixXd &Hcore,
            const Eigen::MatrixXd &DensityAlpha,
            const Eigen::MatrixXd &DensityBeta,
            const Eigen::MatrixXd &CoeffsAlpha,
            const Eigen::MatrixXd &CoeffsBeta,
            const Eigen::VectorXd &EnergiesAlpha,
            const Eigen::VectorXd &EnergiesBeta,
            const Eigen::VectorXi &OccupationAlpha,
            const Eigen::VectorXi &OccupationBeta,
            const ConfigIntegrals &integrals,
            const double Energy) : FockAlpha(FockAlpha), FockBeta(FockBeta), Hcore(Hcore),
                                   DensityAlpha(DensityAlpha), DensityBeta(DensityBeta), CoeffsAlpha(CoeffsAlpha),
                                   CoeffsBeta(CoeffsBeta), EnergiesAlpha(EnergiesAlpha), EnergiesBeta(EnergiesBeta),
                                   OccupationAlpha(OccupationAlpha), OccupationBeta(OccupationBeta),
                                   HFResult(Energy, integrals) {
        }
        UHFResult(const UHFResult &res) = default;
        Eigen::MatrixXd FockAlpha, FockBeta;
        Eigen::MatrixXd Hcore;
        Eigen::MatrixXd DensityAlpha, DensityBeta;
        Eigen::MatrixXd CoeffsAlpha, CoeffsBeta;
        Eigen::VectorXd EnergiesAlpha, EnergiesBeta;
        Eigen::VectorXi OccupationAlpha, OccupationBeta;
    };
}
