#include "Driver.h"
//
#include <fmt/format.h>
#include <fmt/color.h>
//
#include "Quantum/Integrals/IntegralTransforms.h"
#include "Quantum/Integrals/Integrals.h"
#include "Quantum/Integrals/IntegralsF12.h"
//
#include "Quantum/Molecule.h"
#include "Quantum/Basis/BasisSet.h"
#include "Quantum/DFT/DFT.h"
#include "Quantum/HF/HF.h"
#include "Quantum/CC/CC.h"
#include "Quantum/F12/F12.h"
#include "Quantum/F12_QMC/F12_QMC.h"
//
#include "Quantum/QMC/QMC.h"
#include "Tensor/TensorOperations.h"
//


//
//
namespace SCESS::Driver {
    //
#pragma region Self-Consistent field theories
    //
    //
    void DFT(const char *basisSetName, const Molecule &mol, const char *moleculeName) {
        fmt::print(fmt::emphasis::bold, "--- Performing Kohn-Sham DFT calculation on {} ---\n",
                   moleculeName ? moleculeName : "unknown");

        BasisSet basisSet(basisSetName, mol);

        std::stringstream console;
        DFT::KS_Wfn ksWfn = DFT::KS_DFT(mol,
                                        basisSet,
                                        Quadrature::MolGrid::MEDIUM_RES,
                                        nullptr,
                                        &console);
        fmt::print("{}\n", console.str());

        fmt::print("\n\n--- Kohn-Sham DFT Single Point Energy ---\n");
        fmt::print("    Correlation Functional = VMW\n");
        fmt::print("    Exchange Functional = Slater-X\n");
        fmt::print("  Molecule = {}\n", moleculeName ? moleculeName : "unknown");
        fmt::print("  E(DFT) = {} hartrees\n", ksWfn.E_KS);
    }
    //
    void RHF(const char *basisSetName, const Molecule &mol, const char *moleculeName) {
        fmt::print(fmt::emphasis::bold, "--- Performing RHF calculation on {} ---\n",
                   moleculeName ? moleculeName : "unknown");

        BasisSet basisSet(basisSetName, mol);

        std::stringstream console;
        HF::HF_Wfn rhfWfn = HF::RHF(mol, basisSet, 1E-8, 100,
                                    nullptr, nullptr, &console);
        fmt::print("{}\n", console.str());

        fmt::print("\n\n--- Restricted Hartree-Fock Single Point Energy ---\n");
        fmt::print("  Molecule = {}\n", moleculeName ? moleculeName : "unknown");
        fmt::print("  E(RHF) = {} hartrees\n", rhfWfn.E_HF);
    }
    //
    void UHF(const char *basisSetName, const Molecule &mol, const char *moleculeName) {
        fmt::print(fmt::emphasis::bold, "--- Performing UHF calculation on {} ---\n",
                   moleculeName ? moleculeName : "unknown");

        BasisSet basisSet(basisSetName, mol);

        std::stringstream console;
        HF::HF_Wfn uhfWfn = HF::UHF(mol, basisSet, 1E-8, 100,
                                    nullptr, nullptr, &console);
        fmt::print("{}\n", console.str());

        fmt::print("\n\n--- Unrestricted Hartree-Fock Single Point Energy ---\n");
        fmt::print("  Molecule = {}\n", moleculeName ? moleculeName : "unknown");
        fmt::print("  E(UHF) = {} hartrees\n", uhfWfn.E_HF);
    }
    //
    //
#pragma endregion
    //
    //
#pragma region Moller Plesset Perturbation Theory
    //
    //
    void MP2(const char *basisSetName, const Molecule &mol, const char *moleculeName) {
        fmt::print(fmt::emphasis::bold, "--- Performing MP2 calculation on {} ---\n",
                   moleculeName ? moleculeName : "unknown");

        BasisSet basisSet(basisSetName, mol);

        std::stringstream console;
        Integrals::MolInts mints;
        HF::HF_Wfn rhfWfn = HF::RHF(mol, basisSet, 1E-8, 100,
                                    &mints, nullptr, &console);
        fmt::print("{}\n", console.str());
        console.str("");

        double E_MP2 = HF::MP2_RHF_ClosedShellAO(rhfWfn, mints.eri, nullptr, &console);
        fmt::print("{}\n", console.str());

        fmt::print("\n\n--- Moller-Plesset Second Order Single Point Energy ---\n");
        fmt::print("  Molecule = {}\n", moleculeName ? moleculeName : "unknown");
        fmt::print("  E(RHF) = {} hartrees\n", rhfWfn.E_HF);
        fmt::print("  E(MP2) = {} hartrees\n", E_MP2);
        fmt::print("  E(Total) = {} hartrees\n", rhfWfn.E_HF + E_MP2);
    }
    //
    //
#pragma endregion
    //
    //
#pragma region Configuration Interaction
    //
    //
    void FCI(const char *basisSetName, const Molecule &mol, const char *moleculeName) {
        fmt::print(fmt::emphasis::bold, "--- Performing FCI calculation on {} ---\n",
                   moleculeName ? moleculeName : "unknown");

        BasisSet basisSet(basisSetName, mol);

        std::stringstream console;
        Integrals::MolInts mints;
        HF::HF_Wfn rhfWfn = HF::RHF(mol, basisSet, 1E-8, 100,
                                    &mints, nullptr, &console);
        fmt::print("{}\n", console.str());
        console.str("");

        CI::CI_Wfn ciWfn = CI::CI(rhfWfn, mints, mol.CalculateVnn(), -1, nullptr, &console);
        fmt::print("{}\n", console.str());

        fmt::print("\n\n--- Full Configuration Interaction Single Point Energy ---\n");
        fmt::print("  Molecule = {}\n", moleculeName ? moleculeName : "unknown");
        fmt::print("  E(RHF) = {} hartrees\n", rhfWfn.E_HF);
        fmt::print("  E(FCI correlation) = {} hartrees\n", ciWfn.E_CI(0) - rhfWfn.E_HF);
        fmt::print("  E(Total) = {} hartrees\n", ciWfn.E_CI(0));
    }
    //
    void FCIQMC(const char *basisSetName, const Molecule &mol, const char *moleculeName) {
        fmt::print(fmt::emphasis::bold, "--- Performing FCIQMC calculation on {} ---\n",
                   moleculeName ? moleculeName : "unknown");

        BasisSet basisSet(basisSetName, mol);

        std::stringstream console;
        Integrals::MolInts mints;
        HF::HF_Wfn rhfWfn = HF::RHF(mol, basisSet, 1E-8, 100,
                                    &mints, nullptr, &console);
        fmt::print("{}\n", console.str());
        console.str("");

        double E_FCI = QMC::FCI::FCIQMC(rhfWfn, mints, mol.CalculateVnn(), nullptr, &console);
        fmt::print("{}\n", console.str());

        fmt::print("\n\n--- Full Configuration Interaction Quantum Monte Carlo Single Point Energy ---\n");
        fmt::print("  Molecule = {}\n", moleculeName ? moleculeName : "unknown");
        fmt::print("  E(RHF) = {} hartrees\n", rhfWfn.E_HF);
        fmt::print("  E(FCIQMC correlation) = {} hartrees\n", E_FCI - rhfWfn.E_HF);
        fmt::print("  E(Total) = {} hartrees\n", E_FCI);
    }
    //
    //
#pragma endregion
    //
    //
#pragma region Coupled-Cluster
    //
    //
    void CCSD(const char *basisSetName, const Molecule &mol, const char *moleculeName) {
        fmt::print(fmt::emphasis::bold, "--- Performing CCSD calculation on {} ---\n",
                   moleculeName ? moleculeName : "unknown");

        BasisSet basisSet(basisSetName, mol);

        Integrals::MolInts mints;
        std::stringstream console;
        HF::HF_Wfn rhfWfn = HF::RHF(mol,
                                    basisSet,
                                    1e-08,
                                    100,
                                    &mints,
                                    nullptr,
                                    &console);
        fmt::print("{}\n", console.str());
        console.str("");

        CC::CCSD_Wfn ccsdWfn = CC::CCSD(rhfWfn, mints, nullptr, &console);

        fmt::print("{}\n", console.str());
        console.str("");

        fmt::print("\n\n--- Coupled-Cluster Singles & Doubles Single Point Energy ---\n");
        fmt::print("  Molecule = {}\n", moleculeName ? moleculeName : "unknown");
        fmt::print("  E(RHF) = {} hartrees\n", rhfWfn.E_HF);
        fmt::print("  E(CCSD correlation) = {} hartrees\n", ccsdWfn.E_CCSD);
        fmt::print("  E(Total) = {} hartrees\n", ccsdWfn.E_CCSD + rhfWfn.E_HF);
    }
    //
    void CCSD_T(const char *basisSetName, const Molecule &mol, const char *moleculeName) {
        fmt::print(fmt::emphasis::bold, "--- Performing CCSD(T) calculation on {} ---\n",
                   moleculeName ? moleculeName : "unknown");

        BasisSet basisSet(basisSetName, mol);

        Integrals::MolInts mints;
        std::stringstream console;
        HF::HF_Wfn rhfWfn = HF::RHF(mol,
                                    basisSet,
                                    1e-08,
                                    100,
                                    &mints,
                                    nullptr,
                                    &console);
        fmt::print("{}\n", console.str());
        console.str("");

        CC::CCSD_Wfn ccsdWfn = CC::CCSD_T(rhfWfn, mints, nullptr, &console);

        fmt::print("{}\n", console.str());
        console.str("");

        fmt::print("\n\n--- CCSD(T) Single Point Energy ---\n");
        fmt::print("  Molecule = {}\n", moleculeName ? moleculeName : "unknown");
        fmt::print("  E(RHF) = {} hartrees\n", rhfWfn.E_HF);
        fmt::print("  E(CCSD(T) correlation) = {} hartrees\n", ccsdWfn.E_CCSD + ccsdWfn.E_T);
        fmt::print("  E(Total) = {} hartrees\n", ccsdWfn.E_CCSD + rhfWfn.E_HF + ccsdWfn.E_T);
    }
    //
    void CCQMC(const char *basisSetName, const Molecule &mol, const char *moleculeName) {
        fmt::print(fmt::emphasis::bold, "--- Performing CCQMC calculation on {} ---\n",
                   moleculeName ? moleculeName : "unknown");

        BasisSet basisSet(basisSetName, mol);

        std::stringstream console;
        Integrals::MolInts mints;
        HF::HF_Wfn rhfWfn = HF::RHF(mol, basisSet, 1E-8, 100,
                                    &mints, nullptr, &console);
        fmt::print("{}\n", console.str());
        console.str("");

        double E_CC = QMC::CC::CCQMC(rhfWfn, mints);

        fmt::print("\n\n--- Coupled-Cluster Quantum Monte Carlo Single Point Energy ---\n");
        fmt::print("  Molecule = {}\n", moleculeName ? moleculeName : "unknown");
        fmt::print("  E(RHF) = {} hartrees\n", rhfWfn.E_HF);
        fmt::print("  E(CCQMC correlation) = {} hartrees\n", E_CC - rhfWfn.E_HF);
        fmt::print("  E(Total) = {} hartrees\n", E_CC);
    }
    //
    //
#pragma endregion
    //
    //
#pragma region F12 theories
    //
    //
    namespace {
        //Prepares a F12 calculation, printing its progress
        F12::F12Calculation PrepareF12(const char *obsName, const char *cabsName, const Molecule &mol,
                                       bool buildQMCIntegrals = false) {
            std::stringstream console;
            F12::F12Calculation calc = F12::PrepareF12Calculation(obsName, cabsName, mol, nullptr, &console,
                                                                  buildQMCIntegrals);
            fmt::print("{}", console.str());
            return calc;
        }
    }
    //
    void MP2_F12_Fixed(const char *obsName, const char *cabsName, const Molecule &mol,
                       const char *moleculeName) {
        const F12::F12Calculation calc = PrepareF12(obsName, cabsName, mol);

        auto f12Result = F12::MP2_F12_Fixed(calc.hfWfn, calc.ints, calc.idxs);

        fmt::print("\n\n--- MP2 F12 Single Point Energy ---\n");
        fmt::print(" --- Fixed amplitude approximation used --- \n");
        fmt::print("  Molecule = {}\n", moleculeName ? moleculeName : "unknown");
        fmt::print("  OBS = {}, CABS = {}\n", obsName, cabsName);
        fmt::print("  E(HF) = {} hartrees\n", calc.hfWfn.E_HF);
        fmt::print("  E(CABS singles correction) = {} hartrees\n", calc.E_CABSSingles);
        fmt::print("  E(MP2) = {} hartrees\n", f12Result.E_MP2);
        fmt::print("  E(F12 MP2) = {} hartrees\n", f12Result.E_F12_MP2);
        fmt::print("  E(Total) = {} hartrees\n\n", calc.hfWfn.E_HF + calc.E_CABSSingles + f12Result.E_MP2 +
                                                   f12Result.E_F12_MP2);
    }
    //
    void MP2_F12_Optimise(const char *obsName, const char *cabsName, const Molecule &mol,
                          const char *moleculeName) {
        const F12::F12Calculation calc = PrepareF12(obsName, cabsName, mol);

        auto f12Result = F12::MP2_F12_Optimise(calc.hfWfn, calc.ints, calc.idxs);

        fmt::print("\n\n--- MP2 F12 Single Point Energy ---\n");
        fmt::print(" --- Optimised amplitudes used --- \n");
        fmt::print("  Molecule = {}\n", moleculeName ? moleculeName : "unknown");
        fmt::print("  OBS = {}, CABS = {}\n", obsName, cabsName);
        fmt::print("  E(HF) = {} hartrees\n", calc.hfWfn.E_HF);
        fmt::print("  E(CABS singles correction) = {} hartrees\n", calc.E_CABSSingles);
        fmt::print("  E(MP2) = {} hartrees\n", f12Result.E_MP2);
        fmt::print("  E(F12 MP2) = {} hartrees\n", f12Result.E_F12_MP2);
        fmt::print("  E(Total) = {} hartrees\n\n", calc.hfWfn.E_HF + calc.E_CABSSingles + f12Result.E_MP2 +
                                                   f12Result.E_F12_MP2);
    }
    //
    void CCSD_F12b(const char *obsName, const char *cabsName, const Molecule &mol,
                   const char *moleculeName) {
        fmt::print(fmt::emphasis::bold, "--- Performing CCSD-F12b calculation on {} ---\n",
                   moleculeName ? moleculeName : "unknown");

        const F12::F12Calculation calc = PrepareF12(obsName, cabsName, mol);

        std::stringstream console;
        const F12::CCSD_F12bResult result = F12::CCSD_F12b(calc, nullptr, &console);
        fmt::print("{}\n", console.str());

        if (!result.converged) {
            fmt::print(fg(fmt::color::crimson), "CCSD-F12b failed to converge\n");
            return;
        }

        fmt::print("\n\n--- CCSD-F12b Single Point Energy ---\n");
        fmt::print(" --- 3C(FIX) ansatz with the CABS approximation --- \n");
        fmt::print("  Molecule = {}\n", moleculeName ? moleculeName : "unknown");
        fmt::print("  OBS = {}, CABS = {}\n", obsName, cabsName);
        fmt::print("  E(HF) = {:.10f} hartrees\n", calc.hfWfn.E_HF);
        fmt::print("  E(CABS singles correction) = {:.10f} hartrees\n", calc.E_CABSSingles);
        fmt::print("  E(MP2) = {:.10f} hartrees\n", result.E_MP2);
        fmt::print("  E(CCSD, with geminal coupling) = {:.10f} hartrees\n", result.E_CCSD);
        fmt::print("      of which E(coupling) = {:.10f} hartrees\n", result.E_Coupling);
        fmt::print("  E(F12 geminals) = {:.10f} hartrees\n", result.E_F12);
        fmt::print("  E(CCSD-F12b correlation) = {:.10f} hartrees\n", result.CorrelationEnergy());
        fmt::print(fg(fmt::color::green), "  E(Total) = {:.10f} hartrees\n\n",
                   calc.hfWfn.E_HF + calc.E_CABSSingles + result.CorrelationEnergy());
    }
    //
    //
#pragma endregion
    //
    //
#pragma region F12-QMC theories
    //
    //
    void CCSD_QMC_MP2_F12(const char *obsName, const char *cabsName,
                          const Molecule &mol, const char *moleculeName) {
        const F12::F12Calculation calc = PrepareF12(obsName, cabsName, mol, true);

        double ECCQMC_F12 = F12_QMC::CCSD_QMC_MP2_F12(calc.hfWfn, calc.ints, calc.idxs);
        double ECorrelation = ECCQMC_F12 - calc.hfWfn.E_HF;

        fmt::print(fg(fmt::color::white), "\n\n--- CCSD QMC F12 Single Point Energy ---\n");
        fmt::print(fg(fmt::color::white), " --- Fixed amplitude approximation used --- \n");
        fmt::print(fg(fmt::color::white), "  Molecule = {}\n", moleculeName ? moleculeName : "unknown");
        fmt::print(fg(fmt::color::white), "  E(HF) = {:.8} hartrees\n",
                   calc.hfWfn.E_HF);
        fmt::print(fg(fmt::color::white), "  E(CABS singles correction) = {:.8} hartrees\n",
                   calc.E_CABSSingles);
        fmt::print(fg(fmt::color::white), "  E(CCQMC_F12) = {:.8} hartrees\n",
                   ECorrelation);
        fmt::print(fg(fmt::color::green), "  E(Total) = {:.8} hartrees\n\n",
                   calc.hfWfn.E_HF +
                   calc.E_CABSSingles +
                   ECorrelation);
    }
    //
    //
#pragma endregion
    //
    //
}
