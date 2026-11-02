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
    void MP2_F12_Fixed(const char *obsName, const char *cabsName, const Molecule &mol,
                       const char *moleculeName) {
        auto atoms = mol.GetAtomsLibInt();

        F12BasisSet f12BasisSet(obsName, cabsName, atoms);

        fmt::print("Evaluating overlap, kinetic, nuclear & eri integrals...\n");
        Integrals::MolInts configInts = Integrals::F12::EvaluateMolecularIntegrals(f12BasisSet, atoms);
        fmt::print(fg(fmt::color::green), "Done!\n");

        fmt::print("Evaluating f12, f12^2, f12/r12 & |del f12|^2  integrals...\n");
        Integrals::F12::F12Integrals f12Ints = Integrals::F12::EvaluateF12Integrals(f12BasisSet);
        fmt::print(fg(fmt::color::green), "Done!\n");

        fmt::print("Obtaining UHF reference...\n");
        HF::HF_Wfn uhfWfn = HF::UHF(
            Integrals::F12::ExtractHFIntegrals(configInts,
                                               f12BasisSet.GetHFIndices().first,
                                               f12BasisSet.GetNoHFFuncs()),
            mol.GetNoElectrons() / 2,
            mol.GetNoElectrons() / 2,
            mol.CalculateVnn(),
            std::numeric_limits<double>::epsilon());

        F12::F12Idxs f12Idxs = F12::GenerateF12Idxs(f12BasisSet, uhfWfn);

        fmt::print("Canonicalising f12 integrals...\n");
        Eigen::MatrixXd F;
        Integrals::Transforms::TransformF12Basis(f12BasisSet,
                                                 uhfWfn.D_Alpha,
                                                 uhfWfn.C_Alpha,
                                                 configInts,
                                                 f12Ints,
                                                 true,
                                                 &F);

        fmt::print("Evaluating CABS singles correction...\n");

        double ECABSSingles = F12::EvaluateCABSSinglesCorrection(MatrixToTensor(F), f12Idxs);

        auto f12Result = F12::MP2_F12_Fixed(uhfWfn, configInts, f12Ints, f12Idxs);

        fmt::print("\n\n--- MP2 F12 Single Point Energy ---\n");
        fmt::print(" --- Fixed amplitude approximation used --- \n");
        fmt::print("  Molecule = {}\n", moleculeName ? moleculeName : "unknown");
        fmt::print("  E(HF) = {} hartrees\n", uhfWfn.E_HF);
        fmt::print("  E(CABS singles correction) = {} hartrees\n", ECABSSingles);
        fmt::print("  E(MP2) = {} hartrees\n", f12Result.E_MP2);
        fmt::print("  E(F12 MP2) = {} hartrees\n", f12Result.E_F12_MP2);
        fmt::print("  E(Total) = {} hartrees\n\n", uhfWfn.E_HF + ECABSSingles + f12Result.E_MP2 +
                                                   f12Result.E_F12_MP2);
    }
    //
    void MP2_F12_Optimise(const char *obsName, const char *cabsName, const Molecule &mol,
                          const char *moleculeName) {
        auto atoms = mol.GetAtomsLibInt();

        F12BasisSet f12BasisSet(obsName, cabsName, atoms);

        fmt::print("Evaluating overlap, kinetic, nuclear & eri integrals...\n");
        Integrals::MolInts configInts = Integrals::F12::EvaluateMolecularIntegrals(f12BasisSet, atoms);
        fmt::print(fg(fmt::color::green), "Done!\n");

        fmt::print("Evaluating f12, f12^2, f12/r12 & |del f12|^2  integrals...\n");
        Integrals::F12::F12Integrals f12Ints = Integrals::F12::EvaluateF12Integrals(f12BasisSet);
        fmt::print(fg(fmt::color::green), "Done!\n");

        fmt::print("Obtaining UHF reference...\n");
        HF::HF_Wfn uhfWfn = HF::UHF(
            Integrals::F12::ExtractHFIntegrals(configInts,
                                               f12BasisSet.GetHFIndices().first,
                                               f12BasisSet.GetNoHFFuncs()),
            mol.GetNoElectrons() / 2,
            mol.GetNoElectrons() / 2,
            mol.CalculateVnn(),
            std::numeric_limits<double>::epsilon());

        F12::F12Idxs f12Idxs = F12::GenerateF12Idxs(f12BasisSet, uhfWfn);

        fmt::print("Canonicalising f12 integrals...\n");
        Eigen::MatrixXd F;
        Integrals::Transforms::TransformF12Basis(f12BasisSet,
                                                 uhfWfn.D_Alpha,
                                                 uhfWfn.C_Alpha,
                                                 configInts,
                                                 f12Ints,
                                                 true,
                                                 &F);

        fmt::print("Evaluating CABS singles correction...\n");

        //Need to correct the reference energy with a perturbation by singly exciting into the CABS virtuals.
        //Correct the fact that there is a large basis set error in the non-correlated HF energy
        double ECABSSingles = F12::EvaluateCABSSinglesCorrection(MatrixToTensor(F), f12Idxs);

        auto f12Result = F12::MP2_F12_Optimise(uhfWfn, configInts, f12Ints, f12Idxs);

        fmt::print("\n\n--- MP2 F12 Single Point Energy ---\n");
        fmt::print(" --- Optimised amplitudes used --- \n");
        fmt::print("  Molecule = {}\n", moleculeName ? moleculeName : "unknown");
        fmt::print("  E(HF) = {} hartrees\n", uhfWfn.E_HF);
        fmt::print("  E(CABS singles correction) = {} hartrees\n", ECABSSingles);
        fmt::print("  E(MP2) = {} hartrees\n", f12Result.E_MP2);
        fmt::print("  E(F12 MP2) = {} hartrees\n", f12Result.E_F12_MP2);
        fmt::print("  E(Total) = {} hartrees\n\n", uhfWfn.E_HF + ECABSSingles + f12Result.E_MP2 +
                                                   f12Result.E_F12_MP2);
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
        auto atoms = mol.GetAtomsLibInt();

        F12BasisSet f12BasisSet(obsName, cabsName, atoms);

        fmt::print("Evaluating overlap, kinetic, nuclear & eri integrals...\n");
        Integrals::MolInts configInts = Integrals::F12::EvaluateMolecularIntegrals(f12BasisSet, atoms);
        fmt::print(fg(fmt::color::green), "Done!\n");

        fmt::print("Evaluating f12, f12^2, f12/r12 & |del f12|^2  integrals...\n");
        Integrals::F12::F12Integrals f12Ints = Integrals::F12::EvaluateF12Integrals(f12BasisSet);
        fmt::print(fg(fmt::color::green), "Done!\n");

        fmt::print("Obtaining UHF reference...\n");
        const uint32_t alphaElectrons = mol.GetNoElectrons() / 2;
        const uint32_t &betaElectrons = alphaElectrons;
        HF::HF_Wfn uhfWfn = HF::UHF(
            Integrals::F12::ExtractHFIntegrals(configInts, 0,
                                               f12BasisSet.GetNoHFFuncs()),
            alphaElectrons,
            betaElectrons,
            mol.CalculateVnn(),
            std::numeric_limits<double>::epsilon());

        F12::F12Idxs f12Idxs = F12::GenerateF12Idxs(f12BasisSet, uhfWfn);

        fmt::print("Canonicalising f12 integrals...\n");
        Eigen::MatrixXd F;
        Integrals::Transforms::TransformF12Basis(f12BasisSet,
                                                 uhfWfn.D_Alpha,
                                                 uhfWfn.C_Alpha,
                                                 configInts,
                                                 f12Ints,
                                                 true,
                                                 &F);

        fmt::print("Evaluating CABS singles correction...\n");
        //Need to correct the reference energy with a perturbation by singly exciting into the CABS virtuals.
        //Correct the fact that there is a large basis set error in the non-correlated HF energy
        double ECABSSingles = F12::EvaluateCABSSinglesCorrection(MatrixToTensor(F), f12Idxs);

        double ECCQMC_F12 = F12_QMC::CCSD_QMC_MP2_F12(uhfWfn, configInts, f12Ints, f12Idxs);
        double ECorrelation = ECCQMC_F12 - uhfWfn.E_HF;

        fmt::print(fg(fmt::color::white), "\n\n--- CCSD QMC F12 Single Point Energy ---\n");
        fmt::print(fg(fmt::color::white), " --- Fixed amplitude approximation used --- \n");
        fmt::print(fg(fmt::color::white), "  Molecule = {}\n", moleculeName ? moleculeName : "unknown");
        fmt::print(fg(fmt::color::white), "  E(HF) = {:.8} hartrees\n",
                   uhfWfn.E_HF);
        fmt::print(fg(fmt::color::white), "  E(CABS singles correction) = {:.8} hartrees\n",
                   ECABSSingles);
        fmt::print(fg(fmt::color::white), "  E(CCQMC_F12) = {:.8} hartrees\n",
                   ECorrelation);
        fmt::print(fg(fmt::color::green), "  E(Total) = {:.8} hartrees\n\n",
                   uhfWfn.E_HF +
                   ECABSSingles +
                   ECorrelation);
    }
    //
    //
#pragma endregion
    //
    //
}
