#include "Driver.h"
#include<fmt/format.h>
#include <fmt/color.h>

#include "Quantum/Molecule.h"
#include "Quantum/Basis/BasisSet.h"
#include "Quantum/HF/HF.h"
#include "Quantum/Integrals/IntegralTransforms.h"
#include "Quantum/Integrals/Integrals.h"
#include "Quantum/Integrals/IntegralsF12.h"
#include "Quantum/F12/F12.h"
#include "Quantum/F12_QMC/F12_QMC.h"
#include "Tensor/TensorOperations.h"


//
//
namespace SCESS::Driver {
    //
    void MP2_F12_Fixed_H2O() {
        Molecule water = Molecule::GenerateWaterMolecule();
        MP2_F12_Fixed("sto-3g", "cc-pVDZ", water,
                      "Water (H2O)");
    }
    void MP2_F12_Fixed_LiH() {
        Molecule lithiumHydride = Molecule::GenerateLithiumHydrideMolecule();
        MP2_F12_Fixed("sto-3g", "cc-pVDZ", lithiumHydride,
                      "Lithium Hydride (LiH)");
    }
    //
    void MP2_F12_Fixed(const char *obsName, const char *cabsName, const Molecule &mol,
                       const char *moleculeName) {
        auto atoms = mol.GetAtomsLibInt();

        F12BasisSet f12BasisSet(obsName, cabsName, atoms);

        fmt::print("Evaluating overlap, kinetic, nuclear & eri integrals...\n");
        Integrals::MolInts configInts = Integrals::F12::EvaluateMolecularIntegrals(f12BasisSet, atoms);
        fmt::print("Done!\n");

        fmt::print("Evaluating f12, f12^2, f12/r12 & |del f12|^2  integrals...\n");
        Integrals::F12::F12Integrals f12Ints = Integrals::F12::EvaluateF12Integrals(f12BasisSet);
        fmt::print("Done!\n");

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

        //configInts.eri = eriImported;
        auto f12Result = F12::MP2_F12_Fixed(uhfWfn, configInts, f12Ints, f12Idxs);

        fmt::print("\n\n--- MP2 F12 Single Point Energy ---\n");
        fmt::print(" --- Fixed amplitude approximation used --- \n");
        fmt::print("  Molecule = {}\n", moleculeName[0] == 0 ? "Unknown" : moleculeName);
        fmt::print("  E(HF) = {} hartrees\n", uhfWfn.E_HF);
        fmt::print("  E(CABS singles correction) = {} hartrees\n", ECABSSingles);
        fmt::print("  E(MP2) = {} hartrees\n", f12Result.E_MP2);
        fmt::print("  E(F12 MP2) = {} hartrees\n", f12Result.E_F12_MP2);
        fmt::print("  E(Total) = {} hartrees\n\n", uhfWfn.E_HF + ECABSSingles + f12Result.E_MP2 +
                                                   f12Result.E_F12_MP2);

        fmt::print(fg(fmt::color::dark_green),
                   "\nPress any key to continue..\n");
        std::cin.get();
    }
    //
    void CCQMC_MP2_F12_H20() {
        Molecule water = Molecule::GenerateWaterMolecule();
        CCQMC_MP2_F12("sto-3g", "cc-pVDZ", water,
                      "Water (H2O)");
    }
    void CCQMC_MP2_F12_LiH() {
        Molecule lithiumHydride = Molecule::GenerateLithiumHydrideMolecule();
        CCQMC_MP2_F12("sto-3g", "cc-pVDZ", lithiumHydride,
                      "Lithium Hydride (LiH)");
    }
    //
    void CCQMC_MP2_F12(const char *obsName, const char *cabsName,
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
        fmt::print(fg(fmt::color::white), "  Molecule = {}\n", moleculeName);
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

        fmt::print(fg(fmt::color::dark_green),
                   "\nPress any key to continue..\n");
        std::cin.get();
    }
    //
}
