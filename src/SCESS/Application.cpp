#include <iostream>
#include<string>
//
#include "Tests.h"
#include "Driver.h"
#include "Quantum/Integrals/IntegralsLibInt.h"
#include "Resources/AssetFolder.h"

#include<fmt/format.h>
#include<fmt/color.h>
#include <cppitertools/combinations.hpp>

#include "Quantum/Molecule.h"
#include "Quantum/CI/CI.h"


void SelectCCQMC_MP2_F12();

//
//
int main(int argc, char *argv[]) {
    using namespace SCESS;
    fmt::print(fmt::emphasis::bold, "Schrodinger's Cat's Electronic System Solver; Terminal Edition\n"
               "    By Charlie Wheelwright (cspw2)"
               "\n    Meow!\n\n");
    if (!Res::AssetsFolder::Get().Init())
        return 1;

    Integrals::LibInt::InitLibInt();

    fmt::print(fmt::emphasis::italic, "\n\n\n"
               "--- Command Line Arguments ---\n"
               "\t- \"f12ints\" = Calculate a test set of f12 integrals and transforms the basis set\n"
               "\t- \"mp2f12-opt\" = Perform an MP2-F12 optimised amplitude calculation on a hydrogen molecule.\n"
               "\t- \"mp2f12-fixed\" = Perform an MP2-F12 fixed amplitude calculation on a hydrogen molecule.\n"
               "\t- \"ccsd\" = Perform a CCSD calculation on hydrogen\n"
               "\t- \"dft\" = Perform a test KS-DFT calculation on water with the LDA exchange and VWD correlation functional.\n"
               "\t- \"fci\" = Perform a test full configuration interaction calculation on hydrogen.\n"
               "\t- \"fciqmc\" = Perform a test Full-CI Quantum Monte Carlo procedure on hydrogen.\n"
               "\t- \"ccqmc-mp2-f12\" = Perform a CCSD Quantum Monte Carlo with a second-order geminal excitor operator T2' in the T expansion (new theory)."
               "\n\n\n");

    // Molecule ammonia = Molecule::GenerateAmmoniaMolecule();
    // BasisSet basisSet("sto-3g", ammonia);
    //
    // auto noSpinOrbs = 2 * basisSet.GetNoCGFs();
    // auto dets1 = CI::EnumerateDet(noSpinOrbs,
    //                               ammonia.GetNoElectrons() / 2,
    //                               ammonia.GetNoElectrons() / 2,
    //                               2);
    // auto dets2 = CI::EnumerateAllDet(10, 2);

    // volatile bool equal = true;
    //
    // for (auto d1: dets1) {
    //     bool hasFound = false;
    //     for (auto d2: dets2) {
    //         if (d1.GetAlphaOrbitals() == d2.GetAlphaOrbitals() && d1.GetBetaOrbitals() == d2.GetBetaOrbitals()) {
    //             hasFound = true;
    //             break;
    //         }
    //     }
    //
    //     if (!hasFound) {
    //         equal = false;
    //         break;
    //     }
    //}

    //Tests::TestUHF();
    //Tests::TestMP2_F12_Fixed();
    //Tests::TestFCIQMC();
    //Tests::TestCCQMC();
    //Tests::TestCCQMC_MP2_F12();
    //Tests::TestFullCI();
    //Tests::TestMolecularInts();
    //Tests::TestF12Integrals();
    //Tests::TestTensorEinsum();
    SelectCCQMC_MP2_F12();

    //Driver::MP2_F12_Fixed_LiH();
    //Driver::CCQMC_MP2_F12_LiH();

    if (argc > 1) {
        std::string arg1 = argv[1];
        if (arg1 == "f12ints")
            Tests::TestF12Integrals();
        else if (arg1 == "mp2f12-opt")
            Tests::TestMP2_F12_Optimise();
        else if (arg1 == "mp2f12-fixed")
            Tests::TestMP2_F12_Fixed();
        else if (arg1 == "ccsd")
            Tests::TestCoupledCluster();
        else if (arg1 == "dft")
            Tests::TestDFT();
        else if (arg1 == "fci")
            Tests::TestFullCI();
        else if (arg1 == "fciqmc")
            Tests::TestFCIQMC();
        else if (arg1 == "ccqmc-mp2-f12")
            SelectCCQMC_MP2_F12();
        else
            std::cout << "Command \"" << arg1 << "\" not recognised\n";
    } else {
        std::cout << "--- No command given ---\n";
    }

    Integrals::LibInt::CloseLibInt();

    return 0;
}
//
void SelectCCQMC_MP2_F12() {
    using namespace SCESS;

    fmt::print(fmt::emphasis::bold, "\nPlease select an option:\n");
    fmt::print(fmt::emphasis::italic,
               "    [0] Hydrogen Single Point\n"
               "    [1] Hydrogen Dissociation Curve\n"
               "    [2] Lithium Hydride Single Point \n"
               "    [3] Water Single Point\n"
               "    [4] Ammonia Single Point\n"
               "    [5] Methane Single Point\n"
               "    [6] Carbon Monoxide Single Point\n\n");

    auto GetNumericEntry = [](int maxValue) {
        int option;
        while (true) {
            fmt::print("[~] ");
            try {
                std::cin >> option;

                if (option <= maxValue)
                    return option;
            } catch (...) {
            }
            fmt::print(fg(fmt::color::red), "\nPlease select a valid option\n");
        }
    };

    int option = GetNumericEntry(6);

    std::string moleculeName;
    Molecule molecule;

    switch (option) {
        case 0: {
            molecule = Molecule::GenerateHydrogenMolecule();
            moleculeName = "Hydrogen (H2)";
            break;
        }
        case 1: return;
        case 2: {
            molecule = Molecule::GenerateLithiumHydrideMolecule();
            moleculeName = "Lithium Hydride (LiH)";
            break;
        }
        case 3: {
            molecule = Molecule::GenerateWaterMolecule();
            moleculeName = "Water (H2O)";
            break;
        }
        case 4: {
            molecule = Molecule::GenerateAmmoniaMolecule();
            moleculeName = "Ammonia (NH3)";
            break;
        }
        case 5: {
            molecule = Molecule::GenerateMethaneMolecule();
            moleculeName = "Methane (CH4)";
            break;
        }
        case 6: {
            molecule = Molecule::GenerateCOMolecule();
            moleculeName = "Carbon Monoxide (CO)";
            break;
        }
        default:
            return;
    }


    std::string obsName, cabsName;

    auto OrdinalToBasisSetName = [](int ordinal) {
        switch (ordinal) {
            case 0:
                return "sto-3g";
            case 1:
                return "cc-pVDZ";
            case 2:
                return "cc-pVTZ";
            default:
                return "";
        }
    };

    fmt::print(fmt::emphasis::bold, "\nPlease select a basis set for the OBS space:\n");
    fmt::print(fmt::emphasis::italic,
               "    [0] sto-3g\n"
               "    [1] cc-pVDZ\n"
               "    [2] cc-pVTZ \n\n");

    option = GetNumericEntry(2);
    obsName = OrdinalToBasisSetName(option);

    fmt::print(fmt::emphasis::bold, "\nPlease select a basis set for the CABS space:\n");
    fmt::print(fmt::emphasis::italic,
               "    [0] sto-3g\n"
               "    [1] cc-pVDZ\n"
               "    [2] cc-pVTZ \n\n");

    option = GetNumericEntry(2);
    cabsName = OrdinalToBasisSetName(option);

    Driver::CCQMC_MP2_F12(obsName.c_str(), cabsName.c_str(), molecule, moleculeName.c_str());
}
