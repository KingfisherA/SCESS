#include <string>
#include <fmt/format.h>
#include <fmt/color.h>
#include <unordered_map>
#include <boost/algorithm/string/case_conv.hpp>
//
#include "Driver.h"
#include "Tests.h"
//
#include "Quantum/Molecule.h"
#include "Quantum/Integrals/IntegralsLibInt.h"
//
#include <boost/algorithm/string/trim.hpp>

#include "Resources/AssetFolder.h"
//


//
//
void PrintHelpDialog();
//
int main(int argc, char *argv[]) {
    using namespace SCESS;
    fmt::print(fmt::emphasis::bold, "Schrodinger's Cat's Electronic System Solver; Terminal Edition\n"
               "    By Charlie Wheelwright (cspw2)"
               "\n    Meow!\n\n");

    if (argc - 1 == 0) {
        PrintHelpDialog();
        return 0;
    }

    std::vector<std::string> args;
    args.reserve(argc - 1);

    for (int i = 1; i < argc; i++)
        args.push_back(argv[i]);

    for (auto &arg: args) {
        if (arg == "--help") {
            PrintHelpDialog();
            return 0;
        }
        boost::algorithm::trim(arg);
    }

    if (args.size() < 3) {
        fmt::print(fg(fmt::color::crimson),
                   "SCESS requires at least 3 arguments however only {} were provided\n", args.size());
        return 1;
    }

    if (!Res::AssetsFolder::Get().Init()) {
        fmt::print(fg(fmt::color::crimson),
                   "Required assets directory could not be located\n");
        return 1;
    }

    enum class MethodEnum :int {
        RHF,
        UHF,
        DFT,
        MP2,
        CCSD,
        CCSD_T,
        FCI,
        FCIQMC,
        CCQMC,
        MP2_F12_Fixed,
        MP2_F12_Opt,
        CCSD_QMC_MP2_F12
    } selectedMethod;
    enum class MoleculeEnum :int {
        H2,
        LiH,
        NH3,
        H2O
    } selectedMol;
    std::unordered_map<std::string, MethodEnum> methodMap = {
        {"rhf", MethodEnum::RHF},
        {"uhf", MethodEnum::UHF},
        {"dft", MethodEnum::DFT},
        {"mp2", MethodEnum::MP2},
        {"ccsd", MethodEnum::CCSD},
        {"ccsd_t", MethodEnum::CCSD_T},
        {"fci", MethodEnum::FCI},
        {"fciqmc", MethodEnum::FCIQMC},
        {"ccqmc", MethodEnum::CCQMC},
        {"mp2_f12_fixed", MethodEnum::MP2_F12_Fixed},
        {"mp2_f12_opt", MethodEnum::MP2_F12_Opt},
        {"ccsd_qmc_mp2_f12", MethodEnum::CCSD_QMC_MP2_F12}
    };
    std::unordered_map<std::string, MoleculeEnum> moleculeMap = {
        {"H2", MoleculeEnum::H2},
        {"LiH", MoleculeEnum::LiH},
        {"NH3", MoleculeEnum::NH3},
        {"H2O", MoleculeEnum::H2O},
    };
    std::unordered_set<std::string> allowableAdditionalOps = {"CABS"};
    std::unordered_map<std::string, std::string> options;

    fmt::print(fmt::emphasis::bold, "\n\n--- Selected Launch Options ---\n");

    auto &methodStr = args[0];
    boost::algorithm::to_lower(args[0]);
    if (!methodMap.contains(methodStr)) {
        fmt::print(fg(fmt::color::red), "{} is not a valid theory method driver", methodStr);
        return 1;
    }
    selectedMethod = methodMap.at(methodStr);
    fmt::print("    - Method = ");
    fmt::print(fg(fmt::color::green), "{}\n", methodStr);

    auto &moleculeStr = args[1];
    if (!moleculeMap.contains(moleculeStr)) {
        fmt::print(fg(fmt::color::red), "{} is not an avaliable molecule", moleculeStr);
        return 1;
    }
    selectedMol = moleculeMap.at(moleculeStr);
    fmt::print("    - Molecule = ");
    fmt::print(fg(fmt::color::green), "{}\n", moleculeStr);

    auto basisSets = BasisSet::ListBasisSets();
    for (auto &bs: basisSets)
        basisSets.insert(boost::algorithm::to_lower_copy(bs));

    auto &selectedBasisSet = args[2];
    boost::algorithm::to_lower(selectedBasisSet);
    if (!basisSets.contains(selectedBasisSet)) {
        fmt::print(fg(fmt::color::red), "{} is not an avaliable basis set", selectedBasisSet);
        return 1;
    }
    fmt::print("    - Basis Set = ");
    fmt::print(fg(fmt::color::green), "{}\n", selectedBasisSet);

    for (int i = 3; i < args.size(); ++i) {
        std::string &op = args[i];

        auto markerStart = op.find("--");
        if (markerStart != std::string::npos && (i + 1) < args.size()) {
            op = op.substr(markerStart + 2, op.length() - (markerStart + 2));
            if (allowableAdditionalOps.contains(op))
                options.emplace(op, args[i + 1]);
            else {
                fmt::print(fg(fmt::color::crimson), "{} is not a valid additional option\n", op);
                return 1;
            }
            ++i;
        }
    }

    if (!options.empty()) {
        fmt::print(fmt::emphasis::bold, "--- Selected Additional Options ---\n");
        for (auto &[k1, v1]: options) {
            fmt::print("    - {} = ", k1);
            fmt::print(fg(fmt::color::green), "{}\n", v1);
        }
    }


    Integrals::LibInt::InitLibInt();

    Molecule mol;
    const char *moleculeName = nullptr;

    switch (selectedMol) {
        case MoleculeEnum::H2:
            mol = Molecule::GenerateHydrogenMolecule();
            moleculeName = "Hydrogen (H2)";
            break;
        case MoleculeEnum::LiH:
            mol = Molecule::GenerateLithiumHydrideMolecule();
            moleculeName = "Lithium Hydride (LiH)";
            break;
        case MoleculeEnum::H2O:
            mol = Molecule::GenerateWaterMolecule();
            moleculeName = "Water (H2O)";
            break;
        case MoleculeEnum::NH3:
            mol = Molecule::GenerateAmmoniaMolecule();
            moleculeName = "Ammonia (NH3)";
            break;
    }

    const char *cabsBasisSet = nullptr;
    if (selectedMethod == MethodEnum::MP2_F12_Fixed ||
        selectedMethod == MethodEnum::MP2_F12_Opt ||
        selectedMethod == MethodEnum::CCSD_QMC_MP2_F12) {
        if (!options.contains("CABS")) {
            fmt::print(fg(fmt::color::crimson), "F12 method selected but no CABS basis set was specified\n");
            return 1;
        }
        boost::algorithm::to_lower(options.at("CABS"));
        if (!basisSets.contains(options.at("CABS"))) {
            fmt::print(fg(fmt::color::crimson), "CABS basis set {} is invalid\n", options.at("CABS"));
            return 1;
        }
        cabsBasisSet = options.at("CABS").c_str();
    }

    fmt::print("\n\n");

    switch (selectedMethod) {
        case MethodEnum::RHF:
            Driver::RHF(selectedBasisSet.c_str(), mol, moleculeName);
            break;
        case MethodEnum::UHF:
            Driver::UHF(selectedBasisSet.c_str(), mol, moleculeName);
            break;
        case MethodEnum::DFT:
            Driver::DFT(selectedBasisSet.c_str(), mol, moleculeName);
            break;
        case MethodEnum::MP2:
            Driver::MP2(selectedBasisSet.c_str(), mol, moleculeName);
            break;
        case MethodEnum::CCSD:
            Driver::CCSD(selectedBasisSet.c_str(), mol, moleculeName);
            break;
        case MethodEnum::CCSD_T:
            Driver::CCSD_T(selectedBasisSet.c_str(), mol, moleculeName);
            break;
        case MethodEnum::FCI:
            Driver::FCI(selectedBasisSet.c_str(), mol, moleculeName);
            break;
        case MethodEnum::FCIQMC:
            Driver::FCIQMC(selectedBasisSet.c_str(), mol, moleculeName);
            break;
        case MethodEnum::CCQMC:
            Driver::CCQMC(selectedBasisSet.c_str(), mol, moleculeName);
            break;
        case MethodEnum::MP2_F12_Fixed:
            Driver::MP2_F12_Fixed(selectedBasisSet.c_str(), cabsBasisSet, mol, moleculeName);
            break;
        case MethodEnum::MP2_F12_Opt:
            Driver::MP2_F12_Optimise(selectedBasisSet.c_str(), cabsBasisSet, mol, moleculeName);
            break;
        case MethodEnum::CCSD_QMC_MP2_F12:
            Driver::CCSD_QMC_MP2_F12(selectedBasisSet.c_str(), cabsBasisSet, mol, moleculeName);
            break;
    }

    Integrals::LibInt::CloseLibInt();

    return 0;
}
//
void PrintHelpDialog() {
    fmt::print("SCESS ");
    fmt::print(fmt::emphasis::italic, "[method] [molecule] [basis set] [additional options]\n");
    fmt::print(
        "\n"
        "--- Available Theory Method Drivers ---\n");
    fmt::print(fg(fmt::color::green) | fmt::emphasis::bold, "    SCF\n");
    fmt::print(""
        "        - RHF = Restricted Hartree-Fock\n"
        "        - UHF = Unrestricted Hartree-Fock\n"
        "        - DFT = Kohn-Sham DFT (VWM + Slater-X)\n"
        "    \n");
    fmt::print(fg(fmt::color::teal) | fmt::emphasis::bold, "    MP\n");
    fmt::print(""
        "        - MP2 = Moller-Plesset Second Order \n"
        "\n");
    fmt::print(fg(fmt::color::violet) | fmt::emphasis::bold, "    Configuration Interaction\n");
    fmt::print(""
        "        - FCI = Full Configuration Interaction\n"
        "\n");
    fmt::print(fg(fmt::color::red) | fmt::emphasis::bold, "    Coupled-Cluster\n");
    fmt::print(""
        "        - CCSD   = Coupled Cluster Singles and Doubles\n"
        "        - CCSD_T = Coupled Cluster Singles and Doubles with perturbative triples\n"
        "\n");
    fmt::print(fg(fmt::color::gold) | fmt::emphasis::bold, "    Quantum Monte Carlo (QMC)\n");
    fmt::print(""
        "        - FCIQMC = Full-CI QMC\n"
        "        - CCQMC = Coupled-Cluster QMC\n"
        "\n");
    fmt::print(fg(fmt::color::silver) | fmt::emphasis::bold, "    F12\n");
    fmt::print(""
        "        - MP2_F12_Fixed = MP2 F12 theory with the fixed amplitude approximation \n"
        "        - MP2_F12_Opt   = MP2 F12 theory with the optimised amplitudes\n"
        "\n");
    fmt::print(fg(fmt::color::crimson) | fmt::emphasis::bold, "    QMC + F12\n");
    fmt::print(""
        "        - CCSD_QMC_MP2_F12 = CCSD-QMC with first-order T2' operator\n"
        "\n"
        "Note!\n"
        "    F12 theories require the --CABS option to be specified and set to one of the available basis sets.\n"
        "\n");
    fmt::print("--- Available Molecules ---\n"
        "    - H2  = Hydrogen \n"
        "    - LiH = Lithium Hydride\n"
        "    - H2O = Water\n"
        "    - NH3 = Ammonia\n"
        "\n");

    if (SCESS::Res::AssetsFolder::Get().Init()) {
        fmt::print("--- Available Basis Sets ---\n");
        auto basisSets = SCESS::BasisSet::ListBasisSets();
        for (const auto &basisSet: basisSets) {
            fmt::print("    - ");
            fmt::print(fg(fmt::color::silver) | fmt::emphasis::italic, "{}\n", basisSet);
        }
        fmt::print("\n");
    } else {
        fmt::print(fg(fmt::color::crimson),
                   "Required assets directory could not be located to list avaliable basis sets\n");
    }

    fmt::print("--- Additional Options ---\n"
        "    --CABS = Sets the complementary auxiliary basis set for a F12 theory. Can be set to any of the available basis sets\n"
        "\n");
}
//
//
