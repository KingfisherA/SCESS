#include<iostream>
#include<string>

#include"GraphicsEngine/Utilities/Memory.h"
#include "Chemistry/Quantum/Molecule.h"
#include "Chemistry/Quantum/HF/HartreeFock.h"
#include "Chemistry/Quantum/DFT/DFT.h"
#include "Chemistry/Quantum/CC/CoupledCluster.h"
#include "Chemistry/Quantum/CI/ConfigurationInteraction.h"
#include "Chemistry/Quantum/QMC/QMC.h"

namespace SCESS::Tests {
    void TestDFT() {
        std::cout << "--- Performing DFT Calculation on water ---\n";

        constexpr const char *BASIS_SET_NAME = "def2-svp.1";

        Molecule waterMolecule = GenerateWaterMolecule();
        waterMolecule.MakeOrbitals(BASIS_SET_NAME);

        DensityFunctional::DFT(waterMolecule);

        std::cout << "Press any key to continue..\n";
        std::cin.get();
    }

    void TestCoupledCluster() {
        std::cout << "---Performing Hartree Fock Calculation on Water---\n";

        constexpr const char *BASIS_SET_NAME = "sto-3g";

        Molecule waterMolecule = GenerateWaterMolecule();
        waterMolecule.MakeOrbitals(BASIS_SET_NAME);

        std::stringstream ss;
        uRef<RHFWavefunction> rhfResult = HartreeFock::RHF(waterMolecule, nullptr, nullptr, 1e-05, 100, &ss,
                                                           BASIS_SET_NAME);
        std::cout << ss.str();
        ss.clear();
        ss.str("");

        uRef<CCSDResult> result = CoupledCluster::CCSD(*rhfResult);

        std::cout << "Press any key to continue..\n";
        std::cin.get();
    }

    void TestFullCI() {
        std::cout << "---Performing Hartree Fock Calculation on Hydrogen ---\n";

        constexpr const char *BASIS_SET_NAME = "sto-3g";

        Molecule hydrogenMolecule = GenerateHydrogenMolecule();
        hydrogenMolecule.MakeOrbitals(BASIS_SET_NAME);

        std::stringstream ss;
        uRef<RHFWavefunction> rhfResult = HartreeFock::RHF(hydrogenMolecule, nullptr, nullptr, 1e-08, 100, &ss,
                                                           BASIS_SET_NAME);
        HartreeFock::MP2ClosedShell(*rhfResult, &ss);

        std::cout << ss.str();
        ss.clear();
        ss.str("");


        uRef<CIWavefunction> ciResult = ConfigurationInteraction::FCI(
            *rhfResult, hydrogenMolecule.CalculateNuclearRepulsion(),
            &ss);
        std::cout << ss.str();
        ss.clear();
        ss.str("");

        std::cout << "Press any key to continue..\n";
        std::cin.get();
    }

    void TestFCIQMC() {
        std::cout << "Performing calculation on hydrogen diatomic\n";

        constexpr const char *BASIS_SET_NAME = "sto-3g";
        Molecule hydrogenMolecule = GenerateHydrogenMolecule();
        hydrogenMolecule.MakeOrbitals(BASIS_SET_NAME);

        std::stringstream ss;
        uRef<RHFWavefunction> rhfResult = HartreeFock::RHF(hydrogenMolecule, nullptr, nullptr, 1e-08, 100, &ss,
                                                           BASIS_SET_NAME);

        std::cout << ss.str() << '\n';

        double energy = QuantumMonteCarlo::SimpleFCIQMC(*rhfResult, hydrogenMolecule.CalculateNuclearRepulsion(), &ss);


        std::cout << "Press any key to continue..\n";
        std::cin.get();
    }
}
