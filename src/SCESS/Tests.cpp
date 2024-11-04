#include<iostream>
#include<string>

#include "GraphicsEngine/Utilities/Memory.h"
#include "Quantum/Molecule.h"
#include "Quantum/HF/HartreeFock.h"
#include "Quantum/DFT/DensityFunctional.h"
#include "Quantum/CC/CoupledCluster.h"
#include "Quantum/CI/ConfigurationInteraction.h"
#include "Quantum/QMC/QMC.h"
#include "Utilities/Utilities.h"
#include "Quantum/Integrals/IntegralsF12.h"

namespace SCESS::Tests {
    void TestDFT() {
        std::cout << "--- Performing DFT Calculation on water ---\n";

        constexpr const char *BASIS_SET_NAME = "sto-3g";

        Molecule waterMolecule = GenerateWaterMolecule();
        waterMolecule.MakeOrbitals(BASIS_SET_NAME);

        uRef<KSWavefunction> ksWavefunction = DensityFunctional::KohnShamDFT(waterMolecule, MolecularGrid::MEDIUM_RES);

        std::cout << "Press any key to continue..\n";
        std::cin.get();
    }

    void TestCoupledCluster() {
        std::cout << "---Performing CCSD Calculation on Hydrogen---\n";

        constexpr const char *BASIS_SET_NAME = "sto-3g";

        Molecule hydrogenMolecule = GenerateHydrogenMolecule();
        hydrogenMolecule.MakeOrbitals(BASIS_SET_NAME);

        std::stringstream ss;
        uRef<RHFWavefunction> rhfResult = HartreeFock::RHF(hydrogenMolecule, nullptr, nullptr, 1e-05, 100, &ss,
                                                           BASIS_SET_NAME);
        std::cout << ss.str();
        ss.clear();
        ss.str("");

        uRef<CCSDWavefunction> ccsd = CoupledCluster::CCSD(*rhfResult, &ss);

        std::cout << ss.str() << '\n';

        Utilities::PrintTensor(ccsd->t1Amps);
        Utilities::PrintTensor(ccsd->t2Amps);

        std::cout << "\n\nE(CCSD Correlation) = " << ccsd->ccsdCorrelationEnergy << '\n';
        std::cout << "E(Total) = " << ccsd->ccsdCorrelationEnergy + rhfResult->Energy << '\n';
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
        std::cout << "Performing FCIQMC calculation on hydrogen diatomic\n";

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

    void TestCCQMC() {
        std::cout << "Performing CCQMC calculation on hydrogen diatomic\n";

        constexpr const char *BASIS_SET_NAME = "sto-3g";
        Molecule hydrogenMolecule = GenerateHydrogenMolecule();
        hydrogenMolecule.MakeOrbitals(BASIS_SET_NAME);

        std::stringstream ss;
        uRef<RHFWavefunction> rhfResult = HartreeFock::RHF(hydrogenMolecule, nullptr, nullptr, 1e-08, 100, &ss,
                                                           BASIS_SET_NAME);

        std::cout << ss.str() << '\n';

        double energy = QuantumMonteCarlo::CCQMC(*rhfResult, hydrogenMolecule.CalculateNuclearRepulsion());


        std::cout << "Press any key to continue..\n";
        std::cin.get();
    }

    void TestHelium() {
        constexpr const char *BASIS_SET_NAME = "def2-TZVP";
        std::cout << "Performing calculation on helium in a" << BASIS_SET_NAME << " basis set\n";

        Molecule heliumMolecule;
        heliumMolecule.InsertAtom(Atom("He", 2, {0, 0, 0}));
        heliumMolecule.MakeOrbitals(BASIS_SET_NAME);

        std::stringstream ss;
        uRef<RHFWavefunction> rhfWfn = HartreeFock::RHF(heliumMolecule, nullptr, nullptr, 1e-08, 100, &ss,
                                                        BASIS_SET_NAME);

        std::cout << ss.str() << '\n';
        ss.str("");
        ss.clear();

        auto ciWfn = ConfigurationInteraction::FCI(*rhfWfn, 0.0, &ss);

        std::cout << ss.str() << '\n';

        std::cout << "Press any key to continue..\n";
        std::cin.get();
    }

    void TestIntegrals() {
        constexpr const char *BASIS_SET_NAME = "sto-3g";
        Molecule hydrogenMolecule = GenerateHydrogenMolecule();
        hydrogenMolecule.MakeOrbitals(BASIS_SET_NAME);

        TEITensor f12Ints = Integrals::F12::TwoElectronSlaterF12Tensor(hydrogenMolecule, 0.5);
        TEITensor f12RepulsionInts = Integrals::F12::TwoElectronSlaterF12RepulsionTensor(hydrogenMolecule, 0.5);

        std::cout << "(ij | r12 * exp(-gamma * r12) | kl) integral tensor:\n";
        Utilities::PrintTensor(f12Ints);

        std::cout << "(ij | r12 * exp(-gamma * r12) * 1/r12 | kl) integral tensor:\n";
        Utilities::PrintTensor(f12RepulsionInts);

        std::cout << "Press any key to continue..\n";
        std::cin.get();
    }
}
