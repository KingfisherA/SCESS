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
#include "Quantum/Integrals/IntegralTransforms.h"
#include "Quantum/Integrals/IntegralsF12.h"
#include "Quantum/Integrals/IntegralsLibInt.h"
#include "Quantum/Basis/BasisSet.h"
#include "Quantum/F12/F12.h"

namespace SCESS::Tests {
    void TestDFT() {
        std::cout << "--- Performing DFT Calculation on water ---\n";

        constexpr const char *BASIS_SET_NAME = "sto-3g";

        Molecule waterMolecule = Molecule::GenerateWaterMolecule();
        BasisSet basisSet(BASIS_SET_NAME, waterMolecule);

        uRef<KSWavefunction> ksWavefunction = DensityFunctional::KohnShamDFT(waterMolecule, basisSet,
                                                                             MolecularGrid::MEDIUM_RES);

        std::cout << "Press any key to continue..\n";
        std::cin.get();
    }

    void TestCoupledCluster() {
        std::cout << "---Performing CCSD Calculation on Hydrogen---\n";

        constexpr const char *BASIS_SET_NAME = "sto-3g";

        Molecule hydrogenMolecule = Molecule::GenerateHydrogenMolecule();
        BasisSet basisSet(BASIS_SET_NAME, hydrogenMolecule);

        std::stringstream ss;
        uRef<RHFWavefunction> rhfResult = HartreeFock::RHF(hydrogenMolecule,
                                                           basisSet,
                                                           1e-05,
                                                           100,
                                                           nullptr,
                                                           &ss,
                                                           BASIS_SET_NAME);
        std::cout << ss.str();
        ss.clear();
        ss.str("");

        uRef<CCSDWavefunction> ccsd = CoupledCluster::CCSD(*rhfResult, nullptr, &ss);

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

        Molecule hydrogenMolecule = Molecule::GenerateHydrogenMolecule();
        BasisSet basisSet(BASIS_SET_NAME, hydrogenMolecule);

        std::stringstream ss;
        uRef<RHFWavefunction> rhfResult = HartreeFock::RHF(hydrogenMolecule,
                                                           basisSet,
                                                           1e-05,
                                                           100,
                                                           nullptr,
                                                           &ss,
                                                           BASIS_SET_NAME);
        HartreeFock::MP2ClosedShell(*rhfResult, nullptr, &ss);

        std::cout << ss.str();
        ss.clear();
        ss.str("");


        uRef<CIWavefunction> ciResult = ConfigurationInteraction::FCI(
            *rhfResult,
            hydrogenMolecule.CalculateNuclearRepulsion(),
            nullptr,
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
        Molecule hydrogenMolecule = Molecule::GenerateHydrogenMolecule();
        BasisSet basisSet(BASIS_SET_NAME, hydrogenMolecule);

        std::stringstream ss;
        uRef<RHFWavefunction> rhfResult = HartreeFock::RHF(hydrogenMolecule,
                                                           basisSet,
                                                           1e-05,
                                                           100,
                                                           nullptr,
                                                           &ss,
                                                           BASIS_SET_NAME);

        std::cout << ss.str() << '\n';

        double energy = QuantumMonteCarlo::FCIQMC(*rhfResult,
                                                  hydrogenMolecule.CalculateNuclearRepulsion(),
                                                  nullptr,
                                                  &ss);


        std::cout << "Press any key to continue..\n";
        std::cin.get();
    }

    void TestCCQMC() {
        std::cout << "Performing CCQMC calculation on hydrogen diatomic\n";

        constexpr const char *BASIS_SET_NAME = "sto-3g";
        Molecule hydrogenMolecule = Molecule::GenerateHydrogenMolecule();
        BasisSet basisSet(BASIS_SET_NAME, hydrogenMolecule);

        std::stringstream ss;
        uRef<RHFWavefunction> rhfResult = HartreeFock::RHF(hydrogenMolecule,
                                                           basisSet,
                                                           1e-08,
                                                           100,
                                                           nullptr,
                                                           &ss,
                                                           BASIS_SET_NAME);

        std::cout << ss.str() << '\n';

        double energy = QuantumMonteCarlo::CCQMC(*rhfResult, hydrogenMolecule.CalculateNuclearRepulsion());

        std::cout << "Final CCQMC energy: " << energy << " hartrees\n";
        std::cout << "Press any key to continue..\n";
        std::cin.get();
    }

    void TestHelium() {
        constexpr const char *BASIS_SET_NAME = "def2-TZVP";
        std::cout << "Performing calculation on helium in a" << BASIS_SET_NAME << " basis set\n";

        Molecule heliumMolecule;
        heliumMolecule.AddAtom(Atom("He", 2, {0, 0, 0}));
        BasisSet basisSet(BASIS_SET_NAME, heliumMolecule);

        std::stringstream ss;
        uRef<RHFWavefunction> rhfWfn = HartreeFock::RHF(heliumMolecule,
                                                        basisSet,
                                                        1e-08,
                                                        100,
                                                        nullptr,
                                                        &ss,
                                                        BASIS_SET_NAME);

        std::cout << ss.str() << '\n';
        ss.str("");
        ss.clear();

        auto ciWfn = ConfigurationInteraction::FCI(*rhfWfn, 0.0, nullptr, &ss);

        std::cout << ss.str() << '\n';

        std::cout << "Press any key to continue..\n";
        std::cin.get();
    }

    void TestF12Integrals() {
        constexpr const char *BASIS_SET_NAME = "sto-3g";
        constexpr const char *CABS_BASIS_SET_NAME = "cc-pVDZ";
        Molecule hydrogenMolecule = Molecule::GenerateHydrogenMolecule();

        auto h2Atoms = hydrogenMolecule.GetAtomsLibInt();

        F12BasisSet f12BasisSet(BASIS_SET_NAME, CABS_BASIS_SET_NAME, h2Atoms);
        Integrals::ConfigInts configInts = Integrals::F12::CalculateConfiguration(f12BasisSet, h2Atoms);
        Integrals::F12::F12Integrals f12Ints = Integrals::F12::EvaluateF12Integrals(f12BasisSet);

        uRef<UHFWavefunction> hfWfn;
        Integrals::Transforms::TransformF12Basis(f12BasisSet, configInts, f12Ints, &hfWfn);

        std::cout << std::setprecision(15);

        std::cout << "\n--- F12 Integral Tensors ---\n";
        std::cout << "--- OBS (Ordinary Basis Set) ---\n";
        std::cout << "  * Name = " << BASIS_SET_NAME << '\n';
        std::cout << "  * Indicies = [" << f12BasisSet.GetHFIndices().first << " ... " << f12BasisSet.GetHFIndices().
                second << "]\n";
        std::cout << "--- CABS (Complementary Auxiliary Basis Set) ---\n";
        std::cout << "  * Name = " << CABS_BASIS_SET_NAME << '\n';
        std::cout << "  * Indicies = [" << f12BasisSet.GetCABSIndices().first << " ... " << f12BasisSet.GetCABSIndices()
                .second << "]\n";

        std::cout << "--All are in physicists' notation--\n";

        std::cout << "------------------------------------\n";
        std::cout << "      <ij | 1/r12 | kl>\n";
        std::cout << "------------------------------------\n";
        Utilities::PrintTensor(Integrals::Transforms::ConvertChemistsToPhysicists(f12Ints.eri));

        std::cout << "\n\n------------------------------------\n";
        std::cout << "      <ij | f(12) | kl>\n";
        std::cout << "------------------------------------\n";
        Utilities::PrintTensor(Integrals::Transforms::ConvertChemistsToPhysicists(f12Ints.f12));

        std::cout << "\n\n------------------------------------\n";
        std::cout << "      <ij | f(12)^2 | kl>\n";
        std::cout << "------------------------------------\n";
        Utilities::PrintTensor(Integrals::Transforms::ConvertChemistsToPhysicists(f12Ints.f12Sq));

        std::cout << "\n\n------------------------------------\n";
        std::cout << "      <ij | |Del 1/r12|^2 | kl>\n";
        std::cout << "------------------------------------\n";
        Utilities::PrintTensor(Integrals::Transforms::ConvertChemistsToPhysicists(f12Ints.f12DelSq));

        std::cout << "\n\n\n------------------------------------\n";
        std::cout << "      <ij | f(12)/r12 | kl>\n";
        std::cout << "------------------------------------\n";
        Utilities::PrintTensor(Integrals::Transforms::ConvertChemistsToPhysicists(f12Ints.f12Repulsion));

        std::cout << "Press any key to continue..\n";
        std::cin.get();
    }

    void TestCCSD_F12() {
        constexpr const char *BASIS_SET_NAME = "sto-3g";
        constexpr const char *CABS_BASIS_SET_NAME = "cc-pVDZ";
        Molecule hydrogenMolecule = Molecule::GenerateHydrogenMolecule();

        auto h2Atoms = hydrogenMolecule.GetAtomsLibInt();

        F12BasisSet f12BasisSet(BASIS_SET_NAME, CABS_BASIS_SET_NAME, h2Atoms);
        Integrals::ConfigInts configInts = Integrals::F12::CalculateConfiguration(f12BasisSet, h2Atoms);
        Integrals::F12::F12Integrals f12Ints = Integrals::F12::EvaluateF12Integrals(f12BasisSet);

        uRef<UHFWavefunction> uhfWfn;
        Integrals::Transforms::TransformF12Basis(f12BasisSet, configInts, f12Ints, &uhfWfn);

        double ccsdF12Energy = F12::CCSD::CCSD_F12(*uhfWfn, f12Ints);

        std::cout << "Press any key to continue..\n";
        std::cin.get();
    }

    void TestMP2_F12() {
        Integrals::F12::F12Integrals f12IntsImported =
                Integrals::F12::ImportF12IntegralsMRCC(
                    Res::AssetsFolder::Get().GetAssetsPathStr() + "/F12 Data/MRCC_F12_Ints.txt", 12);

        constexpr const char *BASIS_SET_NAME = "sto-3g";
        constexpr const char *CABS_BASIS_SET_NAME = "cc-pVDZ";
        Molecule hydrogenMolecule = Molecule::GenerateHydrogenMolecule();

        auto h2Atoms = hydrogenMolecule.GetAtomsLibInt();

        F12BasisSet f12BasisSet(BASIS_SET_NAME, CABS_BASIS_SET_NAME, h2Atoms);
        Integrals::ConfigInts configInts = Integrals::F12::CalculateConfiguration(f12BasisSet, h2Atoms);
        Integrals::F12::F12Integrals f12Ints = Integrals::F12::EvaluateF12Integrals(f12BasisSet);

        uRef<UHFWavefunction> uhfWfn;
        //Need to correct the reference energy with a pertubation. CABS singles.
        //Correct the fact that there is a large basis set error in the non-correlated HF energy
        double cabsSinglesCorrection = -0.009972302139;
        Integrals::Transforms::TransformF12Basis(f12BasisSet, configInts, f12Ints, &uhfWfn, nullptr);

        auto f12Result = F12::MP2::MP2_F12(*uhfWfn, configInts, f12IntsImported);

        std::cout << "\n\n--- MP2 F12 Single Point Energy ---\n";
        std::cout << "  Molecule = Hydrogen (H2)\n";
        //Should be -1.117505880245 hartrees
        std::cout << "  E(HF) = " << uhfWfn->Energy << " hartrees\n";
        //Should be -0.009972302139
        std::cout << "  E(CABS singles correction) = " << cabsSinglesCorrection << " hartrees\n";
        //Should be -0.012503809116
        std::cout << "  E(MP2) = " << f12Result.mp2Energy << " hartrees\n";
        //Should be -0.021009866243 hartrees
        std::cout << "  E(F12 MP2) = " << f12Result.f12MP2Energy << " hartrees\n";
        //Should be -1.160991857744 hartrees
        std::cout << "  E(Total) = " << uhfWfn->Energy + cabsSinglesCorrection + f12Result.mp2Energy + f12Result.
                f12MP2Energy << " hartrees\n\n";

        std::cout << "Press any key to continue..\n";
        std::cin.get();
    }
}
