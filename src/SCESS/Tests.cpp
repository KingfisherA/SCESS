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
        uRef<RHFWavefunction> rhfResult = HartreeFock::RHF(hydrogenMolecule,
                                                           nullptr,
                                                           nullptr,
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

        Molecule hydrogenMolecule = GenerateHydrogenMolecule();
        hydrogenMolecule.MakeOrbitals(BASIS_SET_NAME);

        std::stringstream ss;
        uRef<RHFWavefunction> rhfResult = HartreeFock::RHF(hydrogenMolecule,
                                                           nullptr,
                                                           nullptr,
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
        Molecule hydrogenMolecule = GenerateHydrogenMolecule();
        hydrogenMolecule.MakeOrbitals(BASIS_SET_NAME);

        std::stringstream ss;
        uRef<RHFWavefunction> rhfResult = HartreeFock::RHF(hydrogenMolecule,
                                                           nullptr,
                                                           nullptr,
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
        Molecule hydrogenMolecule = GenerateHydrogenMolecule();
        hydrogenMolecule.MakeOrbitals(BASIS_SET_NAME);

        std::stringstream ss;
        uRef<RHFWavefunction> rhfResult = HartreeFock::RHF(hydrogenMolecule,
                                                           nullptr,
                                                           nullptr,
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
        heliumMolecule.InsertAtom(Atom("He", 2, {0, 0, 0}));
        heliumMolecule.MakeOrbitals(BASIS_SET_NAME);

        std::stringstream ss;
        uRef<RHFWavefunction> rhfWfn = HartreeFock::RHF(heliumMolecule,
                                                        nullptr,
                                                        nullptr,
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

    void TestRepulsionIntegrals() {
        constexpr const char *BASIS_SET_NAME = "sto-3g";
        Molecule waterMolecule = GenerateWaterMolecule();
        waterMolecule.MakeOrbitals(BASIS_SET_NAME);

        libint2::BasisSet basisSet(BASIS_SET_NAME, waterMolecule.GetAtomsLibInt());

        TEITensor eriSCESS = Integrals::TwoElectronRepulsionTensorMultithread(waterMolecule);
        TEITensor eriLibInt = Integrals::LibInt::TwoElectronRepulsionTensor(basisSet);

        std::cout << "SCESS (ij | 1/r12 | kl) integral tensor:\n";
        Utilities::PrintTensor(eriSCESS);

        std::cout << "LibInt (ij | 1/r12 | kl) integral tensor:\n";
        Utilities::PrintTensor(eriLibInt);

        if (eriSCESS.size() == eriLibInt.size()) {
            auto compareDouble = [](double value1, double value2, uint8_t precision)-> bool {
                return std::abs(value1 - value2) < std::pow(10, -precision);
            };

            constexpr uint8_t PRECISION = 4;
            bool equal = true;
            for (size_t i = 0; i < eriSCESS.size(); i++)
                if (compareDouble(eriSCESS.data()[i], eriLibInt.data()[i], PRECISION)) {
                    equal = false;
                    break;
                }

            if (equal)
                std::cout << "Integral tensors are the same to " << PRECISION <<
                        " decimal places. Integrals calculated correctly.\n";
            else
                std::cout << "Integral tensors are not the same. There must be an error in the calculation code.\n";
        } else {
            std::cout << "ERI Tensors are not the same size. Integral comparison failed\n";
        }
        std::cout << "Press any key to continue..\n";
        std::cin.get();
    }

    void TestF12Integrals() {
        constexpr const char *BASIS_SET_NAME = "sto-3g";
        constexpr const char *CABS_BASIS_SET_NAME = "cc-pVDZ";
        Molecule hydrogenMolecule = GenerateHydrogenMolecule();
        hydrogenMolecule.MakeOrbitals(BASIS_SET_NAME);

        libint2::BasisSet hfBasisSet(BASIS_SET_NAME, hydrogenMolecule.GetAtomsLibInt());
        libint2::BasisSet cabsBasisSet(CABS_BASIS_SET_NAME, hydrogenMolecule.GetAtomsLibInt());

        Integrals::F12::LibIntCGF slaterGeminal = Integrals::F12::CalculateSlaterContraction(1.0);
        Integrals::F12::LibIntCGF slaterGeminal2 = Integrals::F12::CalculateSlaterSquaredContraction(slaterGeminal);

        //Calculates (ij | 1/r12 | kl)
        TEITensor repulsionInts = Integrals::LibInt::TwoElectronRepulsionTensor(hfBasisSet);
        //Calculates (ij | f(12) | kl)
        TEITensor f12Ints = Integrals::F12::TwoElectronCGFF12Tensor(hfBasisSet, slaterGeminal);
        //Calculates (ij | f(12)^2 | kl)
        TEITensor f12SquaredInts = Integrals::F12::TwoElectronCGFF12Tensor(hfBasisSet, slaterGeminal2);
        //Calculates (ij | |Del f(12)|^2 | kl)
        TEITensor f12DelInts = Integrals::F12::TwoElectronCGFDelF12SquaredTensor(hfBasisSet, slaterGeminal);
        //Calculates (ij | f(12)/r12 | kl)
        TEITensor f12RepulsionInts = Integrals::F12::TwoElectronCGFF12RepulsionTensor(hfBasisSet, slaterGeminal);

        uRef<UHFWavefunction> uhfWfn = HartreeFock::UHF(hydrogenMolecule, 1E-10);

        repulsionInts = Integrals::Transforms::ConvertChemistsToPhysicists(
            Integrals::Transforms::ToMolecularTEI(repulsionInts, uhfWfn->CAlpha));
        f12Ints = Integrals::Transforms::ConvertChemistsToPhysicists(
            Integrals::Transforms::ToMolecularTEI(f12Ints, uhfWfn->CAlpha));
        f12SquaredInts = Integrals::Transforms::ConvertChemistsToPhysicists(
            Integrals::Transforms::ToMolecularTEI(f12SquaredInts, uhfWfn->CAlpha));
        f12DelInts = Integrals::Transforms::ConvertChemistsToPhysicists(
            Integrals::Transforms::ToMolecularTEI(f12DelInts, uhfWfn->CAlpha));
        f12RepulsionInts = Integrals::Transforms::ConvertChemistsToPhysicists(
            Integrals::Transforms::ToMolecularTEI(f12RepulsionInts, uhfWfn->CAlpha));

        std::cout << std::setprecision(15);

        std::cout << "F12 Integral Tensors\n--All are in physicits' notation--\n";

        std::cout << "<ij | 1/r12 | kl>\n";
        Utilities::PrintTensor(repulsionInts);
        std::cout << "<ij | f(12) | kl>\n";
        Utilities::PrintTensor(f12Ints);
        std::cout << "<ij | f(12)^2 | kl>\n";
        Utilities::PrintTensor(f12SquaredInts);
        std::cout << "<ij | |Del 1/r12|^2 | kl>\n";
        Utilities::PrintTensor(f12DelInts);
        std::cout << "<ij | f(12)/r12 | kl>\n";
        Utilities::PrintTensor(f12RepulsionInts);

        std::cout << "Press any key to continue..\n";
        std::cin.get();
    }
}
