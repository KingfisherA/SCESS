#include<iostream>
#include<string>

#include "Resources/AssetFolder.h"
#define DATADIR SCESS::Res::AssetsFolder::Get().GetResourcePath("Basis").string()
#include <libint2/basis.h>

#include "GraphicsEngine/Utilities/Memory.h"
#include "Quantum/Molecule.h"
#include "Tensor/TensorOperations.h"
#include "Tensor/TensorEinsum.h"
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
    void ProfileMolInts() {
        Molecule hydrogen = Molecule::GenerateWaterMolecule();

        libint2::BasisSet basisSet("sto-3g", hydrogen.GetAtomsLibInt());

        libint2::Engine intEngine = libint2::Engine(libint2::Operator::coulomb,
                                                    basisSet.max_nprim(),
                                                    static_cast<int>(basisSet.max_l()));
        long nInts = 1000000 / (basisSet.size() * basisSet.size() * basisSet.size() * basisSet.size());
        auto tp1 = std::chrono::high_resolution_clock::now();
        for (long i = 0; i < nInts; i++)
            for (int s1 = 0; s1 < basisSet.size(); s1++)
                for (int s2 = 0; s2 < basisSet.size(); s2++)
                    for (int s3 = 0; s3 < basisSet.size(); s3++)
                        for (int s4 = 0; s4 < basisSet.size(); s4++)
                            intEngine.compute(basisSet[s1], basisSet[s2], basisSet[s3], basisSet[s4]);
        auto tp2 = std::chrono::high_resolution_clock::now();

        long timeDifference = std::chrono::duration_cast<std::chrono::microseconds>(tp2 - tp1).count();
        std::cout << "Time to evaluate 1E6 molecular integrals: " << timeDifference << '\n';
    }

    void TestMolecularInts() {
        Molecule hydrogenMolecule = Molecule::GenerateHydrogenMolecule();
        BasisSet basisSet("sto-3g", hydrogenMolecule);

        Integrals::MolInts ints = Integrals::EvaluateMolecularIntegrals(basisSet, hydrogenMolecule.GetAtoms());

        const Gaussian &pgf1 = basisSet.GetCGFs()[0].Get(0);

        std::vector<double> Bx, By, Bz;
        double primIntegral = Integrals::TwoElectronRepulsion(pgf1, pgf1, pgf1, pgf1, Bx, By, Bz);

        double val = boost::math::tgamma_lower(0.5, 1E-8);

        double boys = Integrals::Utils::BoysFunc(0, 0);

        std::cout << "--- Hydrogen molecular integrals ---\n";
        std::cout << "--- < ij | 1/r12 | kl > ---\n";
        Utilities::PrintTensor(ints.eri);

        std::cin.get();
    }

    void TestDFT() {
        std::cout << "--- Performing DFT Calculation on water ---\n";

        constexpr const char *BASIS_SET_NAME = "sto-3g";

        Molecule waterMolecule = Molecule::GenerateWaterMolecule();
        BasisSet basisSet(BASIS_SET_NAME, waterMolecule);

        std::stringstream console;
        uRef<KSWavefunction> ksWavefunction = DensityFunctional::KohnShamDFT(waterMolecule, basisSet,
                                                                             MolecularGrid::MEDIUM_RES, nullptr,
                                                                             &console);
        std::cout << console.str() << '\n';

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

        double energy = QMC::FCIQMC(*rhfResult,
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

        double energy = QMC::CCQMC(*rhfResult, hydrogenMolecule.CalculateNuclearRepulsion());

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
        Integrals::MolInts configInts = Integrals::F12::EvaluateMolecularIntegrals(f12BasisSet, h2Atoms);
        Integrals::F12::F12Integrals f12Ints = Integrals::F12::EvaluateF12Integrals(f12BasisSet);

        uRef<UHFWavefunction> uhfWfn = HartreeFock::UHF(
            Integrals::F12::ExtractHFIntegrals(configInts, 0, 2),
            1,
            1,
            hydrogenMolecule.CalculateNuclearRepulsion(),
            std::numeric_limits<double>::epsilon());
        Integrals::Transforms::TransformF12Basis(f12BasisSet,
                                                 uhfWfn->DensityAlpha,
                                                 uhfWfn->CAlpha,
                                                 configInts,
                                                 f12Ints);

        std::cout << std::setprecision(15);

        std::cout << "\n--- F12 Integral Tensors ---\n";
        std::cout << "--- OBS (Ordinary Basis Set) ---\n";
        std::cout << "  * Name = " << BASIS_SET_NAME << '\n';
        std::cout << "  * Indicies = [" << f12BasisSet.GetHFIndices().first << " ... " << f12BasisSet.GetHFIndices().
                second << "]\n";
        std::cout << "--- CABS (Complementary Auxiliary Basis Set) ---\n";
        std::cout << "  * Name = " << CABS_BASIS_SET_NAME << '\n';
        std::cout << "  * Indicies = [" << f12BasisSet.GetCABSIndices().first
                << " ... " << f12BasisSet.GetCABSIndices().second << "]\n";

        std::cout << "--All are in physicists' notation--\n";

        std::cout << "------------------------------------\n";
        std::cout << "      <ij | 1/r12 | kl>\n";
        std::cout << "------------------------------------\n";
        Utilities::PrintTensor(configInts.eri);

        std::cout << "\n\n------------------------------------\n";
        std::cout << "      <ij | f(12) | kl>\n";
        std::cout << "------------------------------------\n";
        Utilities::PrintTensor(f12Ints.f12);

        std::cout << "\n\n------------------------------------\n";
        std::cout << "      <ij | f(12)^2 | kl>\n";
        std::cout << "------------------------------------\n";
        Utilities::PrintTensor(f12Ints.f12Sq);

        std::cout << "\n\n------------------------------------\n";
        std::cout << "      <ij | |Del 1/r12|^2 | kl>\n";
        std::cout << "------------------------------------\n";
        Utilities::PrintTensor(f12Ints.f12DelSq);

        std::cout << "\n\n\n------------------------------------\n";
        std::cout << "      <ij | f(12)/r12 | kl>\n";
        std::cout << "------------------------------------\n";
        Utilities::PrintTensor(f12Ints.f12Repulsion);

        std::cout << "Press any key to continue..\n";
        std::cin.get();
    }

    void TestCCSD_F12() {
        constexpr const char *BASIS_SET_NAME = "sto-3g";
        constexpr const char *CABS_BASIS_SET_NAME = "cc-pVDZ";
        Molecule hydrogenMolecule = Molecule::GenerateHydrogenMolecule();

        auto h2Atoms = hydrogenMolecule.GetAtomsLibInt();

        F12BasisSet f12BasisSet(BASIS_SET_NAME, CABS_BASIS_SET_NAME, h2Atoms);
        Integrals::MolInts configInts = Integrals::F12::EvaluateMolecularIntegrals(f12BasisSet, h2Atoms);
        Integrals::F12::F12Integrals f12Ints = Integrals::F12::EvaluateF12Integrals(f12BasisSet);

        uRef<UHFWavefunction> uhfWfn = HartreeFock::UHF(configInts,
                                                        1,
                                                        1,
                                                        hydrogenMolecule.CalculateNuclearRepulsion(),
                                                        std::numeric_limits<double>::epsilon());
        double ECABSSingles = -0.009972302139;
        Integrals::Transforms::TransformF12Basis(f12BasisSet,
                                                 uhfWfn->DensityAlpha,
                                                 uhfWfn->CAlpha,
                                                 configInts,
                                                 f12Ints);

        double ccsdF12Energy = F12::CCSD::CCSD_F12(*uhfWfn, f12Ints);

        std::cout << "Press any key to continue..\n";
        std::cin.get();
    }

    void TestMP2_F12_Optimise() {
        constexpr const char *BASIS_SET_NAME = "sto-3g";
        constexpr const char *CABS_BASIS_SET_NAME = "cc-pVDZ";
        Molecule hydrogenMolecule = Molecule::GenerateHydrogenMolecule();

        auto h2Atoms = hydrogenMolecule.GetAtomsLibInt();

        F12BasisSet f12BasisSet(BASIS_SET_NAME, CABS_BASIS_SET_NAME, h2Atoms);
        Integrals::MolInts configInts = Integrals::F12::EvaluateMolecularIntegrals(f12BasisSet, h2Atoms);
        Integrals::F12::F12Integrals f12Ints = Integrals::F12::EvaluateF12Integrals(f12BasisSet);

        uRef<UHFWavefunction> uhfWfn =
                HartreeFock::UHF(
                    Integrals::F12::ExtractHFIntegrals(configInts, 0, 2),
                    1,
                    1,
                    hydrogenMolecule.CalculateNuclearRepulsion(),
                    std::numeric_limits<double>::epsilon());

        Eigen::MatrixXd F;
        Integrals::Transforms::TransformF12Basis(f12BasisSet,
                                                 uhfWfn->DensityAlpha,
                                                 uhfWfn->CAlpha,
                                                 configInts,
                                                 f12Ints,
                                                 true,
                                                 &F);

        //Need to correct the reference energy with a perturbation by singly exciting into the CABS virtuals.
        //Correct the fact that there is a large basis set error in the non-correlated HF energy
        double ECABSSingles = F12::EvaluateSinglesCorrectionEnergy(MatrixToTensor(F));

        //configInts.eri = eriImported;
        auto f12Result = F12::MP2::MP2_F12_Optimise(*uhfWfn, configInts, f12Ints);

        std::cout << "\n\n--- MP2 F12 Single Point Energy ---\n";
        std::cout << " --- Geminal amplitudes are optimised --- \n";
        std::cout << "  Molecule = Hydrogen (H2)\n";
        //Should be -1.117505880245 hartrees
        std::cout << "  E(HF) = " << uhfWfn->Energy << " hartrees\n";
        //Should be -0.009972302139 hartrees
        std::cout << "  E(CABS singles correction) = " << ECABSSingles << " hartrees\n";
        //Should be -0.012503809116 hartrees
        std::cout << "  E(MP2) = " << f12Result.E_MP2 << " hartrees\n";
        //Should be -0.021009866243 hartrees
        std::cout << "  E(F12 MP2) = " << f12Result.E_F12_MP2 << " hartrees\n";
        //Should be -1.160991857744 hartrees
        std::cout << "  E(Total) = " <<
                uhfWfn->Energy +
                ECABSSingles +
                f12Result.E_MP2 +
                f12Result.E_F12_MP2 << " hartrees\n\n";

        std::cout << "Press any key to continue..\n";
        std::cin.get();
    }

    void TestMP2_F12_Fixed() {
        constexpr const char *BASIS_SET_NAME = "sto-3g";
        constexpr const char *CABS_BASIS_SET_NAME = "cc-pVDZ";
        Molecule hydrogenMolecule = Molecule::GenerateHydrogenMolecule();

        auto h2Atoms = hydrogenMolecule.GetAtomsLibInt();

        F12BasisSet f12BasisSet(BASIS_SET_NAME, CABS_BASIS_SET_NAME, h2Atoms);
        Integrals::MolInts configInts = Integrals::F12::EvaluateMolecularIntegrals(f12BasisSet, h2Atoms);
        Integrals::F12::F12Integrals f12Ints = Integrals::F12::EvaluateF12Integrals(f12BasisSet);

        uRef<UHFWavefunction> uhfWfn =
                HartreeFock::UHF(
                    Integrals::F12::ExtractHFIntegrals(configInts, 0, 2),
                    1,
                    1,
                    hydrogenMolecule.CalculateNuclearRepulsion(),
                    std::numeric_limits<double>::epsilon());

        Eigen::MatrixXd F;
        Integrals::Transforms::TransformF12Basis(f12BasisSet,
                                                 uhfWfn->DensityAlpha,
                                                 uhfWfn->CAlpha,
                                                 configInts,
                                                 f12Ints,
                                                 true,
                                                 &F);

        //Need to correct the reference energy with a perturbation by singly exciting into the CABS virtuals.
        //Correct the fact that there is a large basis set error in the non-correlated HF energy
        double ECABSSingles = F12::EvaluateSinglesCorrectionEnergy(MatrixToTensor(F));

        //configInts.eri = eriImported;
        auto f12Result = F12::MP2::MP2_F12_Fixed(*uhfWfn, configInts, f12Ints);

        std::cout << "\n\n--- MP2 F12 Single Point Energy ---\n";
        std::cout << " --- Fixed amplitude approximation used --- \n";
        std::cout << "  Molecule = Hydrogen (H2)\n";
        //Should be -1.117505880245 hartrees
        std::cout << "  E(HF) = " << uhfWfn->Energy << " hartrees\n";
        //Should be -0.009972302139 hartrees
        std::cout << "  E(CABS singles correction) = " << ECABSSingles << " hartrees\n";
        //Should be -0.012503809116 hartrees
        std::cout << "  E(MP2) = " << f12Result.E_MP2 << " hartrees\n";
        //Should be -0.021009866243 hartrees
        std::cout << "  E(F12 MP2) = " << f12Result.E_F12_MP2 << " hartrees\n";
        //Should be -1.160991857744 hartrees
        std::cout << "  E(Total) = " <<
                uhfWfn->Energy +
                ECABSSingles +
                f12Result.E_MP2 +
                f12Result.E_F12_MP2 << " hartrees\n\n";

        std::cout << "Press any key to continue..\n";
        std::cin.get();
    }

    void TestCCQMC_F12() {
        constexpr const char *BASIS_SET_NAME = "sto-3g";
        constexpr const char *CABS_BASIS_SET_NAME = "cc-pVDZ";
        Molecule hydrogenMolecule = Molecule::GenerateHydrogenMolecule();

        auto h2Atoms = hydrogenMolecule.GetAtomsLibInt();

        F12BasisSet f12BasisSet(BASIS_SET_NAME, CABS_BASIS_SET_NAME, h2Atoms);
        Integrals::MolInts configInts = Integrals::F12::EvaluateMolecularIntegrals(f12BasisSet, h2Atoms);
        Integrals::F12::F12Integrals f12Ints = Integrals::F12::EvaluateF12Integrals(f12BasisSet);

        uRef<UHFWavefunction> uhfWfn =
                HartreeFock::UHF(
                    Integrals::F12::ExtractHFIntegrals(configInts, 0, 2),
                    1,
                    1,
                    hydrogenMolecule.CalculateNuclearRepulsion(),
                    std::numeric_limits<double>::epsilon());

        Eigen::MatrixXd F;
        Integrals::Transforms::TransformF12Basis(f12BasisSet,
                                                 uhfWfn->DensityAlpha,
                                                 uhfWfn->CAlpha,
                                                 configInts,
                                                 f12Ints,
                                                 true,
                                                 &F);

        //Need to correct the reference energy with a perturbation by singly exciting into the CABS virtuals.
        //Correct the fact that there is a large basis set error in the non-correlated HF energy

        double ECABSSingles = F12::EvaluateSinglesCorrectionEnergy(MatrixToTensor(F));

        double ECCQMC_F12 = F12::QMC::CCQMC_MP2_F12(*uhfWfn, configInts, f12Ints);
        double ECorrelation = ECCQMC_F12 - uhfWfn->Energy;

        std::cout << "\n\n--- CCQMC F12 Single Point Energy ---\n";
        std::cout << " --- Fixed amplitude approximation used --- \n";
        std::cout << "  Molecule = Hydrogen (H2)\n";
        //Should be -1.117505880245 hartrees
        std::cout << "  E(HF) = " << uhfWfn->Energy << " hartrees\n";
        //Should be -0.009972302139 hartrees
        std::cout << "  E(CABS singles correction) = " << ECABSSingles << " hartrees\n";
        //Should be -0.012503809116 hartrees
        std::cout << "  E(CCQMC_F12) = " << ECorrelation << " hartrees\n";
        //Should be -0.021009866243 hartrees
        //Should be -1.160991857744 hartrees
        std::cout << "  E(Total) = " <<
                uhfWfn->Energy +
                ECABSSingles +
                ECorrelation << " hartrees\n\n";

        std::cout << "\nPress any key to continue..\n";
        std::cin.get();
    }
}
