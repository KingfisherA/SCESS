#include <iostream>
#include <string>

#include "Resources/AssetFolder.h"
#include <fmt/color.h>
#include <fmt/format.h>
#include <libint2/basis.h>

#include "Dahlia/Utilities/Memory.h"
#include "Quantum/Basis/BasisSet.h"
#include "Quantum/Basis/CGF.h"
#include "Quantum/CC/CC.h"
#include "Quantum/CI/CI.h"
#include "Quantum/DFT/DFT.h"
#include "Quantum/F12/F12.h"
#include "Quantum/F12_QMC/F12_QMC.h"
#include "Quantum/HF/HF.h"
#include "Quantum/Integrals/IntegralTransforms.h"
#include "Quantum/Integrals/IntegralsF12.h"
#include "Quantum/Integrals/IntegralsLibInt.h"
#include "Quantum/Molecule.h"
#include "Quantum/QMC/QMC.h"
#include "Tensor/TensorOperations.h"
#include "Utilities/Utilities.h"

//
//
namespace SCESS::Tests {
//
//
// Reference energies (hartrees) for the built-in geometries. Unless noted they
// are from PySCF 2.14; small differences (< 1e-6) are expected for water as
// PySCF stores the STO-3G coefficients to fewer significant figures.
namespace Reference {
// H2 / STO-3G
constexpr double H2_RHF = -1.1175058851;
constexpr double H2_MP2 = -0.0125038087;
constexpr double H2_CCSD = -0.0193419373;
constexpr double H2_FCI = -1.1368478224;
// H2O / STO-3G
constexpr double H2O_RHF = -74.9629832288;
constexpr double H2O_MP2 = -0.0355225667;
constexpr double H2O_CCSD = -0.0494026324;
constexpr double H2O_T = -0.0000673499;
constexpr double H2O_FCI = -75.0125022529;
// LDA (Slater + VWN5) on a converged grid
constexpr double H2O_LDA = -74.7319604693;
// He / def2-TZVP
constexpr double He_RHF = -2.8598954257;
constexpr double He_FCI = -2.8929575257;
// H2 / STO-3G with the cc-pVDZ CABS. From the independent NumPy MP2-F12
// implementation in Scripts/verification/mp2f12_port.py (geminal integrals
// built from PySCF range-separated ERIs).
constexpr double H2_CABS_SINGLES = -0.0099723;
constexpr double H2_MP2_F12_FIXED = -0.022091052;
constexpr double H2_MP2_F12_OPT = -0.021554970;
// H2 / STO-3G CCSD-F12b with the cc-pVDZ CABS. From PySCF's RCCSD with the
// geminal couplings added in Scripts/verification/ccsdf12b_check.py
constexpr double H2_CCSD_F12B_CCSD = -0.0158468216;
constexpr double H2_CCSD_F12B_COUPLING = 0.0007516386;
// H2O / STO-3G CCSD-F12b with the cc-pVDZ CABS, also from ccsdf12b_check.py
constexpr double H2O_CCSD_F12B_CCSD = -0.0431035721;
constexpr double H2O_CCSD_F12B_COUPLING = 0.0009008077;
constexpr double H2O_CCSD_F12B_F12 = -0.3009608300;
} // namespace Reference
//
static int gFailures = 0;
//
// Compare an energy against its reference and report the result
static bool CheckEnergy(const char *label, double value, double reference,
                        double tolerance) {
  const double diff = std::abs(value - reference);
  const bool pass = diff <= tolerance;
  fmt::print(pass ? fg(fmt::color::green) : fg(fmt::color::crimson),
             "  [{}] {}: {:.10f} (reference {:.10f}, |diff| = {:.2e}, "
             "tolerance {:.0e})\n",
             pass ? "PASS" : "FAIL", label, value, reference, diff, tolerance);
  if (!pass)
    gFailures++;
  return pass;
}
//
void ProfileMolInts() {
  Molecule hydrogen = Molecule::GenerateWaterMolecule();

  libint2::BasisSet basisSet =
      BasisLibrary::LoadLibIntBasisSet("sto-3g", hydrogen.GetAtomsLibInt());

  libint2::Engine intEngine =
      libint2::Engine(libint2::Operator::coulomb, basisSet.max_nprim(),
                      static_cast<int>(basisSet.max_l()));
  long nInts = 1000000 / (basisSet.size() * basisSet.size() * basisSet.size() *
                          basisSet.size());
  auto tp1 = std::chrono::high_resolution_clock::now();
  for (long i = 0; i < nInts; i++)
    for (int s1 = 0; s1 < basisSet.size(); s1++)
      for (int s2 = 0; s2 < basisSet.size(); s2++)
        for (int s3 = 0; s3 < basisSet.size(); s3++)
          for (int s4 = 0; s4 < basisSet.size(); s4++)
            intEngine.compute(basisSet[s1], basisSet[s2], basisSet[s3],
                              basisSet[s4]);
  auto tp2 = std::chrono::high_resolution_clock::now();

  long timeDifference =
      std::chrono::duration_cast<std::chrono::microseconds>(tp2 - tp1).count();
  std::cout << "Time to evaluate 1E6 molecular integrals: " << timeDifference
            << '\n';
}
//
void TestMolecularInts() {
  Molecule hydrogenMolecule = Molecule::GenerateHydrogenMolecule();
  BasisSet basisSet("sto-3g", hydrogenMolecule);

  Integrals::MolInts ints = Integrals::EvaluateMolecularIntegrals(
      basisSet, hydrogenMolecule.GetAtoms());

  const PGF &pgf1 = basisSet.GetCGFs()[0].Get(0);

  std::vector<double> Bx, By, Bz;
  double primIntegral = Integrals::ERI_PGF(pgf1, pgf1, pgf1, pgf1, Bx, By, Bz);

  double val = boost::math::tgamma_lower(0.5, 1E-8);

  double boys = Integrals::Utils::BoysFunc(0, 0);

  std::cout << "--- Hydrogen molecular integrals ---\n";
  std::cout << "--- < ij | 1/r12 | kl > ---\n";
  Utilities::PrintTensor(ints.eri);

  std::cin.get();
}
//
void TestDFT() {
  std::cout << "--- Performing DFT Calculation on water ---\n";

  constexpr const char *BASIS_SET_NAME = "sto-3g";

  Molecule waterMolecule = Molecule::GenerateWaterMolecule();
  BasisSet basisSet(BASIS_SET_NAME, waterMolecule);

  std::stringstream console;
  DFT::KS_Wfn ksWavefunction =
      DFT::KS_DFT(waterMolecule, basisSet, Quadrature::MolGrid::MEDIUM_RES,
                  nullptr, &console);
  std::cout << console.str() << '\n';
  CheckEnergy("E(LDA)", ksWavefunction.E_KS, Reference::H2O_LDA, 1e-5);

  std::cout << "Press any key to continue..\n";
  std::cin.get();
}
//
void TestCoupledCluster() {
  std::cout << "---Performing CCSD Calculation on Hydrogen---\n";

  constexpr const char *BASIS_SET_NAME = "sto-3g";

  Molecule hydrogenMolecule = Molecule::GenerateHydrogenMolecule();
  BasisSet basisSet(BASIS_SET_NAME, hydrogenMolecule);

  Integrals::MolInts mints;

  std::stringstream ss;
  HF::HF_Wfn rhfWfn =
      HF::RHF(hydrogenMolecule, basisSet, 1e-10, 100, &mints, nullptr, &ss);
  std::cout << ss.str();
  ss.clear();
  ss.str("");

  CC::CCSD_Wfn ccsdWfn = CC::CCSD(rhfWfn, mints, nullptr, &ss);

  std::cout << ss.str() << '\n';

  std::cout << "\n\nE(CCSD Correlation) = " << ccsdWfn.E_CCSD << '\n';
  std::cout << "E(Total) = " << ccsdWfn.E_CCSD + rhfWfn.E_HF << '\n';
  // CCSD is exact for two electrons
  CheckEnergy("E(CCSD correlation)", ccsdWfn.E_CCSD, Reference::H2_CCSD, 1e-6);
  std::cout << "Press any key to continue..\n";
  std::cin.get();
}
//
void TestCCSD_T() {
  std::cout << "---Performing CCSD(T) Calculation on Water---\n";

  constexpr const char *BASIS_SET_NAME = "sto-3g";

  Molecule waterMolecule = Molecule::GenerateWaterMolecule();
  BasisSet basisSet(BASIS_SET_NAME, waterMolecule);

  Integrals::MolInts mints;

  std::stringstream ss;
  HF::HF_Wfn rhfWfn =
      HF::RHF(waterMolecule, basisSet, 1e-10, 100, &mints, nullptr, &ss);
  std::cout << ss.str();
  ss.clear();
  ss.str("");

  CC::CCSD_Wfn ccsdWfn = CC::CCSD_T(rhfWfn, mints, nullptr, &ss);

  std::cout << ss.str() << '\n';

  std::cout << "\n\nE(CCSD(T) Correlation) = " << ccsdWfn.E_CCSD + ccsdWfn.E_T
            << '\n';
  std::cout << "E(Total) = " << ccsdWfn.E_CCSD + ccsdWfn.E_T + rhfWfn.E_HF
            << '\n';
  CheckEnergy("E(CCSD correlation)", ccsdWfn.E_CCSD, Reference::H2O_CCSD,
              1e-6);
  CheckEnergy("E(T)", ccsdWfn.E_T, Reference::H2O_T, 1e-7);
  std::cout << "Press any key to continue..\n";
  std::cin.get();
}
//
void TestFullCI() {
  std::cout << "---Performing Full CI Calculation on Hydrogen ---\n";

  constexpr const char *BASIS_SET_NAME = "sto-3g";

  Molecule hydrogenMolecule = Molecule::GenerateHydrogenMolecule();
  BasisSet basisSet(BASIS_SET_NAME, hydrogenMolecule);

  Integrals::MolInts mints;
  std::stringstream ss;
  HF::HF_Wfn rhfWfn =
      HF::RHF(hydrogenMolecule, basisSet, 1e-05, 100, &mints, nullptr, &ss);
  HF::MP2_RHF_ClosedShellAO(rhfWfn, mints.eri, nullptr, &ss);

  std::cout << ss.str();
  ss.clear();
  ss.str("");

  CI::CI_Wfn ciWfn =
      CI::CI(rhfWfn, mints, hydrogenMolecule.CalculateVnn(), -1, nullptr, &ss);
  std::cout << ss.str();
  ss.clear();
  ss.str("");
  CheckEnergy("E(FCI)", ciWfn.E_CI(0), Reference::H2_FCI, 1e-8);

  std::cout << "Press any key to continue..\n";
  std::cin.get();
}
//
void TestFCIQMC() {
  std::cout << "Performing FCI calculation on hydrogen diatomic\n";

  constexpr const char *BASIS_SET_NAME = "sto-3g";
  Molecule hydrogenMolecule = Molecule::GenerateHydrogenMolecule();
  BasisSet basisSet(BASIS_SET_NAME, hydrogenMolecule);

  Integrals::MolInts mints;
  std::stringstream ss;
  HF::HF_Wfn rhfWfn =
      HF::RHF(hydrogenMolecule, basisSet, 1e-05, 100, &mints, nullptr, &ss);

  std::cout << ss.str() << '\n';

  double energy = QMC::FCI::FCIQMC(
      rhfWfn, mints, hydrogenMolecule.CalculateVnn(), nullptr, &ss);
  std::cout << ss.str() << '\n';
  // Stochastic; only a loose check
  CheckEnergy("E(FCIQMC)", energy, Reference::H2_FCI, 5e-3);

  std::cout << "Press any key to continue..\n";
  std::cin.get();
}
//
void TestCCQMC() {
  std::cout << "Performing CC calculation on hydrogen diatomic\n";

  constexpr const char *BASIS_SET_NAME = "sto-3g";
  Molecule hydrogenMolecule = Molecule::GenerateHydrogenMolecule();
  BasisSet basisSet(BASIS_SET_NAME, hydrogenMolecule);

  Integrals::MolInts mints;
  std::stringstream ss;
  HF::HF_Wfn rhfWfn =
      HF::RHF(hydrogenMolecule, basisSet, 1e-05, 100, &mints, nullptr, &ss);

  std::cout << ss.str() << '\n';

  double energy = QMC::CC::CCQMC(rhfWfn, mints);

  std::cout << "Final CC energy: " << energy << " hartrees\n";
  // Stochastic; only a loose check. CC is exact for two electrons
  CheckEnergy("E(CCQMC)", energy, Reference::H2_FCI, 5e-3);
  std::cout << "Press any key to continue..\n";
  std::cin.get();
}
//
void TestHelium() {
  constexpr const char *BASIS_SET_NAME = "def2-TZVP";
  std::cout << "Performing calculation on helium in a" << BASIS_SET_NAME
            << " basis set\n";

  Molecule heliumMolecule;
  heliumMolecule.AddAtom(Atom("He", 2, {0, 0, 0}));
  BasisSet basisSet(BASIS_SET_NAME, heliumMolecule);

  Integrals::MolInts mints;
  std::stringstream ss;
  HF::HF_Wfn rhfWfn =
      HF::RHF(heliumMolecule, basisSet, 1e-10, 100, &mints, nullptr, &ss);

  std::cout << ss.str() << '\n';
  ss.str("");
  ss.clear();

  CI::CI_Wfn ciWfn = CI::CI(rhfWfn, mints, 0.0, -1, nullptr, &ss);

  std::cout << ss.str() << '\n';
  CheckEnergy("E(RHF)", rhfWfn.E_HF, Reference::He_RHF, 1e-6);
  CheckEnergy("E(FCI)", ciWfn.E_CI(0), Reference::He_FCI, 1e-6);

  std::cout << "Press any key to continue..\n";
  std::cin.get();
}
//
void TestF12Integrals() {
  constexpr const char *BASIS_SET_NAME = "sto-3g";
  constexpr const char *CABS_BASIS_SET_NAME = "cc-pVDZ";
  Molecule hydrogenMolecule = Molecule::GenerateHydrogenMolecule();

  const F12::F12Calculation calc = F12::PrepareF12Calculation(
      BASIS_SET_NAME, CABS_BASIS_SET_NAME, hydrogenMolecule, nullptr, nullptr,
      true);
  const F12BasisSet &f12BasisSet = calc.basisSet;
  const F12::F12MOIntegrals &ints = calc.ints;

  std::cout << std::setprecision(15);

  std::cout << "\n--- F12 Integral Tensors ---\n";
  std::cout << "--- OBS (Ordinary Basis Set) ---\n";
  std::cout << "  * Name = " << BASIS_SET_NAME << '\n';
  std::cout << "  * Indicies = [" << f12BasisSet.GetHFIndices().first << " ... "
            << f12BasisSet.GetHFIndices().second << "]\n";
  std::cout << "--- CABS (Complementary Auxiliary Basis Set) ---\n";
  std::cout << "  * Name = " << CABS_BASIS_SET_NAME << '\n';
  std::cout << "  * Indicies = [" << f12BasisSet.GetCABSIndices().first
            << " ... " << f12BasisSet.GetCABSIndices().second << "]\n";

  std::cout << "--All are in physicists' notation, i, j occupied and p, q over "
               "the OBS + CABS--\n";

  std::cout << "------------------------------------\n";
  std::cout << "      <pq | 1/r12 | rs> (OBS)\n";
  std::cout << "------------------------------------\n";
  Utilities::PrintTensor(ints.eriOBS);

  std::cout << "\n\n------------------------------------\n";
  std::cout << "      <ij | 1/r12 | pq>\n";
  std::cout << "------------------------------------\n";
  Utilities::PrintTensor(ints.eri);

  std::cout << "\n\n------------------------------------\n";
  std::cout << "      <ij | f(12) | pq>\n";
  std::cout << "------------------------------------\n";
  Utilities::PrintTensor(ints.f12);

  std::cout << "\n\n------------------------------------\n";
  std::cout << "      <ij | f(12)^2 | pq>\n";
  std::cout << "------------------------------------\n";
  Utilities::PrintTensor(ints.f12Sq);

  std::cout << "\n\n------------------------------------\n";
  std::cout << "      <ij | |Del f(12)|^2 | pq>\n";
  std::cout << "------------------------------------\n";
  Utilities::PrintTensor(ints.f12DelSq);

  std::cout << "\n\n\n------------------------------------\n";
  std::cout << "      <ij | f(12)/r12 | pq>\n";
  std::cout << "------------------------------------\n";
  Utilities::PrintTensor(ints.f12Repulsion);

  std::cout << "Press any key to continue..\n";
  std::cin.get();
}
//
void TestMP2_F12_Optimise() {
  Molecule hydrogenMolecule = Molecule::GenerateHydrogenMolecule();
  const F12::F12Calculation calc =
      F12::PrepareF12Calculation("sto-3g", "cc-pVDZ", hydrogenMolecule);
  auto f12Result = F12::MP2_F12_Optimise(calc.hfWfn, calc.ints, calc.idxs);

  std::cout << "\n\n--- MP2 F12 Single Point Energy ---\n";
  std::cout << " --- Geminal amplitudes are optimised --- \n";
  std::cout << "  Molecule = Hydrogen (H2)\n";
  std::cout << "  E(HF) = " << calc.hfWfn.E_HF << " hartrees\n";
  std::cout << "  E(CABS singles correction) = " << calc.E_CABSSingles
            << " hartrees\n";
  std::cout << "  E(MP2) = " << f12Result.E_MP2 << " hartrees\n";
  std::cout << "  E(F12 MP2) = " << f12Result.E_F12_MP2 << " hartrees\n";
  std::cout << "  E(Total) = "
            << calc.hfWfn.E_HF + calc.E_CABSSingles + f12Result.E_MP2 +
                   f12Result.E_F12_MP2
            << " hartrees\n\n";
  CheckEnergy("E(MP2)", f12Result.E_MP2, Reference::H2_MP2, 1e-8);
  CheckEnergy("E(F12 MP2, optimised)", f12Result.E_F12_MP2,
              Reference::H2_MP2_F12_OPT, 1e-8);

  std::cout << "Press any key to continue..\n";
  std::cin.get();
}
//
void TestMP2_F12_Fixed() {
  Molecule hydrogenMolecule = Molecule::GenerateHydrogenMolecule();
  std::stringstream console;
  const F12::F12Calculation calc = F12::PrepareF12Calculation(
      "sto-3g", "cc-pVDZ", hydrogenMolecule, nullptr, &console);
  std::cout << console.str();
  auto f12Result = F12::MP2_F12_Fixed(calc.hfWfn, calc.ints, calc.idxs);

  std::cout << "\n\n--- MP2 F12 Single Point Energy ---\n";
  std::cout << " --- Fixed amplitude approximation used --- \n";
  std::cout << "  Molecule = Hydrogen (H2)\n";
  std::cout << "  E(HF) = " << calc.hfWfn.E_HF << " hartrees\n";
  std::cout << "  E(CABS singles correction) = " << calc.E_CABSSingles
            << " hartrees\n";
  std::cout << "  E(MP2) = " << f12Result.E_MP2 << " hartrees\n";
  std::cout << "  E(F12 MP2) = " << f12Result.E_F12_MP2 << " hartrees\n";
  std::cout << "  E(Total) = "
            << calc.hfWfn.E_HF + calc.E_CABSSingles + f12Result.E_MP2 +
                   f12Result.E_F12_MP2
            << " hartrees\n\n";
  CheckEnergy("E(HF)", calc.hfWfn.E_HF, Reference::H2_RHF, 1e-8);
  CheckEnergy("E(CABS singles)", calc.E_CABSSingles,
              Reference::H2_CABS_SINGLES, 1e-7);
  CheckEnergy("E(MP2)", f12Result.E_MP2, Reference::H2_MP2, 1e-8);
  CheckEnergy("E(F12 MP2, fixed)", f12Result.E_F12_MP2,
              Reference::H2_MP2_F12_FIXED, 1e-8);

  fmt::print(fg(fmt::color::dark_green), "\nPress any key to continue..\n");
  std::cin.get();
}
//
void TestCCQMC_MP2_F12() {
  Molecule hydrogenMolecule = Molecule::GenerateHydrogenMolecule();
  const F12::F12Calculation calc = F12::PrepareF12Calculation(
      "sto-3g", "cc-pVDZ", hydrogenMolecule, nullptr, nullptr, true);

  double ECCQMC_F12 =
      F12_QMC::CCSD_QMC_MP2_F12(calc.hfWfn, calc.ints, calc.idxs);
  double ECorrelation = ECCQMC_F12 - calc.hfWfn.E_HF;

  fmt::print(fg(fmt::color::white),
             "\n\n--- CCSD QMC F12 Single Point Energy ---\n");
  fmt::print(fg(fmt::color::white),
             " --- Fixed amplitude approximation used --- \n");
  fmt::print(fg(fmt::color::white), "  Molecule = Hydrogen (H2)\n");
  fmt::print(fg(fmt::color::white), "  E(HF) = {:.8} hartrees\n",
             calc.hfWfn.E_HF);
  fmt::print(fg(fmt::color::white),
             "  E(CABS singles correction) = {:.8} hartrees\n",
             calc.E_CABSSingles);
  fmt::print(fg(fmt::color::white), "  E(CCQMC_F12) = {:.8} hartrees\n",
             ECorrelation);
  fmt::print(fg(fmt::color::green), "  E(Total) = {:.8} hartrees\n\n",
             calc.hfWfn.E_HF + calc.E_CABSSingles + ECorrelation);

  fmt::print(fg(fmt::color::dark_green), "\nPress any key to continue..\n");
  std::cin.get();
}
//
void TestUHF() {
  Molecule hydrogenMolecule = Molecule::GenerateHydrogenMolecule();

  auto &hAtom = hydrogenMolecule.GetAtoms().back();
  auto hPos = hAtom.GetPos();
  hPos.x() = 3.0;
  hAtom.SetPos(hPos);

  auto h2Atoms = hydrogenMolecule.GetAtomsLibInt();

  libint2::BasisSet basisSet = BasisLibrary::LoadLibIntBasisSet("cc-pVDZ", h2Atoms);

  Integrals::MolInts mints =
      Integrals::LibInt::EvaluateMolecularIntegrals(basisSet, h2Atoms);

  HF::HF_Wfn uhfWfn = HF::UHF(mints, hydrogenMolecule.GetNoElectrons() / 2,
                              hydrogenMolecule.GetNoElectrons() / 2,
                              hydrogenMolecule.CalculateVnn());

  fmt::print("E(UHF) = {} hartrees\n", uhfWfn.E_HF);
}
//
int RunRegressionTests() {
  gFailures = 0;

  auto RunRHF = [](const Molecule &mol, const BasisSet &basis,
                   Integrals::MolInts &mints) {
    std::stringstream ss;
    return HF::RHF(mol, basis, 1e-10, 100, &mints, nullptr, &ss);
  };

  fmt::print(fmt::emphasis::bold, "\n--- H2 / STO-3G ---\n");
  {
    Molecule mol = Molecule::GenerateHydrogenMolecule();
    BasisSet basis("sto-3g", mol);
    Integrals::MolInts mints;
    HF::HF_Wfn rhf = RunRHF(mol, basis, mints);

    CheckEnergy("E(RHF)", rhf.E_HF, Reference::H2_RHF, 1e-8);
    CheckEnergy("E(MP2 correlation)",
                HF::MP2_RHF_ClosedShellAO(rhf, mints.eri), Reference::H2_MP2,
                1e-8);
    CheckEnergy("E(CCSD correlation)", CC::CCSD_DIIS(rhf, mints).E_CCSD,
                Reference::H2_CCSD, 1e-6);
    CheckEnergy("E(FCI)", CI::CI(rhf, mints, mol.CalculateVnn()).E_CI(0),
                Reference::H2_FCI, 1e-8);
    // Truncated CI at doubles is exact for two electrons
    CheckEnergy("E(CISD)", CI::CI(rhf, mints, mol.CalculateVnn(), 2).E_CI(0),
                Reference::H2_FCI, 1e-8);
  }

  fmt::print(fmt::emphasis::bold, "\n--- H2O / STO-3G ---\n");
  {
    Molecule mol = Molecule::GenerateWaterMolecule();
    BasisSet basis("sto-3g", mol);
    Integrals::MolInts mints;
    HF::HF_Wfn rhf = RunRHF(mol, basis, mints);

    CheckEnergy("E(RHF)", rhf.E_HF, Reference::H2O_RHF, 1e-6);
    CheckEnergy("E(MP2 correlation)",
                HF::MP2_RHF_ClosedShellAO(rhf, mints.eri), Reference::H2O_MP2,
                1e-7);
    CC::CCSD_Wfn ccsd = CC::CCSD_T(rhf, mints);
    CheckEnergy("E(CCSD correlation)", ccsd.E_CCSD, Reference::H2O_CCSD,
                1e-6);
    CheckEnergy("E(T)", ccsd.E_T, Reference::H2O_T, 1e-7);
    CheckEnergy("E(FCI)", CI::CI(rhf, mints, mol.CalculateVnn()).E_CI(0),
                Reference::H2O_FCI, 1e-6);

    DFT::KS_Wfn ks =
        DFT::KS_DFT(mol, basis, Quadrature::MolGrid::MEDIUM_RES, nullptr);
    CheckEnergy("E(LDA)", ks.E_KS, Reference::H2O_LDA, 1e-5);
  }

  fmt::print(fmt::emphasis::bold, "\n--- H2 / STO-3G MP2-F12 (cc-pVDZ CABS) ---\n");
  {
    Molecule mol = Molecule::GenerateHydrogenMolecule();
    const F12::F12Calculation calc =
        F12::PrepareF12Calculation("sto-3g", "cc-pVDZ", mol);
    CheckEnergy("E(CABS singles)", calc.E_CABSSingles,
                Reference::H2_CABS_SINGLES, 1e-7);
    auto RunF12 = [&](bool optimise) {
      return optimise ? F12::MP2_F12_Optimise(calc.hfWfn, calc.ints, calc.idxs)
                      : F12::MP2_F12_Fixed(calc.hfWfn, calc.ints, calc.idxs);
    };

    F12::MP2F12Result fixed = RunF12(false);
    CheckEnergy("E(MP2)", fixed.E_MP2, Reference::H2_MP2, 1e-8);
    CheckEnergy("E(F12 MP2, fixed)", fixed.E_F12_MP2,
                Reference::H2_MP2_F12_FIXED, 1e-8);
    CheckEnergy("E(F12 MP2, optimised)", RunF12(true).E_F12_MP2,
                Reference::H2_MP2_F12_OPT, 1e-8);
  }

  fmt::print(fmt::emphasis::bold,
             "\n--- H2 / STO-3G CCSD-F12b (cc-pVDZ CABS) ---\n");
  {
    Molecule mol = Molecule::GenerateHydrogenMolecule();
    const F12::F12Calculation calc =
        F12::PrepareF12Calculation("sto-3g", "cc-pVDZ", mol);
    const F12::CCSD_F12bResult result = F12::CCSD_F12b(calc);
    CheckEnergy("E(CCSD, with geminal coupling)", result.E_CCSD,
                Reference::H2_CCSD_F12B_CCSD, 1e-8);
    CheckEnergy("E(geminal coupling)", result.E_Coupling,
                Reference::H2_CCSD_F12B_COUPLING, 1e-8);
    CheckEnergy("E(F12 geminals)", result.E_F12, Reference::H2_MP2_F12_FIXED,
                1e-8);
  }

  fmt::print(fmt::emphasis::bold,
             "\n--- H2O / STO-3G CCSD-F12b (cc-pVDZ CABS) ---\n");
  {
    // PySCF uses fewer significant figures for the STO-3G coefficients, which
    // shifts these by up to 1e-7
    Molecule mol = Molecule::GenerateWaterMolecule();
    const F12::F12Calculation calc =
        F12::PrepareF12Calculation("sto-3g", "cc-pVDZ", mol);
    const F12::CCSD_F12bResult result = F12::CCSD_F12b(calc);
    CheckEnergy("E(CCSD, with geminal coupling)", result.E_CCSD,
                Reference::H2O_CCSD_F12B_CCSD, 1e-6);
    CheckEnergy("E(geminal coupling)", result.E_Coupling,
                Reference::H2O_CCSD_F12B_COUPLING, 1e-7);
    CheckEnergy("E(F12 geminals)", result.E_F12, Reference::H2O_CCSD_F12B_F12,
                2e-7);
  }

  fmt::print(fmt::emphasis::bold,
             "\n--- Closed-shell (spatial orbital) CCSD ---\n");
  {
    Molecule h2 = Molecule::GenerateHydrogenMolecule();
    BasisSet h2Basis("sto-3g", h2);
    Integrals::MolInts h2Ints;
    HF::HF_Wfn h2RHF = RunRHF(h2, h2Basis, h2Ints);
    CheckEnergy("H2 E(CCSD correlation)",
                CC::RCCSD_DIIS(h2RHF, h2Ints).E_CCSD, Reference::H2_CCSD,
                1e-8);

    Molecule water = Molecule::GenerateWaterMolecule();
    BasisSet waterBasis("sto-3g", water);
    Integrals::MolInts waterInts;
    HF::HF_Wfn waterRHF = RunRHF(water, waterBasis, waterInts);
    const double E_RCCSD = CC::RCCSD_DIIS(waterRHF, waterInts).E_CCSD;
    CheckEnergy("H2O E(CCSD correlation)", E_RCCSD, Reference::H2O_CCSD, 1e-7);
    // Same integrals and orbitals, so the spin-orbital CCSD must agree to its
    // convergence threshold
    CheckEnergy("H2O E(CCSD), spatial vs spin orbitals", E_RCCSD,
                CC::CCSD_DIIS(waterRHF, waterInts).E_CCSD, 1e-6);
  }

  fmt::print(fmt::emphasis::bold, "\n--- Basis set library ---\n");
  {
    auto CheckName = [](const char *label, const std::string &value,
                        const std::string &reference) {
      const bool pass = value == reference;
      fmt::print(pass ? fg(fmt::color::green) : fg(fmt::color::crimson),
                 "  [{}] {}: {} (reference {})\n", pass ? "PASS" : "FAIL",
                 label, value, reference);
      if (!pass)
        gFailures++;
    };
    auto Found = [](const char *name) -> std::string {
      const BasisSetInfo *info = BasisLibrary::Find(name);
      return info ? info->name : "<not found>";
    };
    CheckName("Find(6-31g**)", Found("6-31g**"), "6-31G(d,p)");
    CheckName("Find(CC-PVDZ)", Found("CC-PVDZ"), "cc-pVDZ");
    CheckName("DefaultCABS(cc-pVTZ-F12)",
              BasisLibrary::DefaultCABS("cc-pVTZ-F12"), "cc-pVTZ-F12-OptRI");
    CheckName("DefaultCABS(cc-pVDZ)", BasisLibrary::DefaultCABS("cc-pVDZ"),
              "aug-cc-pVDZ-OptRI");
    CheckName("DefaultCABS(sto-3g)", BasisLibrary::DefaultCABS("sto-3g"), "");
  }

  if (gFailures == 0)
    fmt::print(fg(fmt::color::green), "\nAll regression tests passed\n");
  else
    fmt::print(fg(fmt::color::crimson), "\n{} regression test(s) failed\n",
               gFailures);

  return gFailures;
}
//
//
} // namespace SCESS::Tests
