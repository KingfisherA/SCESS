#include<Quantum/QuantumPCH.h>

#include<fstream>
#include<map>

#include "Quantum/Molecule.h"
#include "Quantum/Basis/BasisSet.h"
#include "Quantum/Basis/Gaussian.h"
#include "PeriodicTable.h"
#include "Resources/AssetFolder.h"

#include "GraphicsEngine/Utilities/StringOperations.h"

//
//
namespace SCESS {
    using namespace GEng;

    static std::map<std::string, int> atomicSymbolToCharge = {
        {"H", 1},
        {"He", 2},
        {"Li", 3},
        {"Be", 4},
        {"B", 5},
        {"C", 6},
        {"N", 7},
        {"O", 8},
        {"F", 9},
        {"Ne", 10},
        {"Na", 11},
        {"Mg", 12},
        {"Al", 13},
        {"Si", 14},
        {"P", 15},
        {"S", 16},
        {"Cl", 17},
        {"Ar", 18},
        {"K", 19},
        {"Ca", 20},
        {"Sc", 21},
        {"Ti", 22},
        {"V", 23},
        {"Cr", 24},
        {"Mn", 25},
        {"Fe", 26},
        {"Co", 27},
        {"Ni", 28},
        {"Cu", 29},
        {"Zn", 30},
        {"Ga", 31},
        {"Ge", 32},
        {"As", 33},
        {"Se", 34},
        {"Br", 35},
        {"Kr", 36},
        {"Rb", 37},
        {"Sr", 38},
        {"Y", 39},
        {"Zr", 40},
        {"Nb", 41},
        {"Mo", 42},
        {"Tc", 43},
        {"Ru", 44},
        {"Rh", 45},
        {"Pd", 46},
        {"Ag", 47},
        {"Cd", 48},
        {"In", 49},
        {"Sn", 50},
        {"Sb", 51},
        {"Te", 52},
        {"I", 53},
        {"Xe", 54}
    };

    Atom::Atom(std::string atomSymbol, const Eigen::Vector3d &pos): mSymbol(atomSymbol), mPosition(pos),
                                                                    mCharge(atomicSymbolToCharge.at(atomSymbol)) {
    }

    void Atom::SetAtomSymbol(std::string s) {
        mCharge = atomicSymbolToCharge.at(s);
        mSymbol = s;
    }

    //
    //
    void Molecule::AddAtom(const Atom &a) {
        mNoElectrons += a.GetCharge();
        mAtoms.emplace_back(a);
    }

    std::vector<libint2::Atom> Molecule::GetAtomsLibInt() const {
        std::vector<libint2::Atom> atoms(mAtoms.size());

        for (uint32_t i = 0; i < mAtoms.size(); i++) {
            const SCESS::Atom &atm = mAtoms[i];
            atoms[i] = {
                atm.GetCharge(),
                atm.GetPos().x(), atm.GetPos().y(), atm.GetPos().z()
            };
        }

        return atoms;
    }

    void Molecule::ConstructBonding(double tolerance2) {
        std::map<double, std::vector<BondPair> > distanceMap;

        for (size_t i = 0; i < mAtoms.size(); i++) {
            for (size_t j = 0; j < mAtoms.size(); j++) {
                if (i == j) continue;
                double distance2 = (mAtoms[i].GetPos() - mAtoms[j].GetPos()).squaredNorm();

                auto it = distanceMap.find(distance2);
                if (it != distanceMap.end()) {
                    (*it).second.push_back({i, j});
                    continue;
                }
                distanceMap[distance2] = {{i, j}};
            }
        }


        if (distanceMap.size() == 0)
            return;

        for (auto &it: distanceMap) {
            if (it.first > distanceMap.begin()->first + tolerance2)
                continue;

            for (const BondPair &p: it.second) {
                mBonds.insert(p);
            }
        }
    }

    double Molecule::CalculateNuclearRepulsion() const {
        double Vnn = 0.0;
        for (size_t i = 0; i < mAtoms.size(); i++) {
            for (size_t j = i + 1; j < mAtoms.size(); j++) {
                double r = (mAtoms[i].GetPos() - mAtoms[j].GetPos()).norm();
                if (r > 0.0) Vnn += mAtoms[i].GetCharge() * mAtoms[j].GetCharge() / r;
            }
        }
        return Vnn;
    }

    Molecule Molecule::LoadFromXYZ(const std::string &filePath) {
        std::ifstream inStream(filePath);
        Molecule molecule;

        if (inStream.is_open()) {
            std::string line;
            while (std::getline(inStream, line)) {
                Atom a;
                size_t startIdx = line.find("  ", 0);
                if (startIdx == line.npos)
                    continue;
                std::vector<std::string> tokens = StringOperations::Split(line.substr(startIdx), ' ');
                a.SetAtomSymbol(tokens[0]);
                double x = std::stod(tokens[1]);
                double y = std::stod(tokens[2]);
                double z = std::stod(tokens[3]);
                a.SetPos({x, y, z});
                molecule.AddAtom(a);
            }
        }

        return molecule;
    }

    Molecule Molecule::GenerateWaterMolecule() {
        Atom H1("H", Eigen::Vector3d(0.7569520, -0.586349, 0.0) * ANGSTROM_2_BOHR);
        Atom H2("H", Eigen::Vector3d(-0.7569520, -0.586349, 0.0) * ANGSTROM_2_BOHR);
        Atom O("O", Eigen::Vector3d(0.0, 0.0, 0.0));

        Molecule water({H1, O, H2});


        water.ConstructBonding();

        return water;
    }

    Molecule Molecule::GenerateHydrogenMolecule() {
        constexpr double hydrogenBondLength = 0.7122300374 * ANGSTROM_2_BOHR;

        Atom H1("H", {0.0, 0.0, 0.0});
        Atom H2("H", {hydrogenBondLength, 0.0, 0.0});

        Molecule hydrogenMol({H1, H2});
        hydrogenMol.ConstructBonding(hydrogenBondLength * hydrogenBondLength + 0.1f);

        return hydrogenMol;
    }

    Molecule Molecule::GenerateLithiumHydrideMolecule() {
        constexpr double LiHBondLength = 1.60 * ANGSTROM_2_BOHR;

        Atom H("H", {0.0, 0.0, 0.0});
        Atom Li("Li", {LiHBondLength, 0.0, 0.0});

        Molecule LiHMolecule({H, Li});
        LiHMolecule.ConstructBonding(LiHBondLength * LiHBondLength);

        return LiHMolecule;
    }

    Molecule Molecule::GenerateCOMolecule() {
        constexpr double COBondLength = 1.128 * ANGSTROM_2_BOHR;

        Atom C("C", {0.0, 0.0, 0.0});
        Atom O("O", {COBondLength, 0.0, 0.0});

        Molecule COMolecule({C, O});
        COMolecule.ConstructBonding(COBondLength * COBondLength + 0.1f);

        return COMolecule;
    }

    Molecule Molecule::GenerateMethaneMolecule() {
        Atom H1("H", {-1.183771681898, -1.183771681898, 1.183771681898});
        Atom H2("H", {-1.183771681898, 1.183771681898, -1.183771681898});
        Atom H3("H", {1.183771681898, 1.183771681898, 1.183771681898});
        Atom H4("H", {1.183771681898, -1.183771681898, -1.183771681898});
        Atom C("C", {0, 0, 0});

        Molecule methane({H1, H2, H3, H4, C});
        methane.ConstructBonding(5.0f);

        return methane;
    }

    Molecule Molecule::GenerateEtheneMolecule() {
        Atom C1(
            "C", Eigen::Vector3d{-0.66754171391078, -0.57999587657014, -0.98766494519065} * ANGSTROM_2_BOHR);
        Atom C2(
            "C", Eigen::Vector3d{-0.00280161842801, -0.00243411244334, -0.00414503592554} * ANGSTROM_2_BOHR);
        Atom H1(
            "H", Eigen::Vector3d{-0.15803343371248, -1.06435755385651, -1.81247599918377} * ANGSTROM_2_BOHR);
        Atom H2(
            "H", Eigen::Vector3d{-1.75092796040950, -0.59425022738565, -1.01193846765821} * ANGSTROM_2_BOHR);
        Atom H3(
            "H", Eigen::Vector3d{1.08058463979347, 0.01182021009784, 0.02012843810545} * ANGSTROM_2_BOHR);
        Atom H4(
            "H", Eigen::Vector3d{-0.51230991333269, 0.48192756015780, 0.82066600985273} * ANGSTROM_2_BOHR);

        Molecule ethene({C1, C2, H1, H2, H3, H4});
        ethene.ConstructBonding(7.0f);

        return ethene;
    }

    Molecule Molecule::GenerateAmmoniaMolecule() {
        Eigen::Vector3d NPos{-6.6842683864, 1.9896777361, -0.0017590331};
        NPos *= ANGSTROM_2_BOHR;
        Atom N("N", {0, 0, 0});
        Atom H1("H", Eigen::Vector3d{-5.6844990792, 1.9942938061, 0.0072441956} * ANGSTROM_2_BOHR - NPos);
        Atom H2("H", Eigen::Vector3d{-7.0127564806, 2.9323453345, 0.0540465122} * ANGSTROM_2_BOHR - NPos);
        Atom H3("H", Eigen::Vector3d{-7.0127560539, 1.4847331233, 0.7962183253} * ANGSTROM_2_BOHR - NPos);

        Molecule NH3({N, H1, H2, H3});
        NH3.ConstructBonding(4.0f);

        return NH3;
    }
}
