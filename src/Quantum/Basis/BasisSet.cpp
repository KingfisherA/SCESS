#include<Quantum/QuantumPCH.h>

#include "Resources/AssetFolder.h"
#define DATADIR SCESS::Res::AssetsFolder::Get().GetResourcePath("Basis").string()
#include <libint2/basis.h>

#include "Quantum/Molecule.h"
#include "Quantum/Basis/BasisSet.h"
#include "Quantum/HF/HartreeFock.h"
#include <boost/filesystem.hpp>

#include "Utilities/Utilities.h"


//
//
namespace SCESS {
    BasisSet::BasisSet(const std::string &basisSetName, const Molecule &mol) {
        LoadAndAddBasisSet(basisSetName, mol);
    }

    std::unordered_set<std::string> BasisSet::ListBasisSets() {
        namespace fs = boost::filesystem;
        std::unordered_set<std::string> basisSets;
        fs::path basisSetPath = Res::AssetsFolder::Get().GetResourcePath("Basis");

        fs::recursive_directory_iterator itr(basisSetPath);

        fs::recursive_directory_iterator end;

        while (itr != end) {
            if (fs::is_regular_file(itr->path()) && itr->path().extension().string() == ".json") {
                basisSets.insert(itr->path().stem().string());
            }

            boost::system::error_code ec;
            itr.increment(ec);
        }

        return basisSets;
    }

    void BasisSet::LoadAndAddBasisSet(const std::string &basisSetName, const Molecule &mol) {
        std::string jsonFilePath = Res::AssetsFolder::Get().GetResourcePath("Basis")
                .append(basisSetName + ".json").string();

        if (!std::filesystem::exists(jsonFilePath)) {
            throw std::runtime_error("Cannot locate JSON basis set");
        }

        std::ifstream jsonFile(jsonFilePath.data());
        nlohmann::json basisSet;
        try {
            basisSet = nlohmann::json::parse(jsonFile);
        } catch (const std::exception &e) {
            jsonFile.close();
            throw;
        }
        jsonFile.close();

        std::vector<Eigen::Vector3i> angularMomenta;
        auto GenerateAngularMomenta = [&angularMomenta](uint32_t angMom) {
            angularMomenta.clear();
            angularMomenta.reserve(2 * angMom + 1);

            for (uint32_t x = 0; x <= angMom; x++)
                for (uint32_t y = 0; y <= angMom; y++)
                    for (uint32_t z = 0; z <= angMom; z++) {
                        if (x + y + z == angMom)
                            angularMomenta.emplace_back(x, y, z);
                    }
        };


        for (const Atom &a: mol.GetAtoms()) {
            nlohmann::json shells = basisSet["elements"][std::to_string(a.GetCharge())]["electron_shells"];
            for (size_t i = 0; i < shells.size(); i++) {
                nlohmann::json &shell = shells[i];
                nlohmann::json angMomenta = shell["angular_momentum"];
                for (size_t j = 0; j < angMomenta.size(); j++) {
                    int angularMomentum = angMomenta[j].get<int>();
                    GenerateAngularMomenta((uint32_t) angularMomentum);

                    nlohmann::json exponents = shell["exponents"];
                    nlohmann::json coefficients = shell["coefficients"];

                    for (const Eigen::Vector3i &angMom: angularMomenta) {
                        ContractedGaussian cgf(a.GetPos());
                        for (nlohmann::json &coeffCollection: coefficients) {
                            for (size_t k = 0; k < coeffCollection.size(); k++) {
                                double coefficient = std::stod(coeffCollection[k].get<std::string>());
                                if (coefficient == 0.0)
                                    continue;
                                double exponent = std::stod(exponents[k].get<std::string>());
                                cgf.Add(angMom, Gaussian(exponent, coefficient));
                            }
                        }
                        mBasisFuncs.push_back(cgf);
                    }
                }
            }
        }
    }

    void BasisSet::AddCGF(const ContractedGaussian &cgf) {
        mBasisFuncs.push_back(cgf);
    }

    size_t BasisSet::GetNoPrimitives() const {
        size_t numPrimitives = 0;

        for (const ContractedGaussian &cgf: mBasisFuncs)
            numPrimitives += cgf.Size();

        return numPrimitives;
    }


    void F12BasisSet::AddBasisSets(const std::string &hfBasisSetName, const std::string &cabsBasisSetName,
                                   const std::vector<libint2::Atom> &molecule) {
        mHFShellsStartIdx = mShells.size();
        libint2::BasisSet hfBasisSet(hfBasisSetName, molecule, true);
        mShells = hfBasisSet.shells();

        mCABSShellsSetStartIdx = mShells.size();

        libint2::BasisSet cabsBasisSet(cabsBasisSetName, molecule, true);
        auto &basisSetShells = cabsBasisSet.shells();
        mShells.insert(mShells.end(), basisSetShells.begin(), basisSetShells.end());

        mNumOfFuncs = 0;
        mMaxL = 0;
        mNumPrimitives = 0;
        for (auto &&shell: mShells) {
            mNumOfFuncs += shell.size();
            mNumPrimitives = std::max(shell.nprim(), mNumPrimitives);
            for (auto &&c: shell.contr)
                mMaxL = std::max(c.l, mMaxL);
        }
    }

    std::vector<size_t> F12BasisSet::ShellToBasisIdx() const {
        std::vector<size_t> result(mShells.size(), 0);

        size_t counter = 0;
        for (size_t i = 0; i < mShells.size(); i++) {
            result[i] = counter;
            counter += mShells[i].size();
        }

        return result;
    }
}
