#include<pch.h>

#include "Resources/Resources.h"
#define DATADIR SCESS::Res::AssetsFolder::Get().GetResourcePath("Basis").string()
#include <libint2/basis.h>
#include<ranges>

#include "Quantum/Basis/BasisUtils.h"
#include "Quantum/HF/HartreeFock.h"
#include <boost/filesystem.hpp>

#include "Utilities/Utilities.h"


//
//
namespace SCESS::BasisUtils {
    std::unordered_set<std::string> ListBasisSets() {
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

    std::optional<nlohmann::json> LoadBasisSetJSON(
        std::string_view jsonFilePath) {
        if (!std::filesystem::exists(jsonFilePath)) {
            return std::nullopt;
        }

        std::ifstream jsonFile(jsonFilePath.data());

        nlohmann::json inputtedJson;
        try {
            inputtedJson = nlohmann::json::parse(jsonFile);
        } catch (const std::exception &e) {
            std::cout << "Could not parse JSON. \n" << e.what() << '\n';
            return std::nullopt;
        }

        jsonFile.close();

        return inputtedJson;
    }

    std::vector<Eigen::Vector3i> GenerateAngularMomenta(const size_t angMom) {
        std::vector<Eigen::Vector3i> angularMom;
        angularMom.reserve(2 * angMom + 1);

        for (size_t x = 0; x <= angMom; x++)
            for (size_t y = 0; y <= angMom; y++)
                for (size_t z = 0; z <= angMom; z++) {
                    if (x + y + z == angMom)
                        angularMom.emplace_back(x, y, z);
                }

        return angularMom;
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
