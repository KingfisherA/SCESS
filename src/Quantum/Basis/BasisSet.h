#pragma once
#include<filesystem>
#include<string>
#include<vector>
#include<unordered_set>
#include<libint2.hpp>
#include<nlohmann/json.hpp>
//
#include "CGF.h"

namespace SCESS {
    class Molecule;
}

//
//
namespace SCESS {
    //
    //
    //A basis set file in Assets/Basis. Basis sets are organised into one folder per family; the family is the name of
    //the folder and the basis set name is the file name without the .json extension.
    struct BasisSetInfo {
        std::string name;
        std::string family;
        std::filesystem::path path;
    };
    //
    //The basis sets available in Assets/Basis (Basis Set Exchange JSON files)
    class BasisLibrary {
    public:
        //Family of the complementary auxiliary basis sets, which are not used as orbital basis sets
        static constexpr const char *CABS_FAMILY = "CABS";

        //Every basis set, ordered by family (Pople, Dunning, Karlsruhe, F12, CABS, then any other) and then name
        static const std::vector<BasisSetInfo> &List();
        //The families in display order
        static std::vector<std::string> ListFamilies();
        //Case-insensitive lookup which also accepts the 6-31G* style names of the 6-31G(d) style files.
        //Returns nullptr if there is no such basis set.
        static const BasisSetInfo *Find(const std::string &name);
        //As Find() but throws if the basis set does not exist
        static const BasisSetInfo &Get(const std::string &name);
        //The recommended CABS for an orbital basis set (e.g. cc-pVDZ-F12 -> cc-pVDZ-F12-OptRI) or an empty string
        static std::string DefaultCABS(const std::string &obsName);

        static bool IsCABS(const BasisSetInfo &info) { return info.family == CABS_FAMILY; }

        //Loads the parsed JSON file of a basis set
        static nlohmann::json LoadJSON(const std::string &name);
        //Builds the libint2 shells of a basis set for the given atoms. Shells use solid harmonics unless the basis set
        //file marks them as Cartesian (gto_cartesian). Throws if an atom has no basis functions or a shell exceeds the
        //maximum angular momentum libint2 was built with.
        static std::vector<libint2::Shell> LoadShells(const std::string &name, const std::vector<libint2::Atom> &atoms);
        static libint2::BasisSet LoadLibIntBasisSet(const std::string &name, const std::vector<libint2::Atom> &atoms) {
            return libint2::BasisSet(LoadShells(name, atoms));
        }
    };
    //
    //
    class BasisSet {
    public:
        BasisSet() = default;
        BasisSet(const std::string &basisSetName, const Molecule &mol);
        BasisSet(const BasisSet &) = default;
        BasisSet(BasisSet &&) = default;

        BasisSet &operator=(const BasisSet &) = default;
        BasisSet &operator=(BasisSet &&) = default;

        //Names of every basis set in Assets/Basis, see BasisLibrary for families and ordering
        static std::vector<std::string> ListBasisSets();

        void LoadAndAddBasisSet(const std::string &basisSetName, const Molecule &mol);
        void AddCGF(const CGF &cgf);

        template<typename Iterator>
        void AddCGFs(const Iterator begin, const Iterator end) {
            //Hand-built functions have no shells so the basis set becomes CGF-only
            mShells.clear();
            mCartesianToSpherical.resize(0, 0);
            mBasisFuncs.insert(mBasisFuncs.end(), begin, end);
        }

        void Clear() {
            mBasisFuncs.clear();
            mShells.clear();
            mCartesianToSpherical.resize(0, 0);
        }
        //Number of Cartesian contracted Gaussians (the functions evaluated on grids and rendered)
        size_t GetNoCGFs() const { return mBasisFuncs.size(); }
        //Number of basis functions of the integrals and orbitals: the libint2 shells, which are solid harmonic
        //unless the basis set file marks them as Cartesian. For CGF-only basis sets (built with AddCGF) this is the
        //number of CGFs.
        size_t GetNoFunctions() const {
            return HasShells() ? (size_t) mCartesianToSpherical.cols() : mBasisFuncs.size();
        }
        //True if the basis set was loaded by name so that the libint2 shells (and the CGFs) describe it
        bool HasShells() const { return !mShells.empty(); }
        const std::vector<libint2::Shell> &GetShells() const { return mShells; }
        //T (no. CGFs x no. functions) expresses every basis function in the Cartesian CGFs,
        //phi_s(r) = sum_c T(c, s) phi_c(r). Use it to evaluate orbitals (C_CGF = T C) or basis functions on a grid.
        const Eigen::MatrixXd &CartesianToSpherical() const { return mCartesianToSpherical; }
        size_t GetNoPrimitives() const;

        CGF &Get(size_t i) { return mBasisFuncs[i]; }
        const CGF &Get(size_t i) const { return mBasisFuncs[i]; }
        std::vector<CGF> &GetCGFs() { return mBasisFuncs; }
        const std::vector<CGF> &GetCGFs() const { return mBasisFuncs; }

    private:
        void BuildCartesianToSpherical();

        std::vector<CGF> mBasisFuncs;
        std::vector<libint2::Shell> mShells;
        Eigen::MatrixXd mCartesianToSpherical;
    };
    //
    //
    class F12BasisSet {
    public:
        using Shells_T = std::vector<libint2::Shell>;
        F12BasisSet() = default;
        F12BasisSet(const F12BasisSet &) = default;

        F12BasisSet(const std::string &hfBasisSetName,
                    const std::string &cabsBasisSetName, const std::vector<libint2::Atom> &molecule) {
            AddBasisSets(hfBasisSetName, cabsBasisSetName, molecule);
        }

        void AddBasisSets(const std::string &hfBasisSetName,
                          const std::string &cabsBasisSetName, const std::vector<libint2::Atom> &molecule);

        void Clear() { mShells.clear(); }

        const libint2::Shell &GetShell(size_t idx) const {
            return mShells[idx];
        }

        libint2::Shell &GetShell(size_t idx) {
            return mShells[idx];
        }

        libint2::Shell &operator[](size_t idx) { return mShells[idx]; }
        const libint2::Shell &operator[](size_t idx) const { return mShells[idx]; }

        size_t GetNoFuncsInShell(size_t idx) const { return mShells[idx].size(); }


        [[nodiscard]] const Shells_T &GetShells() const { return mShells; }
        Shells_T &GetShells() { return mShells; }


        [[nodiscard]] std::pair<size_t, size_t> GetHFShellIndices() const {
            return {mHFShellsStartIdx, mCABSShellsSetStartIdx - 1};
        }

        [[nodiscard]] std::pair<size_t, size_t> GetCABSShellIndices() const {
            return {mCABSShellsSetStartIdx, mShells.size() - 1};
        }

        [[nodiscard]] std::pair<size_t, size_t> GetHFIndices() const {
            return {mHFShellsStartIdx, GetNoHFFuncs() - 1};
        }

        [[nodiscard]] std::pair<size_t, size_t> GetCABSIndices() const {
            return {GetNoHFFuncs(), mNumOfFuncs - 1};
        }

        std::vector<size_t> ShellToBasisIdx() const;

        [[nodiscard]] size_t GetNoHFShells() const {
            return mCABSShellsSetStartIdx - mHFShellsStartIdx;
        }

        [[nodiscard]] size_t GetNoHFFuncs() const {
            size_t result = 0;
            for (size_t i = 0; i < mCABSShellsSetStartIdx; i++)
                result += mShells[i].size();
            return result;
        }

        [[nodiscard]] size_t GetNoCABSShells() const { return mShells.size() - mCABSShellsSetStartIdx; }

        [[nodiscard]] size_t GetNoCABSFuncs() const {
            size_t result = 0;
            for (size_t i = mCABSShellsSetStartIdx; i < mShells.size(); i++)
                result += mShells[i].size();
            return result;
        }

        [[nodiscard]] size_t GetNumShells() const { return mShells.size(); }
        [[nodiscard]] size_t GetNumOfFuncs() const { return mNumOfFuncs; }
        [[nodiscard]] size_t GetNumOfPrimitives() const { return mNumPrimitives; }
        [[nodiscard]] int GetMaxL() const { return mMaxL; }

    private:
        Shells_T mShells;
        size_t mHFShellsStartIdx = 0,
                mCABSShellsSetStartIdx = 0,
                mNumOfFuncs = 0,
                mNumPrimitives = 0;
        int mMaxL = 0;
    };
    //
    //
}
