#pragma once

#include<string_view>
#include<vector>
#include<nlohmann/json.hpp>

#include "Quantum/Basis/Gaussian.h"
#include "Quantum/Integrals/IntegralsLibInt.h"

namespace SCESS {
    class Molecule;
}

namespace SCESS {


    class BasisSet {
    public:
        BasisSet() = default;
        BasisSet(const std::string &basisSetName, const Molecule &mol);
        BasisSet(const BasisSet&) = default;
        BasisSet(BasisSet&&) = default;

        BasisSet& operator=(const BasisSet&) = default;
        BasisSet& operator=(BasisSet&&) = default;

        static std::unordered_set<std::string> ListBasisSets();

        void LoadAndAddBasisSet(const std::string &basisSetName, const Molecule &mol);
        void AddCGF(const ContractedGaussian &cgf);
        template<typename Iterator>
        void AddCGFs(const Iterator begin, const Iterator end) {
            mBasisFuncs.insert(mBasisFuncs.end(), begin, end);
        }
        void Clear() { mBasisFuncs.clear(); }
        size_t GetNoCGFs() const { return mBasisFuncs.size(); }
        size_t GetNoPrimitives() const;

        ContractedGaussian& Get(size_t i){return mBasisFuncs[i];}
        const ContractedGaussian& Get(size_t i) const{return mBasisFuncs[i];}
        std::vector<ContractedGaussian> &GetCGFs() { return mBasisFuncs; }
        const std::vector<ContractedGaussian> &GetCGFs() const { return mBasisFuncs; }

    private:
        std::vector<ContractedGaussian> mBasisFuncs;
    };


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
            return {GetNoHFFuncs(), mNumOfFuncs};
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
}
