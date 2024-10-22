#pragma once

#include "Chemistry/Quantum/Basis/Gaussian.h"

#include<glm/vec3.hpp>
#include<vector>
#include<string_view>
#include<utility>
#include<unordered_map>
#include<unordered_set>
#include<Eigen/Core>

namespace SCESS {
    struct BondPair {
        size_t atom1 = 0;
        size_t atom2 = 0;
    };

    inline bool operator==(const BondPair &a, const BondPair &b) {
        return (a.atom1 == b.atom1 && a.atom2 == b.atom2) ||
               (a.atom1 == b.atom2 && a.atom2 == b.atom1);
    }
}

template<>
struct std::hash<SCESS::BondPair> {
    std::size_t operator()(const SCESS::BondPair &k) const {
        std::hash<size_t> hashObj;
        return hashObj(k.atom1) + hashObj(k.atom2);
    }
};

namespace SCESS {
    constexpr double ANGSTROM_2_BOHR = 1.8897259886;

    class Atom {
    private:
        Eigen::Vector3d mPosition;
        std::string mSymbol;
        int mCharge;

    public:
        Atom(const Atom &a) : mSymbol(a.mSymbol),
                              mCharge(a.mCharge), mPosition(a.mPosition) {
        }

        Atom(Atom &&a) noexcept : mSymbol(a.mSymbol), mCharge(a.mCharge), mPosition(std::move(a.mPosition)) {
        }

        Atom(const char *atomSymbol,
             const Eigen::Vector3d &pos = {0., 0., 0.})
            : mSymbol(atomSymbol), mPosition(pos), mCharge(0) {
        }


        Atom(const char *atomSymbol,
             int charge, const Eigen::Vector3d &pos = {0., 0., 0.})
            : mSymbol(atomSymbol), mCharge(charge), mPosition(pos) {
        }

        void SetPos(const Eigen::Vector3d &pos) {
            mPosition = pos;
        }

        const Eigen::Vector3d GetPos() const {
            return mPosition;
        }

        void SetAtomSymbol(const std::string &s) {
            mSymbol = s;
        }

        const std::string GetAtomSymbol() const {
            return mSymbol;
        }

        void SetCharge(int charge) {
            mCharge = charge;
        }

        const int GetCharge() const {
            return mCharge;
        }
    };

    class Molecule {
    private:
        std::vector<Atom> mAtoms;
        std::unordered_set<BondPair> mBonds;
        std::vector<ContractedGaussian> mOrbitals;
        size_t mNoElectrons;

    public:
        Molecule() : mNoElectrons(0) {
        }
        Molecule(const Molecule &m) : mAtoms(m.mAtoms), mBonds(m.mBonds), mOrbitals(m.mOrbitals),
                                      mNoElectrons(m.mNoElectrons) {
        }
        Molecule(Molecule &&m) noexcept : mAtoms(std::move(m.mAtoms)), mBonds(std::move(m.mBonds)),
                                          mOrbitals(std::move(m.mOrbitals)), mNoElectrons(m.mNoElectrons) {
        }
        explicit Molecule(const std::vector<Atom> &atomsIn) : mAtoms(atomsIn), mNoElectrons(0) {
        }
        explicit Molecule(std::vector<Atom> &&atomsIn) : mAtoms(std::move(atomsIn)), mNoElectrons(0) {
        }

        static Molecule LoadFromXYZ(const std::string &filePath);

        void InsertAtom(const Atom &a) {
            mAtoms.emplace_back(a);
        }

        void MakeOrbitals(const std::string &basisSetName);

        const Atom &GetAtom(size_t index) const {
            return mAtoms[index];
        }

        Atom &GetAtom(size_t index) {
            return mAtoms[index];
        }

        const std::vector<Atom> &GetAtoms() const {
            return mAtoms;
        }

        const size_t GetNoAtoms() const {
            return mAtoms.size();
        }

        std::vector<Atom> &GetAtoms() {
            return mAtoms;
        }

        /// <summary>
        ///
        /// </summary>
        /// <param name="tolerance2">The squared distance between atoms where a bond is considered in bohr radii.</param>
        void ConstructBonding(double tolerance2 = 0.5);

        std::unordered_set<BondPair> &GetBonds() {
            return mBonds;
        }

        const std::unordered_set<BondPair> &GetBonds() const {
            return mBonds;
        }

        ContractedGaussian &GetOrbital(uint32_t i) {
            return mOrbitals[i];
        }

        const ContractedGaussian &GetOrbital(uint32_t i) const {
            return mOrbitals[i];
        }

        std::vector<ContractedGaussian> &GetOrbitals() {
            return mOrbitals;
        }

        const std::vector<ContractedGaussian> &GetOrbitals() const {
            return mOrbitals;
        }

        size_t GetNoOrbitals() const {
            return mOrbitals.size();
        }

        size_t GetNoElectrons() const {
            return mNoElectrons;
        }

        double CalculateNuclearRepulsion() const;
    };

    struct AtomStyle {
        //RGB colour
        glm::vec3 mColour;
        //In Bohr radii
        float mRadius;

        static std::unordered_map<std::string_view, AtomStyle> sAtomStyle;


        inline static AtomStyle GetStyle(const std::string_view atom) {
            auto it = sAtomStyle.find(atom);
            if (it != sAtomStyle.end())
                return it->second;
            return {{0.1, 0.1, 0.1}, 1.0f};
        }
    };

    Molecule GenerateWaterMolecule();

    Molecule GenerateHydrogenMolecule();

    Molecule GenerateLithiumHydrideMolecule();

    Molecule GenerateCOMolecule();

    Molecule GenerateMethaneMolecule();

    Molecule GenerateEtheneMolecule();

    Molecule GenerateAmmoniaMolecule();

    Molecule GenerateCarbeneMolecule();
}
