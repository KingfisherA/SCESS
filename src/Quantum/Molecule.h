#pragma once

#include "Quantum/Basis/Gaussian.h"

#include<vector>
#include<string_view>
#include<utility>
#include<unordered_set>
#include<Eigen/Core>
#include<libint2.hpp>

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
    public:
        Atom() = default;
        Atom(const Atom &a) : mSymbol(a.mSymbol),
                              mCharge(a.mCharge), mPosition(a.mPosition) {
        }
        Atom(Atom &&a) noexcept : mSymbol(a.mSymbol), mCharge(a.mCharge), mPosition(std::move(a.mPosition)) {
        }
        Atom(std::string atomSymbol,
             const Eigen::Vector3d &pos);
        Atom(std::string atomSymbol, int charge, const Eigen::Vector3d &pos): mSymbol(atomSymbol), mCharge(charge),
                                                                              mPosition(pos) {
        }

        Atom &operator=(const Atom &a) {
            mSymbol = a.mSymbol;
            mCharge = a.mCharge;
            mPosition = a.mPosition;
            return *this;
        }

        void SetPos(const Eigen::Vector3d &pos) {
            mPosition = pos;
        }
        const Eigen::Vector3d &GetPos() const {
            return mPosition;
        }
        void SetAtomSymbol(std::string s);
        const std::string &GetAtomSymbol() const {
            return mSymbol;
        }
        void SetCharge(int charge) {
            mCharge = charge;
        }
        int GetCharge() const {
            return mCharge;
        }

    private:
        Eigen::Vector3d mPosition;
        std::string mSymbol;
        int mCharge = 0;
    };

    class Molecule {
    public:
        Molecule() = default;
        Molecule(const Molecule &m) : mAtoms(m.mAtoms), mBonds(m.mBonds),
                                      mNoElectrons(m.mNoElectrons) {
        }
        Molecule(Molecule &&m) noexcept : mAtoms(std::move(m.mAtoms)), mBonds(std::move(m.mBonds)),
                                          mNoElectrons(m.mNoElectrons) {
        }
        explicit Molecule(const std::vector<Atom> &atomsIn): mAtoms(atomsIn) {
            for (const Atom &a: mAtoms)
                mNoElectrons += a.GetCharge();
        }
        explicit Molecule(std::vector<Atom> &&atomsIn) : mAtoms(std::move(atomsIn)) {
            for (const Atom &a: mAtoms)
                mNoElectrons += a.GetCharge();
        }

        void AddAtom(const Atom &a);
        const Atom &GetAtom(size_t index) const {
            return mAtoms[index];
        }
        Atom &GetAtom(size_t index) {
            return mAtoms[index];
        }
        const std::vector<Atom> &GetAtoms() const {
            return mAtoms;
        }
        std::vector<Atom> &GetAtoms() {
            return mAtoms;
        }
        const size_t GetNoAtoms() const {
            return mAtoms.size();
        }
        std::vector<libint2::Atom> GetAtomsLibInt() const;

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
        size_t GetNoElectrons() const {
            return mNoElectrons;
        }

        double CalculateNuclearRepulsion() const;

        static Molecule LoadFromXYZ(const std::string &filePath);
        static Molecule GenerateWaterMolecule();
        static Molecule GenerateHydrogenMolecule();
        static Molecule GenerateLithiumHydrideMolecule();
        static Molecule GenerateCOMolecule();
        static Molecule GenerateMethaneMolecule();
        static Molecule GenerateEtheneMolecule();
        static Molecule GenerateAmmoniaMolecule();
        static Molecule GenerateCarbeneMolecule();

    private:
        std::vector<Atom> mAtoms;
        std::unordered_set<BondPair> mBonds;
        size_t mNoElectrons = 0;
    };
}
