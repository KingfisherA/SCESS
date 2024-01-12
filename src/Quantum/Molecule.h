#pragma once

#include "Basis/Gaussian.h"

#include<vector>
#include<string_view>
#include<Eigen/Core>
#include<utility>

namespace SCESS 
{
    constexpr double ANGSTROM_2_BOHR = 1.8897259886;

    class Atom {
    private:
        Eigen::Vector3d mPosition;
        std::string mSymbol;
        int mCharge;
    public:

        Atom(const char* atomID) noexcept :mSymbol(mSymbol){}
        
        Atom(const Atom& a) : mSymbol(a.mSymbol), 
            mCharge(a.mCharge),  mPosition(a.mPosition)
        {

        }

        Atom(Atom&& a) noexcept : mSymbol(a.mSymbol), mCharge(a.mCharge), mPosition(std::move(a.mPosition))
        {}

        Atom(const char* atomSymbol,
             const Eigen::Vector3d& pos = { 0., 0., 0. })
            : mSymbol(atomSymbol), mPosition(pos)
        {

        }


        Atom(const char* atomSymbol,
            int charge, const Eigen::Vector3d& pos = { 0., 0., 0. })
            : mSymbol(atomSymbol), mCharge(charge), mPosition(pos)
        {
    
        }


        const Eigen::Vector3d GetPos()const
        {
            return mPosition;
        }

        const std::string GetAtomSymbol() const
        {
            return mSymbol;
        }

        void SetCharge(int charge)
        {
            mCharge = charge;
        }

        const int GetCharge() const
        {
            return mCharge;
        }
    };

    class Molecule
    {
    private:
        std::vector<Atom> mAtoms;
        std::vector<GaussianCollection> mOrbitals;
        size_t mNoElectrons;
    public:

        Molecule():mNoElectrons(0){}
        Molecule(const Molecule& m):mAtoms(m.mAtoms), mOrbitals(m.mOrbitals), mNoElectrons(m.mNoElectrons){}
        Molecule(Molecule&& m) noexcept:mAtoms(std::move(m.mAtoms)), mOrbitals(std::move(m.mOrbitals)), mNoElectrons(m.mNoElectrons) {}
        Molecule(const std::vector<Atom>& atomsIn) :mAtoms(atomsIn), mNoElectrons(0)
        {}
        Molecule(std::vector<Atom>&& atomsIn) :mAtoms(std::move(atomsIn)), mNoElectrons(0)
        {}

        void AddAtom(const Atom& a)
        {
            mAtoms.emplace_back(a);
        }

        void MakeOrbitals();

        const Atom& GetAtom(size_t index) const
        {
            return mAtoms[index];
        }

        Atom& GetAtom(size_t index)
        {
            return mAtoms[index];
        }

        const std::vector<Atom>& GetAtoms() const
        {
            return mAtoms;
        }

        std::vector<Atom>& GetAtoms()
        {
            return mAtoms;
        }

        GaussianCollection& GetOrbital(uint32_t i)
        {
            return mOrbitals[i];
        }

        const GaussianCollection& GetOrbital(uint32_t i) const
        {
            return mOrbitals[i];
        }

        std::vector<GaussianCollection>& GetOrbitals()
        {
            return mOrbitals;
        }

        const std::vector<GaussianCollection>& GetOrbitals() const
        {
            return mOrbitals;
        }

        size_t NoOrbitals() const {
            return mOrbitals.size();
        }

        size_t GetNoElectrons() const
        {
            return mNoElectrons;
        }
    };

    Molecule GenerateWaterMolecule();

    Molecule GenerateHydrogenMolecule();

    Molecule GenerateLithiumHydrideMolecule();

    Molecule GenerateCOMolecule();

    Molecule GenerateMethaneMolecule();

    Molecule GenerateEtheneMolecule();

    Molecule GenerateAmmoniaMolecule();
}