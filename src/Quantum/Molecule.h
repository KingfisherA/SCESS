#pragma once

#include "Basis/Gaussian.h"

#include<vector>
#include<string_view>
#include<Eigen/Core>
#include<utility>
#include<unordered_map>

namespace SCESS
{
	constexpr double ANGSTROM_2_BOHR = 1.8897259886;

	class Atom {
	private:
		Eigen::Vector3d mPosition;
		std::string mSymbol;
		int mCharge;
	public:

		Atom(const char* atomID) noexcept :mSymbol(mSymbol) {}

		Atom(const Atom& a) : mSymbol(a.mSymbol),
			mCharge(a.mCharge), mPosition(a.mPosition)
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
		std::vector<std::pair<unsigned int, unsigned int>> mBonds;
		std::vector<GaussianCollection> mOrbitals;
		size_t mNoElectrons;
	public:

		Molecule() :mNoElectrons(0) { }
		Molecule(const Molecule& m) :mAtoms(m.mAtoms), mBonds(m.mBonds), mOrbitals(m.mOrbitals), mNoElectrons(m.mNoElectrons) {}
		Molecule(Molecule&& m) noexcept :mAtoms(std::move(m.mAtoms)), mBonds(std::move(m.mBonds)), mOrbitals(std::move(m.mOrbitals)), mNoElectrons(m.mNoElectrons) {}
		Molecule(const std::vector<Atom>& atomsIn) : mAtoms(atomsIn), mNoElectrons(0)
		{
			ConstructBonding();
		}
		Molecule(std::vector<Atom>&& atomsIn) :mAtoms(std::move(atomsIn)), mNoElectrons(0)
		{
			ConstructBonding();
		}

		void InsertAtom(const Atom& a)
		{
			mAtoms.emplace_back(a);
		}

		void MakeOrbitals(const std::string& basisSetName);

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

		const size_t GetNoAtoms() const
		{
			return mAtoms.size();
		}

		std::vector<Atom>& GetAtoms()
		{
			return mAtoms;
		}

		void ConstructBonding();

		std::vector<std::pair<unsigned int, unsigned int>>& GetBonds()
		{
			return mBonds;
		}

		const std::vector<std::pair<unsigned int, unsigned int>>& GetBonds() const
		{
			return mBonds;
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

		size_t GetNoOrbitals() const {
			return mOrbitals.size();
		}

		size_t GetNoElectrons() const
		{
			return mNoElectrons;
		}

		double CalculateNuclearRepulsion() const;
	};

	struct AtomStyle
	{
		//RGB colour
		glm::vec3 mColour;
		//In Bohr radii
		float mRadius;

		static std::unordered_map<std::string_view, AtomStyle> sAtomStyle;



		inline static const AtomStyle& GetStyle(const std::string_view atom)
		{
			auto it = sAtomStyle.find(atom);
			if (it != sAtomStyle.end())
				return it->second;
			return { {0.1, 0.1, 0.1}, 1.0f };
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