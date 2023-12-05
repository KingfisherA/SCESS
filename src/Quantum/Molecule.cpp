#include<pch.h>

#include"Molecule.h"
#include"Basis/BasisUtils.h"
#include"Basis/Gaussian.h"

//
//
namespace SCESS
{
	void Molecule::MakeOrbitals()
	{
		std::optional<nlohmann::json> basisSetOpt = BasisUtils::LoadBasisSetJSON(R"(C:\dev\Chemistry\SCESS\Data\sto-3g.json)");
		if (!basisSetOpt.has_value())
		{
			throw std::exception("Basis set could not be loaded succesfully");
		}
		nlohmann::json& basisSet = *basisSetOpt;

		mOrbitals.clear();
		for (Atom& a : mAtoms)
		{
			nlohmann::json& atomBasis = basisSet[a.GetAtomSymbol()];
			a.SetCharge(atomBasis["atomic_number"].get<int>());
			nlohmann::json& cgfs = atomBasis["cgfs"];
			for (size_t i = 0; i < cgfs.size(); i++)
			{
				nlohmann::json& AO = cgfs[i];
				std::string AOType;
				AO["type"].get_to(AOType);
				
				std::vector<SCESS::OrbitalType> types;
				if (AOType == "S")
				{
					types.push_back(SCESS::OrbitalType::S);
				}
				else if (AOType == "P")
				{
					types.push_back(SCESS::OrbitalType::P_X);
					types.push_back(SCESS::OrbitalType::P_Y);
					types.push_back(SCESS::OrbitalType::P_Z);
				}
				else if (AOType == "D")
				{
					types.push_back(SCESS::OrbitalType::D_X2);
					types.push_back(SCESS::OrbitalType::D_Y2);
					types.push_back(SCESS::OrbitalType::D_Z2);
					types.push_back(SCESS::OrbitalType::D_XY);
					types.push_back(SCESS::OrbitalType::D_XZ);
					types.push_back(SCESS::OrbitalType::D_YZ);
				}
				else if (AOType == "F")
				{
					types.push_back(SCESS::OrbitalType::F_X3);
					types.push_back(SCESS::OrbitalType::F_Y3);
					types.push_back(SCESS::OrbitalType::F_Z3);
					types.push_back(SCESS::OrbitalType::F_X2_Y);
					types.push_back(SCESS::OrbitalType::F_X_Y2);
					types.push_back(SCESS::OrbitalType::F_X2_Z);
					types.push_back(SCESS::OrbitalType::F_X_Z2);
					types.push_back(SCESS::OrbitalType::F_Y2_Z);
					types.push_back(SCESS::OrbitalType::F_Y_Z2);
					types.push_back(SCESS::OrbitalType::F_XYZ);
				}

				nlohmann::json& gtos = AO["gtos"];
				for (SCESS::OrbitalType type : types)
				{
					GaussianCollection collection(a.GetPos());
					for (nlohmann::json& gto : gtos)
					{
						double exponent, coeff=0.0;
						gto["alpha"].get_to(exponent);
						gto["coeff"].get_to(coeff);
						collection.Add(type, Gaussian(exponent, coeff));
					}
					mOrbitals.push_back(collection);
				}
			}
		}
	}

	Molecule GenerateWaterMolecule()
	{
		constexpr double waterBondLength = 0.957*ANGSTROM_2_BOHR;
		constexpr double waterBondAngle = 104.5*std::numbers::pi/180.0;
		SCESS::Atom H1("H", Eigen::Vector3d{ waterBondLength * std::sin(waterBondAngle * 0.5),waterBondLength*std::cos(waterBondAngle*0.5),0.0});
		SCESS::Atom H2("H", Eigen::Vector3d{ -waterBondLength * std::sin(waterBondAngle * 0.5),waterBondLength * std::cos(waterBondAngle * 0.5),0.0 });
		SCESS::Atom O("O", Eigen::Vector3d{ 0.0,0.0, 0.0});

		SCESS::Molecule water({ H1, O, H2 });
		water.MakeOrbitals();

		return water;
	}

	Molecule GenerateHydrogenMolecule()
	{
		constexpr double hydrogenBondLength = 0.7122300374 * ANGSTROM_2_BOHR;

		SCESS::Atom H1("H", { 0.0,0.0,0.0 });
		SCESS::Atom H2("H", { hydrogenBondLength,0.0,0.0 });

		SCESS::Molecule hydrogenMol({ H1, H2 });
		hydrogenMol.MakeOrbitals();

		return hydrogenMol;
	}

	Molecule GenerateLithiumHydrideMolecule()
	{
		constexpr double LiHBondLength = 1.60 * ANGSTROM_2_BOHR;

		SCESS::Atom H("H", { 0.0,0.0,0.0 });
		SCESS::Atom Li("Li", { LiHBondLength,0.0,0.0 });

		SCESS::Molecule LiHMolecule({ H, Li });
		LiHMolecule.MakeOrbitals();

		return LiHMolecule;
	}


	Molecule GenerateCOMolecule()
	{
		constexpr double COBondLength = 1.128 * ANGSTROM_2_BOHR;

		SCESS::Atom C("C", { 0.0,0.0,0.0 });
		SCESS::Atom O("O", { COBondLength,0.0,0.0 });

		SCESS::Molecule COMolecule({ C,O });
		COMolecule.MakeOrbitals();

		return COMolecule;
	}
}