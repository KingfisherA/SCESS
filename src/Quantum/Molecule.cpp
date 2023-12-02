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
}