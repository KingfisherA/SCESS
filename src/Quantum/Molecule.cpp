#include<pch.h>
#include<boost/locale/conversion.hpp>

#include"Molecule.h"
#include"Basis/BasisUtils.h"
#include"Basis/Gaussian.h"
#include"Utilities/PathUtils.h"
#include"Chemistry/PeriodicTable.h"

//
//
namespace SCESS
{
	std::unordered_map<std::string_view, AtomStyle> AtomStyle::sAtomStyle =
	{
		{"H", {{1.0f, 1.0f, 1.0f}, 1.0f}},
		{"C", {{0.1f, 0.1f, 0.1f}, 2.0f}},
		{"Li", {{0.8f, 0.9f, 1.0f}, 3.15f}},
		{"N", {{0.0f, 0.5f, 0.5f}, 1.8f}},
		{"O", {{1.0f, 0.0f, 0.0f}, 1.7f}}
	};

	void Molecule::MakeOrbitals(const std::string& basisSetName)
	{
		std::optional<nlohmann::json> basisSetOpt = BasisUtils::LoadBasisSetJSON(PathUtils::GetAssetPath("Basis\\" + basisSetName + ".json").string());
		if (!basisSetOpt.has_value())
		{
			throw std::exception("Basis set could not be loaded succesfully");
		}
		nlohmann::json& basisSet = *basisSetOpt;
		std::vector<SCESS::OrbitalType> types;
		mOrbitals.clear();
		mNoElectrons = 0;

		if (basisSetName == "sto-3g")
		{
			for (Atom& a : mAtoms)
			{
				nlohmann::json atomBasis = basisSet[a.GetAtomSymbol()];
				a.SetCharge(atomBasis["atomic_number"].get<int>());
				mNoElectrons += a.GetCharge();
				nlohmann::json& cgfs = atomBasis["cgfs"];
				for (size_t i = 0; i < cgfs.size(); i++)
				{
					nlohmann::json& AO = cgfs[i];
					std::string AOType;
					AO["type"].get_to(AOType);

					types.clear();
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
							double exponent, coeff = 0.0;
							gto["alpha"].get_to(exponent);
							gto["coeff"].get_to(coeff);
							collection.Add(type, Gaussian(exponent, coeff));
						}
						mOrbitals.push_back(collection);
					}
				}
			}
			return;
		}

		PeriodicTable::OpenQuery();

		for (Atom& a : mAtoms)
		{
			nlohmann::json atomData = PeriodicTable::QueryElementSym(a.GetAtomSymbol());
			a.SetCharge(atomData["Atomic_Number"].get<int>());
			mNoElectrons += a.GetCharge();

			nlohmann::json shells = basisSet["elements"][std::to_string(a.GetCharge())]["electron_shells"];
			for (size_t i = 0; i < shells.size(); i++)
			{
				nlohmann::json& shell = shells[i];
				nlohmann::json angMomenta = shell["angular_momentum"];
				for (size_t j = 0; j < angMomenta.size(); j++)
				{
					int angularMomentum = angMomenta[j].get<int>();

					types.clear();
					switch (angularMomentum)
					{
					case 0:
						types.push_back(SCESS::OrbitalType::S);
						break;
					case 1:
						types.push_back(SCESS::OrbitalType::P_X);
						types.push_back(SCESS::OrbitalType::P_Y);
						types.push_back(SCESS::OrbitalType::P_Z);
						break;
					case 2:
						types.push_back(SCESS::OrbitalType::D_X2);
						types.push_back(SCESS::OrbitalType::D_Y2);
						types.push_back(SCESS::OrbitalType::D_Z2);
						types.push_back(SCESS::OrbitalType::D_XY);
						types.push_back(SCESS::OrbitalType::D_XZ);
						types.push_back(SCESS::OrbitalType::D_YZ);
						break;
					case 3:
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
						break;
					}

					for (SCESS::OrbitalType type : types)
					{
						GaussianCollection collection(a.GetPos());
						nlohmann::json exponents = shell["exponents"];
						nlohmann::json coefficients = shell["coefficients"];
						for (nlohmann::json& coeffCollection : coefficients)
						{
							for (size_t k = 0; k < coeffCollection.size(); k++)
							{
								std::string exponentS, coeffS;
								coeffCollection[k].get_to(coeffS);
								exponents[k].get_to(exponentS);
								collection.Add(type, Gaussian(std::stod(exponentS), std::stod(coeffS)));
							}
						}
						mOrbitals.push_back(collection);
					}
				}
			}
		}

		PeriodicTable::CloseQuery();
	}

	void Molecule::ConstructBonding()
	{
		for (size_t i = 0; i < mAtoms.size()-1; i++)
		{
			std::pair<unsigned int, unsigned int> bond;
			bond.first = i;
			bond.second = 0;

			double lastDist2 = 1E100;
			for (size_t j = i + 1; j < mAtoms.size(); j++)
			{
				double distance2 = (mAtoms[i].GetPos() - mAtoms[j].GetPos()).squaredNorm();
				if (distance2 < lastDist2)
				{
					bond.second = j;
					lastDist2 = distance2;
				}
			}
			mBonds.push_back(bond);
		}
	}

	double Molecule::CalculateNuclearRepulsion() const
	{
		double Vnn = 0.0;
		for (size_t i = 0; i < mAtoms.size(); i++) {
			for (size_t j = i + 1; j < mAtoms.size(); j++)
			{
				double r = (mAtoms[i].GetPos() - mAtoms[j].GetPos()).norm();
				if (r > 0.0) Vnn += mAtoms[i].GetCharge() * mAtoms[j].GetCharge() / r;
			}
		}
		return Vnn;
	}

	Molecule GenerateWaterMolecule()
	{
		constexpr double waterBondLength = 0.957 * ANGSTROM_2_BOHR;
		constexpr double waterBondAngle = 104.5 * std::numbers::pi / 180.0;
		SCESS::Atom H1("H", Eigen::Vector3d{ waterBondLength * std::sin(waterBondAngle * 0.5),waterBondLength * std::cos(waterBondAngle * 0.5),0.0 });
		SCESS::Atom H2("H", Eigen::Vector3d{ -waterBondLength * std::sin(waterBondAngle * 0.5),waterBondLength * std::cos(waterBondAngle * 0.5),0.0 });
		SCESS::Atom O("O", Eigen::Vector3d{ 0.0,0.0, 0.0 });

		SCESS::Molecule water({ H1, O, H2 });

		return water;
	}

	Molecule GenerateHydrogenMolecule()
	{
		constexpr double hydrogenBondLength = 0.7122300374 * ANGSTROM_2_BOHR;

		SCESS::Atom H1("H", { 0.0,0.0,0.0 });
		SCESS::Atom H2("H", { hydrogenBondLength,0.0,0.0 });

		SCESS::Molecule hydrogenMol({ H1, H2 });

		return hydrogenMol;
	}

	Molecule GenerateLithiumHydrideMolecule()
	{
		constexpr double LiHBondLength = 1.60 * ANGSTROM_2_BOHR;

		SCESS::Atom H("H", { 0.0,0.0,0.0 });
		SCESS::Atom Li("Li", { LiHBondLength,0.0,0.0 });

		SCESS::Molecule LiHMolecule({ H, Li });

		return LiHMolecule;
	}

	Molecule GenerateCOMolecule()
	{
		constexpr double COBondLength = 1.128 * ANGSTROM_2_BOHR;

		SCESS::Atom C("C", { 0.0,0.0,0.0 });
		SCESS::Atom O("O", { COBondLength,0.0,0.0 });

		SCESS::Molecule COMolecule({ C,O });

		return COMolecule;
	}

	Molecule GenerateMethaneMolecule()
	{
		SCESS::Atom H1("H", { -1.183771681898, -1.183771681898 ,  1.183771681898 });
		SCESS::Atom H2("H", { -1.183771681898 ,  1.183771681898, -1.183771681898 });
		SCESS::Atom H3("H", { 1.183771681898,   1.183771681898,   1.183771681898 });
		SCESS::Atom H4("H", { 1.183771681898 ,-1.183771681898, -1.183771681898 });
		SCESS::Atom C("C", { 0,0,0 });

		SCESS::Molecule methane({ H1, H2, H3, H4, C });

		return methane;
	}

	Molecule GenerateEtheneMolecule()
	{
		SCESS::Atom C1("C", Eigen::Vector3d{ -0.66754171391078,-0.57999587657014 ,-0.98766494519065 }*SCESS::ANGSTROM_2_BOHR);
		SCESS::Atom C2("C", Eigen::Vector3d{ -0.00280161842801,-0.00243411244334 ,-0.00414503592554 }*SCESS::ANGSTROM_2_BOHR);
		SCESS::Atom H1("H", Eigen::Vector3d{ -0.15803343371248,-1.06435755385651 ,-1.81247599918377 }*SCESS::ANGSTROM_2_BOHR);
		SCESS::Atom H2("H", Eigen::Vector3d{ -1.75092796040950,-0.59425022738565 ,-1.01193846765821 }*SCESS::ANGSTROM_2_BOHR);
		SCESS::Atom H3("H", Eigen::Vector3d{ 1.08058463979347 ,  0.01182021009784 ,  0.02012843810545 }*SCESS::ANGSTROM_2_BOHR);
		SCESS::Atom H4("H", Eigen::Vector3d{ -0.51230991333269,  0.48192756015780 ,  0.82066600985273 }*SCESS::ANGSTROM_2_BOHR);

		SCESS::Molecule ethene({ C1,C2,H1,H2,H3,H4 });

		return ethene;
	}

	Molecule GenerateAmmoniaMolecule()
	{
		Eigen::Vector3d NPos{ -6.6842683864,   1.9896777361, -0.0017590331 };
		NPos *= ANGSTROM_2_BOHR;
		Atom N("N", { 0,0,0 });
		Atom H1("H", Eigen::Vector3d{ -5.6844990792,   1.9942938061 ,  0.0072441956 }*ANGSTROM_2_BOHR - NPos);
		Atom H2("H", Eigen::Vector3d{ -7.0127564806 ,  2.9323453345 ,  0.0540465122 }*ANGSTROM_2_BOHR - NPos);
		Atom H3("H", Eigen::Vector3d{ -7.0127560539 ,  1.4847331233  , 0.7962183253 }*ANGSTROM_2_BOHR - NPos);

		Molecule NH3({ N,H1,H2,H3 });

		return NH3;
	}
}