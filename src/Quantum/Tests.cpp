#include<pch.h>

#include"Basis/BasisUtils.h"
#include"Molecule.h"
#include"Integrals/IntegralHelper.h"
#include"Basis/Gaussian.h"
#include"HF/HartreeFock.h"


void PrintMoleculeGaussians(const SCESS::Molecule& mol)
{
    std::stringstream ss;
    for (const SCESS::GaussianCollection& collection : mol.GetOrbitals())
    {
        for (const SCESS::Gaussian& g : collection.GetGaussians())
        {
            ss << "-----\n";
            ss << "Position: " << g.mPosition.x() << ", " << g.mPosition.y() << ", " << g.mPosition.z() << '\n';
            ss << "Angular: " << g.mAngular.x() << ", " << g.mAngular.y() << ", " << g.mAngular.z() << '\n';
            ss << "Exponent: " << g.mExponent << " Coefficient: " << g.mCoeff << " Normalisation: " << g.mNorm << '\n';
            ss << "-----\n";
        }
    }
    puts(ss.str().c_str());
}

void CalculateHydrogenMoleculeRHF()
{
    constexpr double hydrogenBondLength = 0.7122300374 * SCESS::ANGSTROM_2_BOHR;

    SCESS::Atom H1("H", { 0.0,0.0,0.0 });
    SCESS::Atom H2("H", { hydrogenBondLength,0.0,0.0 });

    SCESS::Molecule hydrogenMol({ H1, H2 });
    hydrogenMol.MakeOrbitals();
  
    SCESS::HF::ReducedHF(hydrogenMol);
}

//Should be -74.9659012
void CalculateWaterMoleculeRHF()
{
    SCESS::Atom H1("H", { 0., 1.43233673, -0.96104039 });
    SCESS::Atom H2("H", { 0.0,  -1.43233673, -0.96104039 });
    SCESS::Atom O  ("O", { 0.0, 0.0, 0.24026010});

    SCESS::Molecule water({H1, H2, O});
    water.MakeOrbitals();

    PrintMoleculeGaussians(water);

    SCESS::HF::ReducedHF(water);
}

//Should be -39.7268637
void CalculateMethaneMoleculeRHF()
{
    SCESS::Atom H1("H", Eigen::Vector3d{ -9.3880270305,   1.9905798101, - 0.0000169136 }*SCESS::ANGSTROM_2_BOHR);
    SCESS::Atom H2("H", Eigen::Vector3d{ -8.3050312953 ,  1.9905766754 ,  0.0000244554 }*SCESS::ANGSTROM_2_BOHR);
    SCESS::Atom H3("H", Eigen::Vector3d{ -9.7489816731,   2.8099628719 ,  0.6092400251 }*SCESS::ANGSTROM_2_BOHR);
    SCESS::Atom H4("H", Eigen::Vector3d{ -9.7489800239 ,  1.0532708184 ,  0.4049933956 }*SCESS::ANGSTROM_2_BOHR);
    SCESS::Atom C ("C", Eigen::Vector3d{ -9.3880270305,   1.9905798101, - 0.0000169136 }*SCESS::ANGSTROM_2_BOHR);

    SCESS::Molecule methane({ H1, H2, H3, H4, C });
    methane.MakeOrbitals();
 
    SCESS::HF::ReducedHF(methane);

}

//Should be -77.0739544
void CalculateEtheneRHF()
{
    SCESS::Atom C1("C", Eigen::Vector3d{-0.6530176758, 0.0000000000, 0.0000000000}*SCESS::ANGSTROM_2_BOHR);
    SCESS::Atom C2("C", Eigen::Vector3d{ 0.6530176758, 0.0000000000, 0.0000000000 }*SCESS::ANGSTROM_2_BOHR);
    SCESS::Atom H1("H", Eigen::Vector3d{ -1.2288875372, -0.9156191261, 0.0000000000 }*SCESS::ANGSTROM_2_BOHR);
    SCESS::Atom H2("H", Eigen::Vector3d{ -1.2288875372, 0.9156191261, 0.0000000000 }*SCESS::ANGSTROM_2_BOHR);
    SCESS::Atom H3("H", Eigen::Vector3d{ 1.2288875372, 0.9156191261, 0.0000000000 }*SCESS::ANGSTROM_2_BOHR);
    SCESS::Atom H4("H", Eigen::Vector3d{ 1.2288875372, -0.9156191261, 0.0000000000 }*SCESS::ANGSTROM_2_BOHR);

    SCESS::Molecule ethene({ C1,C2,H1,H2,H3,H4 });
    ethene.MakeOrbitals();

    SCESS::HF::ReducedHF(ethene);
}

void CalculateHeliumRHF()
{
    SCESS::Atom He("He", Eigen::Vector3d{ 0,0,0 });

    SCESS::Molecule HeM({ He });
    HeM.MakeOrbitals();

    SCESS::HF::ReducedHF(HeM);
}