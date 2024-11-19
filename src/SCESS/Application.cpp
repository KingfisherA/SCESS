#include<pch.h>
#include <iostream>
#include "Tests.h"
#include "Quantum/Integrals/IntegralsLibInt.h"
#include "Resources/Resources.h"

int main(int argc, char *argv[]) {
    using namespace SCESS;
    std::cout << "SCESS Command Line \n";

    if (!Res::AssetsFolder::Get().Init())
        return 1;


    Integrals::LibInt::InitLibInt();
    //
    // //SCESS::Tests::TestRepulsionIntegrals();
    // //SCESS::Tests::TestF12Integrals();
    // //SCESS::Tests::TestDFT();
    //
    // Molecule waterMolecule = GenerateLithiumHydrideMolecule();
    // waterMolecule.MakeOrbitals("sto-3g");
    //
    // //TEITensor tei = Integrals::TwoElectronRepulsionTensor(waterMolecule);
    // TEITensor tei2 = Integrals::TwoElectronRepulsionTensorLibInt(waterMolecule, "sto-3g");
    //
    // // std::cout <<"SCESS tensor\n";
    // // Utilities::PrintTensor(tei);
    //
    // std::cout <<"LibInt Tensor\n";
    // Utilities::PrintTensor(tei2);
    //

    //Tests::TestFCIQMC();
    //Tests::TestCCQMC();

    Integrals::LibInt::CloseLibInt();

    return 0;
}
