#include<pch.h>
#include <iostream>
#include "Tests.h"
#include "Quantum/Molecule.h"
#include "Quantum/Integrals/IntegralsLibInt.h"
#include "Resources/Resources.h"

int main(int argc, char *argv[]) {
    using namespace SCESS;
    std::cout << "Schrodinger's Cat's Electronic System Solver; Terminal Edition\n    By Charlie Wheelwright (cspw2) \n    Meow!\n";
    if (!Res::AssetsFolder::Get().Init())
        return 1;

    Integrals::LibInt::InitLibInt();

    Tests::TestF12Integrals();
    //Tests::TestCCSD_F12();

    Integrals::LibInt::CloseLibInt();

    return 0;
}
