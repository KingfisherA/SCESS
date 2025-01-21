#include <iostream>
#include "Quantum/Integrals/IntegralsF12.h"

#include "Tests.h"
#include "Quantum/Molecule.h"
#include "Quantum/Integrals/IntegralsLibInt.h"
#include "Resources/AssetFolder.h"

int main(int argc, char *argv[]) {
    using namespace SCESS;
    std::cout << "Schrodinger's Cat's Electronic System Solver; Terminal Edition\n    By Charlie Wheelwright (cspw2) \n    Meow!\n";
    if (!Res::AssetsFolder::Get().Init())
        return 1;

    Integrals::LibInt::InitLibInt();

    Tests::TestMP2_F12();
    //Tests::TestF12Integrals();
    //Tests::TestCCSD_F12();

    Integrals::LibInt::CloseLibInt();

    return 0;
}
