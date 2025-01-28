#include <iostream>
#include "Quantum/Integrals/IntegralsF12.h"

#include<string>
#include "Tests.h"
#include "Quantum/Molecule.h"
#include "Quantum/Integrals/IntegralsLibInt.h"
#include "Resources/AssetFolder.h"


int main(int argc, char *argv[]) {
    using namespace SCESS;
    std::cout <<
            "Schrodinger's Cat's Electronic System Solver; Terminal Edition\n    By Charlie Wheelwright (cspw2) \n    Meow!\n";
    if (!Res::AssetsFolder::Get().Init())
        return 1;

    Integrals::LibInt::InitLibInt();

    //Tests::TestTensorEinsum();
    //Tests::TestCCQMC_F12();
    Tests::TestMP2_F12();

    if (argc > 1) {
        std::string arg1 = argv[1];
        if (arg1 == "f12ints")
            Tests::TestF12Integrals();
        else if (arg1 == "mp2f12")
            Tests::TestMP2_F12();
        else if (arg1 == "ccsd")
            Tests::TestCoupledCluster();
        else if (arg1 == "dft")
            Tests::TestDFT();
        else if (arg1 == "fci")
            Tests::TestFullCI();
        else
            std::cout << "Command \"" << arg1 << "\" not recognised\n";
    } else {
        std::cout << "No command given.\n";
    }

    Integrals::LibInt::CloseLibInt();

    return 0;
}
