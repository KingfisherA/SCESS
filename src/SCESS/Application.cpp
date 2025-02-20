#include <iostream>
#include<string>
//
#include "Tests.h"
#include "Quantum/Integrals/IntegralsLibInt.h"
#include "Resources/AssetFolder.h"


//
//
int main(int argc, char *argv[]) {
    using namespace SCESS;
    std::cout <<
            "Schrodinger's Cat's Electronic System Solver; Terminal Edition\n    By Charlie Wheelwright (cspw2) \n    Meow!\n\n";
    if (!Res::AssetsFolder::Get().Init())
        return 1;

    Integrals::LibInt::InitLibInt();

    std::cout << "\n"
            "--- Command Line Arguments ---\n"
            "\t- \"f12ints\" = Calculate a test set of f12 integrals and transform the basis set\n"
            "\t- \"mp2f12-opt\" = Perform an MP2-F12 optimised amplitude calculation on a hydrogen molecule.\n"
            "\t- \"mp2f12-fixed\" = Perform an MP2-F12 fixed amplitude calculation on a hydrogen molecule.\n"
            "\t- \"ccsd\" = Perform a CCSD calculation on hydrogen\n"
            "\t- \"dft\" = Perform a test KS-DFT calculation on water with the LDA exchange and VWD correlation functional.\n"
            "\t- \"fci\" = Perform a test full configuration interaction calculation on hydrogen.\n"
            "\t- \"fciqmc\" = Perform a test Full-CI Quantum Monte Carlo procedure on hydrogen.\n"
            "\t- \"ccqmc-mp2-f12\" = Perform a CCSD Quantum Monte Carlo with a second-order geminal excitor operator T2' in the T expansion (new theory). \n";

    //Tests::TestMP2_F12_Fixed();
    //Tests::TestFCIQMC();
    Tests::TestCCQMC_MP2_F12();
    //Tests::TestFullCI();
    //Tests::TestMolecularInts();
    //Tests::TestF12Integrals();
    //Tests::TestTensorEinsum();

    if (argc > 1) {
        std::string arg1 = argv[1];
        if (arg1 == "f12ints")
            Tests::TestF12Integrals();
        else if (arg1 == "mp2f12-opt")
            Tests::TestMP2_F12_Optimise();
        else if (arg1 == "mp2f12-fixed")
            Tests::TestMP2_F12_Fixed();
        else if (arg1 == "ccsd")
            Tests::TestCoupledCluster();
        else if (arg1 == "dft")
            Tests::TestDFT();
        else if (arg1 == "fci")
            Tests::TestFullCI();
        else if (arg1 == "fciqmc")
            Tests::TestFCIQMC();
        else if (arg1 == "ccqmc-mp2-f12")
            Tests::TestCCQMC_MP2_F12();
        else
            std::cout << "Command \"" << arg1 << "\" not recognised\n";
    } else {
        std::cout << "--- No command given ---\n";
    }

    Integrals::LibInt::CloseLibInt();

    return 0;
}
