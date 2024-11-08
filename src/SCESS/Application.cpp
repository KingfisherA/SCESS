#include <iostream>
#include "Tests.h"
#include "Quantum/Integrals/Integrals.h"
#include "Resources/Resources.h"

int main(int argc, char *argv[]) {
    std::cout << "SCESS Command Line \n";

    using SCESS::Res::AssetsFolder;

    if (!AssetsFolder::Get().Init())
        return 1;

    SCESS::Integrals::InitLibInt();

    //SCESS::Tests::TestRepulsionIntegrals();
    SCESS::Tests::TestF12Integrals();
    //SCESS::Tests::TestCCQMC();
    // SCESS::Tests::TestDFT();

    SCESS::Integrals::CloseLibInt();

    return 0;
}
