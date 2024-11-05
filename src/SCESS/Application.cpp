#include <iostream>
#include "Tests.h"
#include "Quantum/Integrals/Integrals.h"

int main(int argc, char* argv[])
{
    std::cout << "SCESS Command Line \n";

    SCESS::Integrals::InitLibInt();

    SCESS::Tests::TestIntegrals();

    SCESS::Integrals::CloseLibInt();

    return 0;
}
