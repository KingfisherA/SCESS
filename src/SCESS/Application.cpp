#include <iostream>
#include "Tests.h"

int main(int argc, char* argv[])
{
    std::cout << "SCESS Command Line \n";

    SCESS::Tests::TestIntegrals();

    std::cin.get();
    return 0;
}
