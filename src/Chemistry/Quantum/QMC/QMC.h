#pragma once
#include "Chemistry/Quantum/HF/HartreeFock.h"


namespace SCESS {
    struct RHFWavefunction;
    struct Excitor {
        public:
        std::array<int, 10> occIndices, virtIndices;
    };

    class QuantumMonteCarlo {
    public:
        //Full CI quantum monte carlo
        static double SimpleFCIQMC(const RHFWavefunction &rhfWfn,
                                   double Vnn, std::stringstream *console = nullptr);

        // Unlinked coupled-cluster quantum monte carlo. The similarity transformed Hamiltonian is not
        // applied here.
        static double CCQMC(const RHFWavefunction &rhfWfn);
    };
} // SCESS
