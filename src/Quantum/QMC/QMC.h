#pragma once
#include  "Quantum/CI/ConfigurationInteraction.h"
#include  "Quantum/HF/HartreeFock.h"


namespace SCESS {
    using Excitor = CIDeterminant;
    struct RHFWavefunction;

    class QuantumMonteCarlo {
    public:
        //Full CI quantum monte carlo
        static double SimpleFCIQMC(const RHFWavefunction &rhfWfn,
                                   double Vnn, std::stringstream *console = nullptr);

        // Unlinked coupled-cluster quantum monte carlo. The similarity transformed Hamiltonian is not
        // applied here.
        static double CCQMC(const RHFWavefunction &rhfWfn, double Vnn);

        static std::vector<Excitor> GenerateAllExcitors(CIDeterminant reference, uint32_t nSpinOrbs, uint32_t nElectrons);
    };
} // SCESS
