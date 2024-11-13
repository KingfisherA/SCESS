#pragma once
#include <boost/circular_buffer.hpp>
#include <glm/vec2.hpp>

#include  "Quantum/CI/ConfigurationInteraction.h"
#include  "Quantum/HF/HartreeFock.h"


namespace SCESS {
    using Excitor = CIDeterminant;
    struct RHFWavefunction;

    class QuantumMonteCarlo {
    public:
        //Full CI quantum monte carlo
        static double SimpleFCIQMC(const RHFWavefunction &rhfWfn,
                                   double Vnn, std::stringstream *console = nullptr,
                                   boost::circular_buffer<glm::vec2> *projEnergyGraph = nullptr);

        // Unlinked coupled-cluster quantum monte carlo. The similarity transformed Hamiltonian is not
        // applied here.
        static double CCQMC(const RHFWavefunction &rhfWfn, double Vnn,
                            boost::circular_buffer<glm::vec2> *projEnergyGraph = nullptr);

        static std::vector<Excitor> GenerateAllExcitors(CIDeterminant reference, uint32_t nSpinOrbs,
                                                        uint32_t nElectrons);
    };
} // SCESS
