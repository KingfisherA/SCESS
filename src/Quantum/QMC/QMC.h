#pragma once
#include <boost/circular_buffer.hpp>
#include <glm/vec2.hpp>

#include "Quantum/CI/ConfigurationInteraction.h"
#include "Quantum/HF/HartreeFock.h"

namespace GEng {
    class Thread;
}

namespace SCESS {
    using Excitor = CIDeterminant;
    struct RHFWavefunction;

    class QuantumMonteCarlo {
    public:
        //Full CI quantum monte carlo
        static double FCIQMC(const RHFWavefunction &rhfWfn,
                             double Vnn,
                             const GEng::Thread *thread = nullptr,
                             std::stringstream *console = nullptr,
                             boost::circular_buffer<glm::vec2> *projEnergyGraph = nullptr);

        // Unlinked coupled-cluster quantum monte carlo. The similarity transformed Hamiltonian is not
        // applied here.
        static double CCQMC(const RHFWavefunction &rhfWfn, double Vnn,
                            const GEng::Thread *thread = nullptr,
                            boost::circular_buffer<glm::vec2> *projEnergyGraph = nullptr);

        static std::vector<Excitor> GenerateAllExcitors(CIDeterminant reference, uint32_t nSpinOrbs,
                                                        uint32_t nElectrons);
    };
} // SCESS
