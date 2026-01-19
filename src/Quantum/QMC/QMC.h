#pragma once
#include <boost/circular_buffer.hpp>
#include <glm/vec2.hpp>
#include <boost/random/mersenne_twister.hpp>
//
#include "Quantum/CI/CI.h"
//


namespace DH {
    class Thread;
}

namespace SCESS::HF {
    struct HF_Wfn;
}

namespace SCESS::CI {
    class SlaterDet;
}

namespace SCESS::Integrals {
    struct MolInts;
}

//
//
namespace SCESS::QMC {
    using SlaterDet = CI::SlaterDet;
    using Excitor = CI::SlaterDet;
    using OrbBitStr = CI::OrbBitStr;
    using HF_Wfn = HF::HF_Wfn;
    using TEITensor = Eigen::Tensor<double, 4>;
    //
    //
    namespace FCI {
        //Full CI quantum monte carlo
        double FCIQMC(const HF_Wfn &rhfWfn,
                      const Integrals::MolInts &mints,
                      double Vnn,
                      const DH::Thread *thread = nullptr,
                      std::stringstream *console = nullptr,
                      boost::circular_buffer<glm::vec2> *projEnergyGraph = nullptr);
    }

    //
    //
    namespace CC {
        // Unlinked coupled-cluster quantum monte carlo. The similarity transformed Hamiltonian is not
        // applied here.
        double CCQMC(const HF_Wfn &rhfWfn,
                     const Integrals::MolInts &mints,
                     const DH::Thread *thread = nullptr,
                     boost::circular_buffer<glm::vec2> *projEnergyGraph = nullptr);
        //
        //
        void CountOperations(Excitor excitor,
                             const uint32_t nAnnihilationBits,
                             const uint32_t nCreationBits,
                             uint32_t &numAnnihAlpha,
                             uint32_t &numCreationAlpha,
                             uint32_t &numAnnihBeta,
                             uint32_t &numCreationBeta);
        //
        uint32_t CountNumElecs(OrbBitStr bitStr, const uint32_t nSpatialOrbitals);
        //
        void DestroyCreateElec(OrbBitStr &sourceBitStr,
                               OrbBitStr &destBitStr,
                               uint32_t orbitalToDelete,
                               uint32_t orbitalToCreate,
                               const uint32_t nAnnihilationBits,
                               const uint32_t nCreationBits);
        //
        SlaterDet RandomSingleExcite(SlaterDet I,
                                     boost::random::mt19937 &randGenerator,
                                     const uint32_t nSpatialOrbitals,
                                     const uint32_t nSpinOrbs,
                                     const uint32_t nElectrons,
                                     const uint32_t nAnnihilationBits,
                                     const uint32_t nCreationBits,
                                     double &pExcit);
        SlaterDet RandomDoubleExcite(SlaterDet I,
                                     boost::random::mt19937 &randGenerator,
                                     const uint32_t nSpatialOrbitals,
                                     const uint32_t nSpinOrbs,
                                     const uint32_t nElectrons,
                                     const uint32_t nAnnihilationBits,
                                     const uint32_t nCreationBits,
                                     double &pExcit);
        //
        uint32_t FindExcitorIdx(Excitor ex, const std::vector<Excitor> &excitors);
        //
        std::vector<Excitor> EnumerateExcitors(SlaterDet reference,
                                               uint32_t nSpinOrbs,
                                               uint32_t nAlphaElectrons,
                                               uint32_t nBetaElectrons,
                                               uint32_t excitationLevel);
        std::vector<Excitor> EnumerateAllExcitors(SlaterDet reference,
                                                  uint32_t nSpinOrbs,
                                                  uint32_t nElectrons);
        //
        void UpdateCumulativePop(std::vector<int64_t> &cumulativeWalkers,
                                 std::vector<int64_t> &excitorWalkers,
                                 int64_t &nExcips);
        //
        std::vector<double> GenerateCumulativeClusterProb(uint32_t nElectrons);
        //
        void PrintCCQMCWavefunction(const std::vector<int64_t> &excitorWalkers);
        //
        //
    }

    //
} // SCESS
