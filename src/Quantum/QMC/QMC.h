#pragma once
#include <functional>
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
    //Called after every report cycle with the cycle number and the projected energy. It is called from the thread
    //running the calculation.
    using ProjEnergyCallback = std::function<void(size_t reportCycle, double projectedEnergy)>;
    //
    //
    namespace FCI {
        //Full CI quantum monte carlo
        double FCIQMC(const HF_Wfn &rhfWfn,
                      const Integrals::MolInts &mints,
                      double Vnn,
                      const DH::Thread *thread = nullptr,
                      std::stringstream *console = nullptr,
                      const ProjEnergyCallback &onReportCycle = {});
    }

    //
    //
    namespace CC {
        // Unlinked coupled-cluster quantum monte carlo. The similarity transformed Hamiltonian is not
        // applied here.
        double CCQMC(const HF_Wfn &rhfWfn,
                     const Integrals::MolInts &mints,
                     const DH::Thread *thread = nullptr,
                     const ProjEnergyCallback &onReportCycle = {});
        //
        //
        uint32_t CountNumElecs(OrbBitStr bitStr, const uint32_t nSpatialOrbitals);
        //
        // Randomly excite one electron of I into an unoccupied orbital of the same spin. pExcit is multiplied by
        // the probability of generating the returned determinant. I is returned unchanged if no excitation exists.
        SlaterDet RandomSingleExcite(SlaterDet I,
                                     boost::random::mt19937 &randGenerator,
                                     uint32_t nSpatialOrbitals,
                                     double &pExcit);
        // Randomly excite a pair of electrons of I into unoccupied orbitals of the same spins. pExcit is
        // multiplied by the probability of generating the returned determinant. I is returned unchanged if no
        // excitation exists.
        SlaterDet RandomDoubleExcite(SlaterDet I,
                                     boost::random::mt19937 &randGenerator,
                                     uint32_t nSpatialOrbitals,
                                     double &pExcit);
        //
        // Apply an excitor (stored as its difference to the reference determinant) to det. The excitor is the
        // operator a+_p1 a+_p2 ... a_m2 a_m1 with spin orbitals ordered with every alpha before every beta orbital,
        // the convention of the Slater-Condon phases in CI. sign is multiplied by the resulting fermionic phase.
        // Returns false if the result vanishes as the excitor overlaps with excitations already in det.
        bool ApplyExcitor(Excitor excitor,
                          SlaterDet reference,
                          SlaterDet &det,
                          int &sign,
                          uint32_t nSpatialOrbitals);
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
                                                  uint32_t nAlphaElectrons,
                                                  uint32_t nBetaElectrons);
        //
        // Total number of excips, sum_i |N_i|
        int64_t CountExcips(const std::vector<int64_t> &excitorWalkers);
        //
        // Cumulative distribution of the cluster size s = 0 ... maxClusterSize with p(s) = 1/2^(s+1); the
        // remaining probability of the infinite tail is given to the largest cluster size.
        std::vector<double> GenerateCumulativeClusterProb(uint32_t maxClusterSize);
        // Sample a cluster size from the cumulative distribution with a uniform random number in [0, 1)
        uint32_t SampleClusterSize(double randomNumber,
                                   const std::vector<double> &cumulativeClusterProb,
                                   double &pSize);
        //
        // Projected correlation energy sum_j <D0|H|Dj> c_j / N_0 where c_j collects the excitor amplitude and the
        // disconnected product of two singles amplitudes on every single and double excitation Dj
        double CalculateProjectedEnergy(const std::vector<Excitor> &excitors,
                                        const std::vector<int64_t> &excitorWalkers,
                                        const Eigen::MatrixXd &H,
                                        SlaterDet reference,
                                        uint32_t nSpatialOrbitals);
        //
        void PrintCCQMCWavefunction(const std::vector<int64_t> &excitorWalkers);
        //
        //
    }

    //
} // SCESS
