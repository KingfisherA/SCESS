#pragma once
//
#include "Quantum/CI/CI.h"
#include "Quantum/F12/F12.h"
//

namespace SCESS::HF {
    struct HF_Wfn;
}

namespace SCESS::Integrals {
    struct MolInts;
};

namespace SCESS::Integrals::F12 {
    struct F12Integrals;
}

namespace SCESS::F12_QMC {
    using Tensor2 = Eigen::Tensor<double, 2>;
    using Tensor4 = Eigen::Tensor<double, 4>;
    using TEITensor = Eigen::Tensor<double, 4>;
    using SlaterDet = CI::SlaterDet;
    using Excitor = SlaterDet;
    using OrbBitStr = CI::OrbBitStr;
    using ElectronPair = std::pair<int, int>;
    //
    //
    double EvaluateGeminalElement(SlaterDet Dj,
                                  ElectronPair gemPair,
                                  const Tensor4 &V_pqxy,
                                  const Tensor4 &h_apxy,
                                  uint32_t nAlphaElectrons,
                                  uint32_t nBetaElectrons,
                                  uint32_t nOcc,
                                  uint32_t numOfAlphaOrbitals,
                                  uint32_t numOfBetaOrbitals);
    inline Excitor MakeGeminalExcitor() {
        Excitor ex = {0, 0};
        ex.GetAlphaOrbitals() |= ((unsigned long) (1) << 63);
        ex.GetBetaOrbitals() |= ((unsigned long) (1) << 63);
        return ex;
    }
    inline bool IsGeminalExcitor(Excitor ex) {
        return (ex.GetAlphaOrbitals() & ((unsigned long) (1) << 63)) && (
                   ex.GetBetaOrbitals() & ((unsigned long) (1) << 63));
    }
    ElectronPair GetGeminalElectronPair(Excitor gemExcitor,
                                        uint32_t nAlphaOrbitals,
                                        uint32_t nBetaOrbitals);
    //
    //
    double CCSD_QMC_MP2_F12(const HF::HF_Wfn &refWfn,
                            const Integrals::MolInts &mints,
                            const Integrals::F12::F12Integrals &f12Ints,
                            F12::F12Idxs f12Idxs);
    //
    //
}
