#pragma once
#include "Quantum/CI/CI.h"
#include "Quantum/F12/F12.h"

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
    //
    //
    double EvaluateGeminalElement(SlaterDet Dj,
                                  const Tensor4 &Vij,
                                  const Tensor4 &Vab);
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
    //
    //
    double CCQMC_MP2_F12(const HF::HF_Wfn &refWfn,
        const Integrals::MolInts &mints,
                         const Integrals::F12::F12Integrals &f12Ints,
        F12::F12Idxs f12Idxs);
    //
    //
}
