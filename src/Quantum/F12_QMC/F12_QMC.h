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
    // The Hamiltonian matrix element <Dj|H|MEG> between the determinant Dj and the multi-electron geminal (MEG)
    // where the alpha orbital gemPair.first and the beta orbital gemPair.second of the reference are replaced by the
    // fixed amplitude geminal. Only determinants which differ from the reference within the geminal pair couple.
    double EvaluateGeminalElement(SlaterDet Dj,
                                  ElectronPair gemPair,
                                  const Tensor4 &V_pqxy,
                                  const Tensor4 &h_pqxy,
                                  uint32_t nAlphaElectrons,
                                  uint32_t nBetaElectrons,
                                  uint32_t numOfAlphaOrbitals,
                                  uint32_t numOfBetaOrbitals);
    inline Excitor MakeGeminalExcitor() {
        Excitor ex = {0, 0};
        ex.GetAlphaOrbitals() |= (1ULL << 63);
        ex.GetBetaOrbitals() |= (1ULL << 63);
        return ex;
    }
    inline bool IsGeminalExcitor(Excitor ex) {
        return (ex.GetAlphaOrbitals() & (1ULL << 63)) && (
                   ex.GetBetaOrbitals() & (1ULL << 63));
    }
    ElectronPair GetGeminalElectronPair(Excitor gemExcitor,
                                        uint32_t nAlphaOrbitals,
                                        uint32_t nBetaOrbitals);
    //
    //
    //ints must include F12MOIntegrals::eri_oxpq (PrepareF12Calculation with buildQMCIntegrals)
    double CCSD_QMC_MP2_F12(const HF::HF_Wfn &refWfn,
                            const F12::F12MOIntegrals &ints,
                            F12::F12Idxs f12Idxs);
    //
    //
}
