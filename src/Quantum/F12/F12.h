#pragma once
#include "Quantum/HF/HF.h"
#include <unsupported/Eigen/CXX11/Tensor>

#include "Quantum/Basis/BasisSet.h"

namespace SCESS::HF {
    struct HF_Wfn;
}

namespace SCESS::Integrals {
    struct MolInts;
};

namespace SCESS::Integrals::F12 {
    struct F12Integrals;
}


//
//
namespace SCESS::F12 {
    //
    //
    using HF_Wfn = HF::HF_Wfn;
    using Tensor2 = Eigen::Tensor<double, 2>;
    using Tensor4 = Eigen::Tensor<double, 4>;
    using TEITensor = Eigen::Tensor<double, 4>;
    //
    //
    struct F12Idxs {
        //Number of orbitals in the ordinary basis set. Also the index for the beginning of the CABS
        int nOBS;
        //Number of orbitals in the complementary auxiliary set
        int nCABS;
        //Only one occupied spatial MO in hydrogen
        int nVirt;
        //The number of spatial orbital occupied
        int nOcc;
    };
    //
    inline F12Idxs GenerateF12Idxs(const F12BasisSet &f12BasisSet, const HF_Wfn &wfn) {
        int alphaOcc = std::distance(wfn.Occ_Alpha.begin(),
                                     std::ranges::find(wfn.Occ_Alpha, 0));
        return F12Idxs{
            (int) f12BasisSet.GetNoHFFuncs(),
            (int) f12BasisSet.GetNoCABSFuncs(),
            alphaOcc,
            (int) f12BasisSet.GetNoHFFuncs() - alphaOcc
        };
    }
    //
    //
    struct MP2F12Result {
        double E_MP2;
        double E_F12_MP2;
    };
    //
    //
    double EvaluateCABSSinglesCorrection(const Tensor2 &F,
                                         F12Idxs f12Idxs,
                                         Tensor2 *t1Amps = nullptr);
    //
    //
    MP2F12Result MP2_F12_Optimise(const HF_Wfn &refWfn,
                                  const Integrals::MolInts &configIntegrals,
                                  const Integrals::F12::F12Integrals &f12Ints,
                                  F12Idxs f12Idxs,
                                  Tensor2 *F_out = nullptr,
                                  TEITensor *V_out = nullptr,
                                  TEITensor *X_out = nullptr,
                                  TEITensor *B_out = nullptr);
    //
    MP2F12Result MP2_F12_Fixed(const HF_Wfn &refWfn,
                               const Integrals::MolInts &configIntegrals,
                               const Integrals::F12::F12Integrals &f12Ints,
                               F12Idxs f12Idxs,
                               Tensor2 *F_pq_out = nullptr,
                               TEITensor *V_ijxy_out = nullptr,
                               TEITensor *C_xyab_out = nullptr,
                               TEITensor *X_xywz_out = nullptr,
                               TEITensor *B_xywz_out = nullptr);
    //
    Tensor4 EvaluateW(const Tensor4 &f12Int, F12Idxs f12Idxs);
    Tensor4 EvaluateVab(const Tensor4 &f12Repulsion,
                        const Tensor4 &f12,
                        const Tensor4 &eri,
                        F12Idxs f12Idxs);
    //
    double CalculateMP2_F12_FixedEnergyElectronPair(long electronI, long electronJ,
                                                    const Tensor2 &F,
                                                    const TEITensor &V,
                                                    const TEITensor &X,
                                                    const TEITensor &B);
    //
    //
}
