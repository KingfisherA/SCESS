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
            (int) f12BasisSet.GetNoHFFuncs() - alphaOcc,
            alphaOcc
        };
    }
    //
    //
    struct MP2F12Result {
        double E_MP2;
        double E_F12_MP2;
    };
    //
    //The intermediates of the MP2-F12/3C approximation in the canonical OBS + CABS MO basis. The four-index
    //quantities are in physicists' notation over the occupied spatial orbitals.
    struct MP2F12Intermediates {
        //Fock matrix over the whole OBS + CABS space
        Tensor2 F;
        //Exchange contribution to the Fock matrix
        Tensor2 K;
        //V_ij^kl = <ij| r12^-1 Q12 f12 |kl>, Coulomb operator on ij and geminal on kl, indexed (i, j, k, l)
        TEITensor V;
        //X_kl^mn = <kl| f12 Q12 f12 |mn>
        TEITensor X;
        //B_kl^mn
        TEITensor B;
        //C_kl^ab, indexed (k, l, a, b) with a, b virtual
        TEITensor C;
    };
    //
    MP2F12Intermediates EvaluateMP2F12Intermediates(const Integrals::MolInts &mints,
                                                    const Integrals::F12::F12Integrals &f12Ints,
                                                    F12Idxs f12Idxs);
    //
    //Conventional closed-shell MP2 energy from physicists' notation MO integrals
    double EvaluateMP2Energy(const TEITensor &eri, const Eigen::VectorXd &epsilon, F12Idxs f12Idxs);
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
                               const Integrals::MolInts &mints,
                               const Integrals::F12::F12Integrals &f12Ints,
                               F12Idxs f12Idxs,
                               Tensor2 *F_pq_out = nullptr,
                               TEITensor *V_ijxy_out = nullptr,
                               TEITensor *C_xyab_out = nullptr,
                               TEITensor *X_xywz_out = nullptr,
                               TEITensor *B_xywz_out = nullptr);
    //
    Tensor4 EvaluateW(const Tensor4 &f12Int, F12Idxs f12Idxs);
    //
    //Intermediates for the MEG Hamiltonian matrix elements, both indexed (p, q, i, j) for all OBS p, q and
    //occupied i, j: V_pq^ij = <pq| r12^-1 Q12 f12 |ij> and h_pq^ij = <pq| (h1 + h2) Q12 f12 |ij>
    void EvaluateAdditionalMP2Intermediates(
        const Integrals::MolInts &configIntegrals,
        const Integrals::F12::F12Integrals &f12Ints,
        F12Idxs f12Idxs,
        TEITensor *V_pqxy_out = nullptr,
        TEITensor *h_pqxy_out = nullptr
    );
    //
    double CalculateMP2_F12_FixedEnergyElectronPair(long electronI, long electronJ,
                                                    const Tensor2 &F,
                                                    const TEITensor &V,
                                                    const TEITensor &X,
                                                    const TEITensor &B);
    //
    //
}
