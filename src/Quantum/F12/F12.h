#pragma once
#include "Quantum/HF/HF.h"
#include <unsupported/Eigen/CXX11/Tensor>

#include "Quantum/Basis/BasisSet.h"
#include "Quantum/Integrals/Integrals.h"
#include "Quantum/Integrals/IntegralsF12.h"

namespace SCESS::HF {
    struct HF_Wfn;
}

namespace DH {
    class Thread;
}

namespace SCESS {
    class Molecule;
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
    using F12MOIntegrals = Integrals::F12::F12MOIntegrals;
    //
    //Everything a F12 calculation needs once the reference has been obtained and the CABS canonicalised
    struct F12Calculation {
        F12BasisSet basisSet;
        //Closed-shell (UHF with equal alpha and beta occupations) reference in the OBS
        HF_Wfn hfWfn;
        F12Idxs idxs;
        //The integral blocks in the canonical OBS + CABS MO basis, physicists' notation
        F12MOIntegrals ints;
        double E_CABSSingles = 0.0;
    };
    //
    /**
     * Builds the OBS + CABS basis, evaluates the integrals, solves for the reference in the OBS, constructs and
     * canonicalises the CABS and evaluates the CABS singles correction. Shared by every F12 method.
     * @param obsName Orbital basis set
     * @param cabsName Basis set from which the CABS is built (the OBS is projected out of OBS + cabsName)
     * @param buildQMCIntegrals Also build F12MOIntegrals::eri_oxpq, which only CCSD-QMC[MP2-F12] needs
     */
    F12Calculation PrepareF12Calculation(const std::string &obsName,
                                         const std::string &cabsName,
                                         const Molecule &mol,
                                         const DH::Thread *thread = nullptr,
                                         std::stringstream *console = nullptr,
                                         bool buildQMCIntegrals = false);
    //
    /**
     * Extracts the F12MOIntegrals blocks from full (N^4) integral tensors which are already in the canonical
     * OBS + CABS MO basis and physicists' notation, e.g. integrals imported from MRCC
     * @param F Fock matrix in the same MO basis
     */
    F12MOIntegrals ExtractF12MOIntegrals(const Integrals::MolInts &moInts,
                                         const Integrals::F12::F12Integrals &moF12Ints,
                                         const Eigen::MatrixXd &F,
                                         F12Idxs f12Idxs);
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
    MP2F12Intermediates EvaluateMP2F12Intermediates(const F12MOIntegrals &ints, F12Idxs f12Idxs);
    //
    //Conventional closed-shell MP2 energy from the <ij|pq> MO integrals (physicists' notation)
    double EvaluateMP2Energy(const TEITensor &eri, const Eigen::VectorXd &epsilon, F12Idxs f12Idxs);
    //
    //
    double EvaluateCABSSinglesCorrection(const Tensor2 &F,
                                         F12Idxs f12Idxs,
                                         Tensor2 *t1Amps = nullptr);
    //
    //
    MP2F12Result MP2_F12_Optimise(const HF_Wfn &refWfn,
                                  const F12MOIntegrals &ints,
                                  F12Idxs f12Idxs,
                                  Tensor2 *F_out = nullptr,
                                  TEITensor *V_out = nullptr,
                                  TEITensor *X_out = nullptr,
                                  TEITensor *B_out = nullptr);
    //
    MP2F12Result MP2_F12_Fixed(const HF_Wfn &refWfn,
                               const F12MOIntegrals &ints,
                               F12Idxs f12Idxs,
                               Tensor2 *F_pq_out = nullptr,
                               TEITensor *V_ijxy_out = nullptr,
                               TEITensor *C_xyab_out = nullptr,
                               TEITensor *X_xywz_out = nullptr,
                               TEITensor *B_xywz_out = nullptr);
    //
    //W_ab^ij = <ab| f12 (3/8 + 1/8 P_ij) |ij> over the virtual orbitals of the OBS + CABS, indexed (a, b, i, j).
    //f12 holds the <ij|pq> integrals (F12MOIntegrals::f12)
    Tensor4 EvaluateW(const Tensor4 &f12, F12Idxs f12Idxs);
    //
    //Intermediates for the MEG Hamiltonian matrix elements, both indexed (p, q, i, j) for all OBS p, q and
    //occupied i, j: V_pq^ij = <pq| r12^-1 Q12 f12 |ij> and h_pq^ij = <pq| (h1 + h2) Q12 f12 |ij>.
    //Needs F12MOIntegrals::eri_oxpq (see PrepareF12Calculation)
    void EvaluateAdditionalMP2Intermediates(
        const F12MOIntegrals &ints,
        F12Idxs f12Idxs,
        TEITensor *V_pqxy_out = nullptr,
        TEITensor *h_pqxy_out = nullptr
    );
    //
    struct CCSD_F12bResult {
        //Conventional MP2 correlation energy in the OBS (for reference)
        double E_MP2 = 0.0;
        //CCSD correlation energy including the coupling of the doubles to the geminals
        double E_CCSD = 0.0;
        //The part of E_CCSD due to the coupling of the doubles to the geminals
        double E_Coupling = 0.0;
        //Energy of the fixed amplitude geminals, 2 c V + c (B - (e_i + e_j) X) c
        double E_F12 = 0.0;
        bool converged = false;

        [[nodiscard]] double CorrelationEnergy() const { return E_CCSD + E_F12; }
    };
    //
    /**
     * Closed-shell CCSD-F12b/3C(FIX) with the CABS approximation and the SP (fixed amplitude) geminals
     * c_ij^kl = 3/8 delta_ik delta_jl + 1/8 delta_il delta_jk.
     *
     * The doubles residual gains the geminal coupling sum_kl c_ij^kl (C_kl^ab + V_ab^kl) and the energy is
     *   E = E_CCSD(t1, t2) + sum_ijkl (2 c_ij^kl - c_ij^lk) sum_ab C_kl^ab t_ij^ab + E_MP2-F12(3C, FIX)
     * where, as in CCSD-F12x, V_ab^kl = <ab| r12^-1 f12 |kl> - sum_rs <ab| r12^-1 |rs><rs| f12 |kl> has no CABS
     * summations. Couplings of the geminals to the singles are neglected.
     * @param calc A calculation prepared with PrepareF12Calculation
     */
    CCSD_F12bResult CCSD_F12b(const F12Calculation &calc,
                              const DH::Thread *thread = nullptr,
                              std::stringstream *console = nullptr);
    //
    double CalculateMP2_F12_FixedEnergyElectronPair(long electronI, long electronJ,
                                                    const Tensor2 &F,
                                                    const TEITensor &V,
                                                    const TEITensor &X,
                                                    const TEITensor &B);
    //
    //
}
