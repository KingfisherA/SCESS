#include <Quantum/Q_PCH.h>
#include "F12.h"
//
#include "Dahlia/Tools/Threading.h"
//
#include "Quantum/CC/CC.h"
#include "Quantum/HF/HF.h"
#include "Quantum/Integrals/IntegralTransforms.h"
#include "Quantum/Integrals/Integrals.h"
#include "Quantum/Integrals/IntegralsF12.h"
#include "Tensor/TensorEinsum.h"
#include "Tensor/TensorOperations.h"
#include "Tensor/TensorSlice.h"

//
//
namespace SCESS::F12 {
    using namespace Eigen;
    //
    //
    namespace {
        enum {
            //Occupied orbitals
            I, J, K, L,
            //OBS virtual orbitals
            A, B,
            //All OBS orbitals
            R, S
        };
        //
        //Contracts a pair quantity X_kl^ab, indexed (k, l, a, b), with the fixed SP geminal amplitudes:
        //sum_kl c_ij^kl X_kl^ab = 3/8 X_ij^ab + 1/8 X_ji^ab, indexed (i, j, a, b)
        Tensor4 ContractFixedAmplitudes(const Tensor4 &X_klab) {
            Tensor4 result = (3.0 / 8.0) * X_klab;
            result += (1.0 / 8.0) * Einsum<Indices<J, I, A, B>, Indices<I, J, A, B> >(X_klab);
            return result;
        }
    }
    //
    //
    CCSD_F12bResult CCSD_F12b(const F12Calculation &calc, const DH::Thread *thread, std::stringstream *console) {
        // * The reference is a closed-shell determinant in canonical orbitals
        // * The integrals have been canonicalised and transformed into the OBS + CABS MO basis set and are in
        //   physicists' notation
        const F12Idxs f12Idxs = calc.idxs;
        const Index nOcc = f12Idxs.nOcc;
        const Index nVirt = f12Idxs.nVirt;
        const Index nOBS = f12Idxs.nOBS;

        TensorSliceIndices idxs = {
            {I, {0, nOcc}},
            {J, {0, nOcc}},
            {K, {0, nOcc}},
            {L, {0, nOcc}},
            {A, {nOcc, nVirt}},
            {B, {nOcc, nVirt}},
            {R, {0, nOBS}},
            {S, {0, nOBS}}
        };

        CCSD_F12bResult result;

        //-- MP2-F12/3C intermediates --
        if (console)
            *console << "Evaluating the MP2-F12/3C intermediates...\n";
        const F12MOIntegrals &ints = calc.ints;
        const MP2F12Intermediates intermediates = EvaluateMP2F12Intermediates(ints, f12Idxs);
        result.E_MP2 = EvaluateMP2Energy(ints.eri, calc.hfWfn.Ep_Alpha, f12Idxs);

        //Energy of the fixed amplitude geminals, 2 c V + c (B - (e_i + e_j) X) c
        for (Index i = 0; i < nOcc; i++)
            for (Index j = 0; j < nOcc; j++)
                result.E_F12 += CalculateMP2_F12_FixedEnergyElectronPair(i, j, intermediates.F, intermediates.V,
                                                                         intermediates.X, intermediates.B);

        if (thread && thread->ShouldStop())
            return result;

        //-- Coupling of the conventional doubles to the geminals --
        //V_ab^kl = <ab| r12^-1 f12 |kl> - sum_rs <ab| r12^-1 |rs><rs| f12 |kl>. As in CCSD-F12x the projector has
        //no CABS summations here, unlike the V_ij^kl intermediate of the geminal energy. It is formed indexed
        //(k, l, a, b) using <ab|O|kl> = <kl|O|ab>.
        Tensor4 V_klab = Slice<Indices<K, L, A, B> >(idxs, ints.f12Repulsion);
        //abrs,klrs->klab
        V_klab -= Einsum<Indices<A, B, R, S>, Indices<K, L, R, S>, Indices<K, L, A, B> >(
            Slice<Indices<A, B, R, S> >(idxs, ints.eriOBS),
            Slice<Indices<K, L, R, S> >(idxs, ints.f12));

        //Doubles residual: sum_kl c_ij^kl (C_kl^ab + V_ab^kl); energy: sum_kl c_ij^kl C_kl^ab
        CC::RCCSD_F12Coupling coupling;
        coupling.residual = ContractFixedAmplitudes(intermediates.C + V_klab);
        coupling.energy = ContractFixedAmplitudes(intermediates.C);

        if (thread && thread->ShouldStop())
            return result;

        //-- Closed-shell CCSD in the OBS, solved with the geminal couplings --
        const CC::RCCSD_Wfn ccsd = CC::SolveRCCSD_DIIS(ints.eriOBS, calc.hfWfn.Ep_Alpha.head(nOBS), nOcc,
                                                       &coupling, thread, console);
        result.converged = ccsd.converged;
        result.E_CCSD = ccsd.E_CCSD;
        result.E_Coupling = ccsd.E_F12_Coupling;

        if (console && result.converged)
            *console << std::setprecision(12)
                    << "E(MP2) = " << result.E_MP2 << " hartrees\n"
                    << "E(F12 geminals) = " << result.E_F12 << " hartrees\n"
                    << "E(CCSD-F12b correlation) = " << result.CorrelationEnergy() << " hartrees\n"
                    << "E(Total) = E(HF) + E(CABS singles) + E(CCSD-F12b correlation) = "
                    << calc.hfWfn.E_HF + calc.E_CABSSingles + result.CorrelationEnergy() << " hartrees\n\n";

        return result;
    }
    //
    //
}
