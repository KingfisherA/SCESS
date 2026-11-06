#include <Quantum/Q_PCH.h>
//
#include "Dahlia/Tools/Threading.h"
//
#include <cfloat>

#include "Quantum/CC/CC.h"
#include "Quantum/HF/HF.h"
#include "Quantum/Integrals/IntegralTransforms.h"
#include "Tensor/TensorEinsum.h"
#include "Tensor/TensorSlice.h"

//Closed-shell CCSD in spatial orbitals. The amplitude equations are those of Hirata et al., J. Chem. Phys. 120, 2581
//(2004), Eqs. (35)-(45), in the form of PySCF's rccsd.update_amps and rintermediates (from which the einsum strings in
//the comments are taken), specialised to canonical orbitals (f_ov = 0, diagonal f_oo and f_vv).
namespace SCESS::CC {
    using namespace Eigen;
    //
    //
    namespace {
        using Tensor2 = Tensor<double, 2>;
        using Tensor4 = Tensor<double, 4>;

        enum {
            //Occupied orbitals
            I, J, K, L,
            //Virtual orbitals
            A, B, C, D
        };

        //The integral blocks in chemists' notation, named as in PySCF, e.g. ovov(k, c, l, d) = (kc|ld)
        struct RCCSDIntegrals {
            Tensor4 ovov, oovv, ovvo, ovvv, ovoo;
            //<kl|ij> = (ki|lj), indexed (k, l, i, j)
            Tensor4 oooo_phys;
            //<ab|cd> = (ac|bd), indexed (a, b, c, d)
            Tensor4 vvvv_phys;
            //<ij|ab>, indexed (i, j, a, b)
            Tensor4 oovv_phys;
            //2 <ij|ab> - <ij|ba>
            Tensor4 L_oovv_phys;
            //2 (kc|ld) - (kd|lc)
            Tensor4 L_ovov;
        };

        //Slices the blocks out of the physicists' notation integrals <pq|rs> = (pr|qs)
        RCCSDIntegrals MakeIntegralBlocks(const TEITensor &eri, Index nOcc, Index nVirt) {
            const TensorSliceIndices idxs = {
                {I, {0, nOcc}}, {J, {0, nOcc}}, {K, {0, nOcc}}, {L, {0, nOcc}},
                {A, {nOcc, nVirt}}, {B, {nOcc, nVirt}}, {C, {nOcc, nVirt}}, {D, {nOcc, nVirt}}
            };

            RCCSDIntegrals ints;
            //(kc|ld) = <kl|cd>
            ints.ovov = Einsum<Indices<K, L, C, D>, Indices<K, C, L, D> >(Slice<Indices<K, L, C, D> >(idxs, eri));
            //(ki|ac) = <ka|ic>
            ints.oovv = Einsum<Indices<K, A, I, C>, Indices<K, I, A, C> >(Slice<Indices<K, A, I, C> >(idxs, eri));
            //(kc|ai) = <ka|ci>
            ints.ovvo = Einsum<Indices<K, A, C, I>, Indices<K, C, A, I> >(Slice<Indices<K, A, C, I> >(idxs, eri));
            //(kd|ac) = <ka|dc>
            ints.ovvv = Einsum<Indices<K, A, D, C>, Indices<K, D, A, C> >(Slice<Indices<K, A, D, C> >(idxs, eri));
            //(lc|ki) = <lk|ci>
            ints.ovoo = Einsum<Indices<L, K, C, I>, Indices<L, C, K, I> >(Slice<Indices<L, K, C, I> >(idxs, eri));
            ints.oooo_phys = Slice<Indices<K, L, I, J> >(idxs, eri);
            ints.vvvv_phys = Slice<Indices<A, B, C, D> >(idxs, eri);
            ints.oovv_phys = Slice<Indices<I, J, A, B> >(idxs, eri);

            ints.L_oovv_phys = 2.0 * ints.oovv_phys -
                               Einsum<Indices<I, J, B, A>, Indices<I, J, A, B> >(ints.oovv_phys);
            ints.L_ovov = 2.0 * ints.ovov - Einsum<Indices<K, D, L, C>, Indices<K, C, L, D> >(ints.ovov);
            return ints;
        }

        //X_ij^ab + X_ji^ba
        Tensor4 PlusTranspose(const Tensor4 &X) {
            return X + Einsum<Indices<J, I, B, A>, Indices<I, J, A, B> >(X);
        }

        //tau_ij^ab = t_ij^ab + t_i^a t_j^b
        Tensor4 MakeTau(const Tensor2 &t1, const Tensor4 &t2) {
            //ia,jb->ijab
            return t2 + Einsum<Indices<I, A>, Indices<J, B>, Indices<I, J, A, B> >(t1, t1);
        }

        double CorrelationEnergy(const RCCSDIntegrals &ints, const Tensor2 &t1, const Tensor4 &t2) {
            //f_ia = 0 so E = sum_ijab (2 <ij|ab> - <ij|ba>) tau_ij^ab
            const Tensor<double, 0> e = (ints.L_oovv_phys * MakeTau(t1, t2)).sum();
            return e(0);
        }

        double CouplingEnergy(const RCCSD_F12Coupling *f12, const Tensor4 &t2) {
            if (!f12)
                return 0.0;
            //sum_ijab (2 e_ij^ab - e_ij^ba) t_ij^ab
            const Tensor<double, 0> e = ((2.0 * f12->energy -
                                          Einsum<Indices<I, J, B, A>, Indices<I, J, A, B> >(f12->energy)) * t2).sum();
            return e(0);
        }

        //One Jacobi update of the amplitudes, t <- R(t) / D
        void UpdateAmplitudes(const RCCSDIntegrals &ints,
                              const Tensor2 &eia,
                              const Tensor4 &eijab,
                              const Tensor4 *t2Driver,
                              Tensor2 &t1,
                              Tensor4 &t2) {
            const Tensor4 tau = MakeTau(t1, t2);

            //-- Foo, Fvv & Fov without the diagonal Fock matrix, which cancels against the denominators --
            //kcld,ilcd->ki (the t1 t1 part of tau gives the kcld,ic,ld->ki terms)
            const Tensor2 Foo = Einsum<Indices<K, C, L, D>, Indices<I, L, C, D>, Indices<K, I> >(ints.L_ovov, tau);
            //-kcld,klad->ac
            const Tensor2 Fvv = -Einsum<Indices<K, C, L, D>, Indices<K, L, A, D>, Indices<A, C> >(ints.L_ovov, tau);
            //kcld,ld->kc
            const Tensor2 Fov = Einsum<Indices<K, C, L, D>, Indices<L, D>, Indices<K, C> >(ints.L_ovov, t1);

            //-- T1 equation --
            //ac,ic->ia
            Tensor2 t1new = Einsum<Indices<A, C>, Indices<I, C>, Indices<I, A> >(Fvv, t1);
            //-ki,ka->ia
            t1new -= Einsum<Indices<K, I>, Indices<K, A>, Indices<I, A> >(Foo, t1);
            //2 kc,kica->ia - kc,ikca->ia
            t1new += 2.0 * Einsum<Indices<K, C>, Indices<K, I, C, A>, Indices<I, A> >(Fov, t2);
            t1new -= Einsum<Indices<K, C>, Indices<I, K, C, A>, Indices<I, A> >(Fov, t2);
            {
                //kc,ic,ka->ia
                const Tensor2 Fov_t1 = Einsum<Indices<K, C>, Indices<I, C>, Indices<K, I> >(Fov, t1);
                t1new += Einsum<Indices<K, I>, Indices<K, A>, Indices<I, A> >(Fov_t1, t1);
            }
            //2 kcai,kc->ia - kiac,kc->ia
            t1new += 2.0 * Einsum<Indices<K, C, A, I>, Indices<K, C>, Indices<I, A> >(ints.ovvo, t1);
            t1new -= Einsum<Indices<K, I, A, C>, Indices<K, C>, Indices<I, A> >(ints.oovv, t1);
            //2 kdac,ikcd->ia - kcad,ikcd->ia (and the kdac,kd,ic / kcad,kd,ic terms through tau)
            Parallel(t1new) += 2.0 * Einsum<Indices<K, D, A, C>, Indices<I, K, C, D>, Indices<I, A> >(ints.ovvv, tau);
            Parallel(t1new) -= Einsum<Indices<K, C, A, D>, Indices<I, K, C, D>, Indices<I, A> >(ints.ovvv, tau);
            //-2 lcki,klac->ia + kcli,klac->ia (and the lcki,lc,ka / kcli,lc,ka terms through tau)
            t1new -= 2.0 * Einsum<Indices<L, C, K, I>, Indices<K, L, A, C>, Indices<I, A> >(ints.ovoo, tau);
            t1new += Einsum<Indices<K, C, L, I>, Indices<K, L, A, C>, Indices<I, A> >(ints.ovoo, tau);

            //-- T2 equation --
            Tensor4 t2new;
            {
                //kibc,ka->abic with -t1, plus (ia|cb) indexed (a, b, i, c)
                Tensor4 tmp2 = Evaluate(-Einsum<Indices<K, I, B, C>, Indices<K, A>, Indices<A, B, I, C> >(ints.oovv, t1));
                tmp2 += Einsum<Indices<I, A, C, B>, Indices<A, B, I, C> >(ints.ovvv);
                //abic,jc->ijab
                t2new = PlusTranspose(Evaluate(
                    Einsum<Indices<A, B, I, C>, Indices<J, C>, Indices<I, J, A, B> >(tmp2, t1)));
            }
            {
                //kcai,jc->akij plus (ia|jk) indexed (a, k, i, j)
                Tensor4 tmp2 = Einsum<Indices<K, C, A, I>, Indices<J, C>, Indices<A, K, I, J> >(ints.ovvo, t1);
                tmp2 += Einsum<Indices<I, A, J, K>, Indices<A, K, I, J> >(ints.ovoo);
                //akij,kb->ijab
                t2new -= PlusTranspose(Einsum<Indices<A, K, I, J>, Indices<K, B>, Indices<I, J, A, B> >(tmp2, t1));
            }
            //Driving term (ia|jb) = <ij|ab>
            t2new += ints.oovv_phys;
            if (t2Driver)
                t2new += *t2Driver;

            //Loo = Foo + 2 lcki,lc->ki - kcli,lc->ki
            Tensor2 Loo = Foo;
            Loo += 2.0 * Einsum<Indices<L, C, K, I>, Indices<L, C>, Indices<K, I> >(ints.ovoo, t1);
            Loo -= Einsum<Indices<K, C, L, I>, Indices<L, C>, Indices<K, I> >(ints.ovoo, t1);
            //Lvv = Fvv + 2 kdac,kd->ac - kcad,kd->ac
            Tensor2 Lvv = Fvv;
            Lvv += 2.0 * Einsum<Indices<K, D, A, C>, Indices<K, D>, Indices<A, C> >(ints.ovvv, t1);
            Lvv -= Einsum<Indices<K, C, A, D>, Indices<K, D>, Indices<A, C> >(ints.ovvv, t1);

            //Woooo(k, l, i, j) = lcki,jc->klij + kclj,ic->klij + kcld,ijcd->klij (with tau) + <kl|ij>
            Tensor4 Woooo = ints.oooo_phys;
            Woooo += Einsum<Indices<L, C, K, I>, Indices<J, C>, Indices<K, L, I, J> >(ints.ovoo, t1);
            Woooo += Einsum<Indices<K, C, L, J>, Indices<I, C>, Indices<K, L, I, J> >(ints.ovoo, t1);
            Parallel(Woooo) += Einsum<Indices<K, C, L, D>, Indices<I, J, C, D>, Indices<K, L, I, J> >(ints.ovov, tau);

            //Wvvvv(a, b, c, d) = -kdac,kb->abcd - kcbd,ka->abcd + <ab|cd>
            Tensor4 Wvvvv = ints.vvvv_phys;
            Parallel(Wvvvv) -= Einsum<Indices<K, D, A, C>, Indices<K, B>, Indices<A, B, C, D> >(ints.ovvv, t1);
            Parallel(Wvvvv) -= Einsum<Indices<K, C, B, D>, Indices<K, A>, Indices<A, B, C, D> >(ints.ovvv, t1);

            //Wvoov(a, k, i, c)
            Tensor4 Wvoov = Einsum<Indices<K, C, A, I>, Indices<A, K, I, C> >(ints.ovvo);
            //kcad,id->akic - kcli,la->akic
            Parallel(Wvoov) += Einsum<Indices<K, C, A, D>, Indices<I, D>, Indices<A, K, I, C> >(ints.ovvv, t1);
            Wvoov -= Einsum<Indices<K, C, L, I>, Indices<L, A>, Indices<A, K, I, C> >(ints.ovoo, t1);
            //-1/2 ldkc,ilda->akic - 1/2 lckd,ilad->akic + ldkc,ilad->akic
            Parallel(Wvoov) -= 0.5 * Einsum<Indices<L, D, K, C>, Indices<I, L, D, A>, Indices<A, K, I, C> >(
                ints.ovov, t2);
            Parallel(Wvoov) -= 0.5 * Einsum<Indices<L, C, K, D>, Indices<I, L, A, D>, Indices<A, K, I, C> >(
                ints.ovov, t2);
            Parallel(Wvoov) += Einsum<Indices<L, D, K, C>, Indices<I, L, A, D>, Indices<A, K, I, C> >(ints.ovov, t2);
            {
                //-ldkc,id,la->akic
                const Tensor4 tmp = Einsum<Indices<L, D, K, C>, Indices<I, D>, Indices<L, K, C, I> >(ints.ovov, t1);
                Wvoov -= Einsum<Indices<L, K, C, I>, Indices<L, A>, Indices<A, K, I, C> >(tmp, t1);
            }

            //Wvovo(a, k, c, i)
            Tensor4 Wvovo = Einsum<Indices<K, I, A, C>, Indices<A, K, C, I> >(ints.oovv);
            //kdac,id->akci - lcki,la->akci
            Parallel(Wvovo) += Einsum<Indices<K, D, A, C>, Indices<I, D>, Indices<A, K, C, I> >(ints.ovvv, t1);
            Wvovo -= Einsum<Indices<L, C, K, I>, Indices<L, A>, Indices<A, K, C, I> >(ints.ovoo, t1);
            //-1/2 lckd,ilda->akci
            Parallel(Wvovo) -= 0.5 * Einsum<Indices<L, C, K, D>, Indices<I, L, D, A>, Indices<A, K, C, I> >(
                ints.ovov, t2);
            {
                //-lckd,id,la->akci
                const Tensor4 tmp = Einsum<Indices<L, C, K, D>, Indices<I, D>, Indices<L, C, K, I> >(ints.ovov, t1);
                Wvovo -= Einsum<Indices<L, C, K, I>, Indices<L, A>, Indices<A, K, C, I> >(tmp, t1);
            }

            //klij,klab->ijab
            Parallel(t2new) += Einsum<Indices<K, L, I, J>, Indices<K, L, A, B>, Indices<I, J, A, B> >(Woooo, tau);
            //abcd,ijcd->ijab: the particle-particle ladder, the most expensive term
            Parallel(t2new) += Einsum<Indices<A, B, C, D>, Indices<I, J, C, D>, Indices<I, J, A, B> >(Wvvvv, tau);
            //ac,ijcb->ijab
            t2new += PlusTranspose(Evaluate(Einsum<Indices<A, C>, Indices<I, J, C, B>, Indices<I, J, A, B> >(Lvv, t2)));
            //ki,kjab->ijab
            t2new -= PlusTranspose(Einsum<Indices<K, I>, Indices<K, J, A, B>, Indices<I, J, A, B> >(Loo, t2));
            {
                //2 akic,kjcb->ijab - akci,kjcb->ijab
                Tensor4 tmp = Evaluate(2.0 * Einsum<Indices<A, K, I, C>, Indices<K, J, C, B>, Indices<I, J, A, B> >(
                    Wvoov, t2));
                Parallel(tmp) -= Einsum<Indices<A, K, C, I>, Indices<K, J, C, B>, Indices<I, J, A, B> >(Wvovo, t2);
                t2new += PlusTranspose(tmp);
                //akic,kjbc->ijab
                tmp = Evaluate(Einsum<Indices<A, K, I, C>, Indices<K, J, B, C>, Indices<I, J, A, B> >(Wvoov, t2));
                t2new -= PlusTranspose(tmp);
                //bkci,kjac->ijab
                tmp = Evaluate(Einsum<Indices<B, K, C, I>, Indices<K, J, A, C>, Indices<I, J, A, B> >(Wvovo, t2));
                t2new -= PlusTranspose(tmp);
            }

            t1 = t1new / eia;
            t2 = t2new / eijab;
        }

        //Direct inversion in the iterative subspace over a flattened amplitude vector
        class AmplitudeDIIS {
        public:
            explicit AmplitudeDIIS(size_t maxVectors) : mMaxVectors(maxVectors) {
            }

            //Stores the amplitudes and their change in the last iteration, returning the extrapolated amplitudes
            VectorXd Extrapolate(const VectorXd &amplitudes, const VectorXd &error) {
                mAmplitudes.push_back(amplitudes);
                mErrors.push_back(error);
                if (mErrors.size() > mMaxVectors) {
                    mAmplitudes.erase(mAmplitudes.begin());
                    mErrors.erase(mErrors.begin());
                }

                const Index n = (Index) mErrors.size();
                MatrixXd B = MatrixXd::Constant(n + 1, n + 1, -1.0);
                B(n, n) = 0.0;
                for (Index i = 0; i < n; i++)
                    for (Index j = 0; j < n; j++)
                        B(i, j) = mErrors[i].dot(mErrors[j]);

                //Scale for the conditioning of the solve
                const double maxVal = B.topLeftCorner(n, n).cwiseAbs().maxCoeff();
                if (maxVal > 0.0)
                    B.topLeftCorner(n, n) /= maxVal;

                VectorXd rhs = VectorXd::Zero(n + 1);
                rhs(n) = -1.0;
                const VectorXd coefficients = B.fullPivLu().solve(rhs);

                VectorXd result = VectorXd::Zero(amplitudes.size());
                for (Index i = 0; i < n; i++)
                    result += coefficients(i) * mAmplitudes[i];
                return result;
            }

        private:
            size_t mMaxVectors;
            std::vector<VectorXd> mAmplitudes, mErrors;
        };

        VectorXd Flatten(const Tensor2 &t1, const Tensor4 &t2) {
            VectorXd v(t1.size() + t2.size());
            v.head(t1.size()) = Map<const VectorXd>(t1.data(), t1.size());
            v.tail(t2.size()) = Map<const VectorXd>(t2.data(), t2.size());
            return v;
        }

        void Unflatten(const VectorXd &v, Tensor2 &t1, Tensor4 &t2) {
            Map<VectorXd>(t1.data(), t1.size()) = v.head(t1.size());
            Map<VectorXd>(t2.data(), t2.size()) = v.tail(t2.size());
        }
    }
    //
    //
    RCCSD_Wfn SolveRCCSD_DIIS(const TEITensor &eri,
                              const VectorXd &orbitalEnergies,
                              const Index nOcc,
                              const RCCSD_F12Coupling *f12,
                              const DH::Thread *thread,
                              std::stringstream *console) {
        //Convergence thresholds on the change in energy and the root mean square change in the amplitudes
        constexpr double ENERGY_CONVERGENCE = 1.0E-10;
        constexpr double AMPLITUDE_CONVERGENCE = 1.0E-8;
        constexpr size_t MAX_ITERATIONS = 100;
        //DIIS starts once a few Jacobi iterations have built up a subspace
        constexpr size_t DIIS_START = 2;
        constexpr size_t MAX_DIIS_VECS = 8;

        const Index nOrbs = eri.dimension(0);
        const Index nVirt = nOrbs - nOcc;

        if (console)
            *console << std::setprecision(12)
                    << (f12
                            ? "--- Closed-shell CCSD with F12 geminals (CCSD-F12b) with DIIS acceleration ---\n"
                            : "--- Closed-shell Coupled Cluster Singles & Doubles (CCSD) with DIIS acceleration ---\n");

        const RCCSDIntegrals ints = MakeIntegralBlocks(eri, nOcc, nVirt);
        const Tensor4 *t2Driver = f12 ? &f12->residual : nullptr;

        //Denominators e_i - e_a and e_i + e_j - e_a - e_b
        Tensor2 eia(nOcc, nVirt);
        for (Index i = 0; i < nOcc; i++)
            for (Index a = 0; a < nVirt; a++)
                eia(i, a) = orbitalEnergies(i) - orbitalEnergies(nOcc + a);
        Tensor4 eijab(nOcc, nOcc, nVirt, nVirt);
        for (Index i = 0; i < nOcc; i++)
            for (Index j = 0; j < nOcc; j++)
                for (Index a = 0; a < nVirt; a++)
                    for (Index b = 0; b < nVirt; b++)
                        eijab(i, j, a, b) = eia(i, a) + eia(j, b);

        //MP2 (plus geminal coupling) initial guess
        Tensor2 t1(nOcc, nVirt);
        t1.setZero();
        Tensor4 t2 = ints.oovv_phys;
        if (t2Driver)
            t2 += *t2Driver;
        t2 = t2 / eijab;

        AmplitudeDIIS diis(MAX_DIIS_VECS);
        double E_last = DBL_MAX;
        RCCSD_Wfn wfn;

        for (size_t iteration = 0; iteration < MAX_ITERATIONS; iteration++) {
            if (thread && thread->ShouldStop())
                return {};

            const VectorXd previous = Flatten(t1, t2);
            UpdateAmplitudes(ints, eia, eijab, t2Driver, t1, t2);
            VectorXd updated = Flatten(t1, t2);
            const VectorXd change = updated - previous;

            if (iteration >= DIIS_START) {
                if (iteration == DIIS_START && console)
                    *console << "Turning on DIIS\n";
                Unflatten(diis.Extrapolate(updated, change), t1, t2);
            } else
                diis.Extrapolate(updated, change);

            const double E = CorrelationEnergy(ints, t1, t2) + CouplingEnergy(f12, t2);
            const double dE = std::abs(E - E_last);
            const double rmsChange = change.norm() / std::sqrt((double) change.size());
            E_last = E;

            if (console)
                *console << "   CCSD(" << iteration << "): E(corr) = " << E << "      dE = " << dE
                        << "      rms(dt) = " << rmsChange << '\n';

            if (dE < ENERGY_CONVERGENCE && rmsChange < AMPLITUDE_CONVERGENCE) {
                wfn.converged = true;
                break;
            }
        }

        if (!wfn.converged) {
            if (console)
                *console << "!!! CCSD procedure failed to converge !!!\n\n\n";
            return wfn;
        }

        wfn.E_CCSD = E_last;
        wfn.E_F12_Coupling = CouplingEnergy(f12, t2);
        wfn.t1 = std::move(t1);
        wfn.t2 = std::move(t2);

        if (console) {
            *console << "---- CCSD procedure converged! ----\n"
                    << "E(CCSD correlation) = " << wfn.E_CCSD << " hartrees\n";
            if (f12)
                *console << "   of which E(T2-geminal coupling) = " << wfn.E_F12_Coupling << " hartrees\n";
        }
        return wfn;
    }
    //
    RCCSD_Wfn RCCSD_DIIS(const HFWfn &rhfWfn, const Integrals::MolInts &mints,
                         const DH::Thread *thread, std::stringstream *console) {
        const Index nOcc = std::distance(rhfWfn.Occ_Alpha.begin(), std::ranges::find(rhfWfn.Occ_Alpha, 0));
        const TEITensor eri = Integrals::Transforms::TransformTEI(mints.eri, rhfWfn.C_Alpha, true);

        RCCSD_Wfn wfn = SolveRCCSD_DIIS(eri, rhfWfn.Ep_Alpha, nOcc, nullptr, thread, console);
        if (console && wfn.converged)
            *console << "E(Total) = E(Hartree Fock) + E(CCSD) = " << rhfWfn.E_HF + wfn.E_CCSD << " hartrees\n\n\n";
        return wfn;
    }
    //
    //
}
