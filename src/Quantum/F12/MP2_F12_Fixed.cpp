#include<Quantum/Q_PCH.h>
#include "F12.h"
//
#include "Tensor/TensorOperations.h"
#include "Tensor/TensorSlice.h"
#include "Tensor/TensorEinsum.h"
#include "Quantum/Integrals/Integrals.h"
#include "Quantum/Integrals/IntegralsF12.h"
#include "Quantum/Integrals/IntegralTransforms.h"
#include "Quantum/HF/HF.h"

//
//
namespace SCESS::F12 {
    using namespace Eigen;
    //-- Indices --
    enum {
        //Complete virtual space orbitals
        alpha, beta,
        //Occupied orbitals
        I, J, K, L, M, N, O,
        //OBS Virtual orbitals
        A, B, C, D,
        //All OBS orbitals
        P, Q, R, S, T, U,
        //CABS virtual orbitals
        X, Y, Z, W
    };
    //
    //
    Tensor4 EvaluateW(const Tensor4 &f12, F12Idxs f12Idxs) {
        TensorSliceIndices idxs = {
            //All Occupied orbital indices
            {I, {0, f12Idxs.nOcc}},
            {J, {0, f12Idxs.nOcc}},
            //
            {alpha, {f12Idxs.nOcc, f12Idxs.nVirt + f12Idxs.nCABS}},
            {beta, {f12Idxs.nOcc, f12Idxs.nVirt + f12Idxs.nCABS}}
        };

        Tensor4 W = (3.0 / 8.0) * Slice<Indices<alpha, beta, I, J> >(idxs, f12);
        W -= (1.0 / 8.0) * Slice<Indices<alpha, beta, I, J> >(idxs, f12);

        return W;
    }
    //
    void EvaluateAdditionalMP2Intermediates(const Integrals::MolInts &mints,
                                            const Integrals::F12::F12Integrals &f12Ints,
                                            F12Idxs f12Idxs,
                                            TEITensor *V_pqxy_out,
                                            TEITensor *h_apxy_out
    ) {
        TensorSliceIndices idxs = {
            //All Occupied orbital indices
            {I, {0, f12Idxs.nOcc}},
            {J, {0, f12Idxs.nOcc}},
            {M, {0, f12Idxs.nOcc}},
            //All virtual orbital indices
            {A, {f12Idxs.nOcc, f12Idxs.nVirt}},
            {B, {f12Idxs.nOcc, f12Idxs.nVirt}},
            {C, {f12Idxs.nOcc, f12Idxs.nVirt}},
            {D, {f12Idxs.nOcc, f12Idxs.nVirt}},
            //All OBS indices
            {P, {0, f12Idxs.nOBS}},
            {Q, {0, f12Idxs.nOBS}},
            {R, {0, f12Idxs.nOBS}},
            {S, {0, f12Idxs.nOBS}},
            {T, {0, f12Idxs.nOBS}},
            {U, {0, f12Idxs.nOBS}},
            //CABS indices
            {X, {f12Idxs.nOBS, f12Idxs.nCABS}},
        };

        const TEITensor &eri = mints.eri;
        const TEITensor &f12 = f12Ints.f12;
        const TEITensor &f12Sq = f12Ints.f12Sq;
        const TEITensor &f12DelSq = f12Ints.f12DelSq;
        const TEITensor &f12Repulsion = f12Ints.f12Repulsion;

        //Make the V_pqxy intermediate
        //Note that i=x and j=y as fixed amplitude approximation employed
        Tensor4 V_pqij = Slice<Indices<P, Q, I, J> >(idxs, f12Repulsion);

        //abmx,cdmx->abcd
        V_pqij -= Einsum<Indices<P, Q, M, X>, Indices<M, X, I, J>, Indices<P, Q, I, J> >(
            Slice<Indices<P, Q, M, X> >(idxs, eri),
            Slice<Indices<M, X, I, J> >(idxs, f12));

        //abxm,cdxm->abcd
        V_pqij -= Einsum<Indices<P, Q, X, M>, Indices<X, M, I, J>, Indices<P, Q, I, J> >(
            Slice<Indices<P, Q, X, M> >(idxs, eri),
            Slice<Indices<X, M, I, J> >(idxs, f12));


        V_pqij -= Einsum<Indices<P, Q, R, S>, Indices<R, S, I, J>, Indices<P, Q, I, J> >(
            Slice<Indices<P, Q, R, S> >(idxs, eri),
            Slice<Indices<R, S, I, J> >(idxs, f12));

        Tensor2 h_xp = Slice<Indices<X, A> >(idxs, MatrixToTensor(mints.T + mints.V));

        Tensor4 h_apxy =
                Einsum<Indices<A, X, I, J>, Indices<X, P>, Indices<A, P, I, J> >(
                    Slice<Indices<A, X, I, J> >(idxs, f12),
                    h_xp
                );

        h_apxy += Einsum<Indices<X, A, I, J>, Indices<X, P>, Indices<A, P, I, J> >(
            Slice<Indices<X, A, I, J> >(idxs, f12),
            h_xp
        );


        if (V_pqxy_out)
            *V_pqxy_out = V_pqij;

        if (h_apxy_out)
            *h_apxy_out = h_apxy;
    }
    //
    //
    MP2F12Result MP2_F12_Fixed(const HF_Wfn &refWfn,
                               const Integrals::MolInts &mints,
                               const Integrals::F12::F12Integrals &f12Ints,
                               F12Idxs f12Idxs,
                               Tensor2 *F_pq_out,
                               TEITensor *V_ijxy_out,
                               TEITensor *C_xyab_out,
                               TEITensor *X_xywz_out,
                               TEITensor *B_xywz_out) {
        // * We assume that the reference wavefunction is actually a closed shell RHF hence alternating
        //      alpha and beta spins (fix this in the future to support open shell :( )
        // * We need | Psi > = |ijk...> |alpha beta alpha beta ... > (alpha and beta are stored alternating)
        // * All calculations are done in SPIN ORBITALS so we need to integrate spin out the spin!
        // * All EVEN indices are spin alpha and ODD are spin beta

        TensorSliceIndices idxs = {
            //All Occupied orbital indices
            {I, {0, f12Idxs.nOcc}},
            {J, {0, f12Idxs.nOcc}},
            {K, {0, f12Idxs.nOcc}},
            {L, {0, f12Idxs.nOcc}},
            {M, {0, f12Idxs.nOcc}},
            {N, {0, f12Idxs.nOcc}},
            {O, {0, f12Idxs.nOcc}},
            //All virtual orbital indices
            {A, {f12Idxs.nOcc, f12Idxs.nVirt}},
            {B, {f12Idxs.nOcc, f12Idxs.nVirt}},
            {C, {f12Idxs.nOcc, f12Idxs.nVirt}},
            {D, {f12Idxs.nOcc, f12Idxs.nVirt}},
            //All OBS indices
            {P, {0, f12Idxs.nOBS}},
            {Q, {0, f12Idxs.nOBS}},
            {R, {0, f12Idxs.nOBS}},
            {S, {0, f12Idxs.nOBS}},
            //CABS indices
            {X, {f12Idxs.nOBS, f12Idxs.nCABS}},
            {Y, {f12Idxs.nOBS, f12Idxs.nCABS}},
            {Z, {f12Idxs.nOBS, f12Idxs.nCABS}},
            {W, {f12Idxs.nOBS, f12Idxs.nCABS}}
        };


        //Extract out all the integral tensors we will need for this calculation
        //These integrals have been canonicalised and transformed into the CABS/OBS MO basis set
        const MatrixXd &TInts = mints.T;
        const MatrixXd &VInts = mints.V;
        const TEITensor &eri = mints.eri;
        const TEITensor &f12 = f12Ints.f12;
        const TEITensor &f12Sq = f12Ints.f12Sq;
        const TEITensor &f12DelSq = f12Ints.f12DelSq;
        const TEITensor &f12Repulsion = f12Ints.f12Repulsion;

        //Make the V intermediate
        //iiii
        TEITensor V_ = Slice<Indices<I, I, I, I> >(idxs, f12Repulsion);

        //ijrs,klrs->ijkl
        V_ -= Einsum<Indices<I, J, R, S>, Indices<K, L, R, S>, Indices<I, J, K, L> >(
            Slice<Indices<I, J, R, S> >(idxs, eri),
            Slice<Indices<K, L, R, S> >(idxs, f12));

        //ijxm,klxm->ijkl
        V_ -= Einsum<Indices<I, J, X, M>, Indices<K, L, X, M>, Indices<I, J, K, L> >(
            Slice<Indices<I, J, X, M> >(idxs, eri),
            Slice<Indices<K, L, X, M> >(idxs, f12));

        //ijmx,klmx->ijkl
        V_ -= Einsum<Indices<I, J, M, X>, Indices<K, L, M, X>, Indices<I, J, K, L> >(
            Slice<Indices<I, J, M, X> >(idxs, eri),
            Slice<Indices<K, L, M, X> >(idxs, f12));


        TEITensor X_ = Slice<Indices<I, I, I, I> >(idxs, f12Sq);
        //ijrs,klrs->ijkl
        X_ -= Einsum<Indices<I, J, R, S>, Indices<K, L, R, S>, Indices<I, J, K, L> >(
            Slice<Indices<I, J, R, S> >(idxs, f12),
            Slice<Indices<K, L, R, S> >(idxs, f12));

        //ijxm,klxm->ijkl
        X_ -= Einsum<Indices<I, J, X, M>, Indices<K, L, X, M>, Indices<I, J, K, L> >(
            Slice<Indices<I, J, X, M> >(idxs, f12),
            Slice<Indices<K, L, X, M> >(idxs, f12));

        //ijmx,klmx->ijkl
        X_ -= Einsum<Indices<I, J, M, X>, Indices<K, L, M, X>, Indices<I, J, K, L> >(
            Slice<Indices<I, J, M, X> >(idxs, f12),
            Slice<Indices<K, L, M, X> >(idxs, f12));

        //Full Fock Matrix
        Tensor<double, 2> F_ = MatrixToTensor(TInts + VInts);
        //Coloumb Fock Terms
        //piqi->pq
        Slice<Indices<P, Q> >(idxs, F_) += 2.0 * Einsum<Indices<P, I, Q, I>, Indices<P, Q> >(
            Slice<Indices<P, I, Q, I> >(idxs, eri));
        //ipix->px
        Slice<Indices<P, X> >(idxs, F_) += 2.0 * Einsum<Indices<I, P, I, X>, Indices<P, X> >(
            Slice<Indices<I, P, I, X> >(idxs, eri));
        //ixiy->xy
        Slice<Indices<X, Y> >(idxs, F_) += 2.0 * Einsum<Indices<I, X, I, Y>, Indices<X, Y> >(
            Slice<Indices<I, X, I, Y> >(idxs, eri));

        Tensor<double, 2> K_(F_.dimensions());
        K_.setZero();
        //Exchange Fock Terms
        //iipq->pq
        Slice<Indices<P, Q> >(idxs, K_) = Einsum<Indices<I, I, P, Q>, Indices<P, Q> >(
            Slice<Indices<I, I, P, Q> >(idxs, eri));
        //iipx->px
        Slice<Indices<P, X> >(idxs, K_) = Einsum<Indices<I, I, P, X>, Indices<P, X> >(
            Slice<Indices<I, I, P, X> >(idxs, eri));
        //iixy->xy
        Slice<Indices<X, Y> >(idxs, K_) = Einsum<Indices<I, I, X, Y>, Indices<X, Y> >(
            Slice<Indices<I, I, X, Y> >(idxs, eri));
        F_ -= K_;

        //Fill in the remaining elements by symmetry
        Slice<Indices<X, P> >(idxs, F_) = Einsum<Indices<P, X>, Indices<X, P> >(Slice<Indices<P, X> >(idxs, F_));
        Slice<Indices<X, P> >(idxs, K_) = Einsum<Indices<P, X>, Indices<X, P> >(Slice<Indices<P, X> >(idxs, K_));

        //Fock matrix without the exchange terms
        Tensor<double, 2> F_K = F_ + K_;

        //Makes the B and C intermediate tensors
        //, &idxs, &f12Idxs.nOcc, &f12Idxs.nOBS, &f12Idxs.nCABS]
        auto MakeB_C = [&f12Sq, &f12, &F_, &F_K, &K_, &idxs]
        (TEITensor &B_, TEITensor &C_)-> void {
            //-- F12 Contributions --
            //pr,lkqr-> plkq -> klpq
            TEITensor tmp = Einsum<Indices<P, R>, Indices<L, K, Q, R>, Indices<K, L, P, Q> >(
                Slice<Indices<P, R> >(idxs, F_K),
                Slice<Indices<L, K, Q, R> >(idxs, f12));
            //mnpq,klpq-> mnkl -> klmn
            B_ -= Einsum<Indices<M, N, P, Q>, Indices<K, L, P, Q>, Indices<K, L, M, N> >(
                Slice<Indices<M, N, P, Q> >(idxs, f12),
                tmp);

            //qr,klpr-> qklp -> klpq
            tmp = Einsum<Indices<Q, R>, Indices<K, L, P, R>, Indices<K, L, P, Q> >(
                Slice<Indices<Q, R> >(idxs, F_K),
                Slice<Indices<K, L, P, R> >(idxs, f12));
            //mnpq,klpq-> mnkl -> klmn
            B_ -= Einsum<Indices<M, N, P, Q>, Indices<K, L, P, Q>, Indices<K, L, M, N> >(
                Slice<Indices<M, N, P, Q> >(idxs, f12),
                tmp);

            //px,lkqx-> plkq -> klpq
            tmp = Einsum<Indices<P, X>, Indices<L, K, Q, X>, Indices<K, L, P, Q> >(
                Slice<Indices<P, X> >(idxs, F_K),
                Slice<Indices<L, K, Q, X> >(idxs, f12));
            //mnpq,klpq-> mnkl -> klmn
            B_ -= Einsum<Indices<M, N, P, Q>, Indices<K, L, P, Q>, Indices<K, L, M, N> >(
                Slice<Indices<M, N, P, Q> >(idxs, f12),
                tmp);

            //qx,klpx-> qklp -> klpq
            tmp = Einsum<Indices<Q, X>, Indices<K, L, P, X>, Indices<K, L, P, Q> >(
                Slice<Indices<Q, X> >(idxs, F_K),
                Slice<Indices<K, L, P, X> >(idxs, f12));
            //mnpq,klpq-> mnkl -> klmn
            B_ -= Einsum<Indices<M, N, P, Q>, Indices<K, L, P, Q>, Indices<K, L, M, N> >(
                Slice<Indices<M, N, P, Q> >(idxs, f12),
                tmp);

            //yp,klip-> ykli ->  kliy
            tmp = Einsum<Indices<Y, P>, Indices<K, L, I, P>, Indices<K, L, I, Y> >(
                Slice<Indices<Y, P> >(idxs, F_K),
                Slice<Indices<K, L, I, P> >(idxs, f12));
            //mniy,kliy-> mnkl -> klmn
            B_ -= Einsum<Indices<M, N, I, Y>, Indices<K, L, I, Y>, Indices<K, L, M, N> >(
                Slice<Indices<M, N, I, Y> >(idxs, f12),
                tmp);

            //xp,lkjp-> xlkj -> klxj
            tmp = Einsum<Indices<X, P>, Indices<L, K, J, P>, Indices<K, L, X, J> >(
                Slice<Indices<X, P> >(idxs, F_K),
                Slice<Indices<L, K, J, P> >(idxs, f12));
            //nmjx,klxj-> nmkl ->  klmn
            B_ -= Einsum<Indices<N, M, J, X>, Indices<K, L, X, J>, Indices<K, L, M, N> >(
                Slice<Indices<N, M, J, X> >(idxs, f12),
                tmp);

            //ip,klpy-> ikly -> kliy
            tmp = Einsum<Indices<I, P>, Indices<K, L, P, Y>, Indices<K, L, I, Y> >(
                Slice<Indices<I, P> >(idxs, F_K),
                Slice<Indices<K, L, P, Y> >(idxs, f12));
            //mniy,kliy-> mnkl -> klmn
            B_ -= Einsum<Indices<M, N, I, Y>, Indices<K, L, I, Y>, Indices<K, L, M, N> >(
                Slice<Indices<M, N, I, Y> >(idxs, f12),
                tmp);

            //yx,klix->kliy
            tmp = Einsum<Indices<Y, X>, Indices<K, L, I, X>, Indices<K, L, I, Y> >(
                Slice<Indices<Y, X> >(idxs, F_K),
                Slice<Indices<K, L, I, X> >(idxs, f12));
            //mniy,kliy->klmn
            B_ -= Einsum<Indices<M, N, I, Y>, Indices<K, L, I, Y>, Indices<K, L, M, N> >(
                Slice<Indices<M, N, I, Y> >(idxs, f12),
                tmp);

            //xy,lkjy->klxj
            tmp = Einsum<Indices<X, Y>, Indices<L, K, J, Y>, Indices<K, L, X, J> >(
                Slice<Indices<X, Y> >(idxs, F_K),
                Slice<Indices<L, K, J, Y> >(idxs, f12));
            //nmjx,klxj->klmn
            B_ -= Einsum<Indices<N, M, J, X>, Indices<K, L, X, J>, Indices<K, L, M, N> >(
                Slice<Indices<N, M, J, X> >(idxs, f12),
                tmp);

            //jp,lkpx->klxj
            tmp = Einsum<Indices<J, P>, Indices<L, K, P, X>, Indices<K, L, X, J> >(
                Slice<Indices<J, P> >(idxs, F_K),
                Slice<Indices<L, K, P, X> >(idxs, f12));
            //nmjx,klxj->klmn
            B_ -= Einsum<Indices<N, M, J, X>, Indices<K, L, X, J>, Indices<K, L, M, N> >(
                Slice<Indices<N, M, J, X> >(idxs, f12),
                tmp);

            //ix,klxy->kliy
            tmp = Einsum<Indices<I, X>, Indices<K, L, X, Y>, Indices<K, L, I, Y> >(
                Slice<Indices<I, X> >(idxs, F_K),
                Slice<Indices<K, L, X, Y> >(idxs, f12));
            //mniy,kliy->klmn
            B_ -= Einsum<Indices<M, N, I, Y>, Indices<K, L, I, Y>, Indices<K, L, M, N> >(
                Slice<Indices<M, N, I, Y> >(idxs, f12),
                tmp);

            //jy,lkyx->klxj
            tmp = Einsum<Indices<J, Y>, Indices<L, K, Y, X>, Indices<K, L, X, J> >(
                Slice<Indices<J, Y> >(idxs, F_K),
                Slice<Indices<L, K, Y, X> >(idxs, f12));
            //nmjx,klxj->klmn
            B_ -= Einsum<Indices<N, M, J, X>, Indices<K, L, X, J>, Indices<K, L, M, N> >(
                Slice<Indices<N, M, J, X> >(idxs, f12),
                tmp);

            //-- Y Contributions --
            //xp,lkbp->klxb
            tmp = Einsum<Indices<X, P>, Indices<L, K, B, P>, Indices<K, L, X, B> >(
                Slice<Indices<X, P> >(idxs, K_),
                Slice<Indices<L, K, B, P> >(idxs, f12));
            //klxb,nmbx->klmn
            B_ -= Einsum<Indices<K, L, X, B>, Indices<N, M, B, X>, Indices<K, L, M, N> >(
                tmp,
                Slice<Indices<N, M, B, X> >(idxs, f12));

            //yp,klap->klay
            tmp = Einsum<Indices<Y, P>, Indices<K, L, A, P>, Indices<K, L, A, Y> >(
                Slice<Indices<Y, P> >(idxs, K_),
                Slice<Indices<K, L, A, P> >(idxs, f12));
            //klay,mnay->klmn
            B_ -= Einsum<Indices<K, L, A, Y>, Indices<M, N, A, Y>, Indices<K, L, M, N> >(
                tmp,
                Slice<Indices<M, N, A, Y> >(idxs, f12));

            //xy,lkby->klxb
            tmp = Einsum<Indices<X, Y>, Indices<L, K, B, Y>, Indices<K, L, X, B> >(
                Slice<Indices<X, Y> >(idxs, K_),
                Slice<Indices<L, K, B, Y> >(idxs, f12));
            //klxb,nmbx->klmn
            B_ -= Einsum<Indices<K, L, X, B>, Indices<N, M, B, X>, Indices<K, L, M, N> >(
                tmp,
                Slice<Indices<N, M, B, X> >(idxs, f12));

            //bp,lkpx->klxb
            tmp = Einsum<Indices<B, P>, Indices<L, K, P, X>, Indices<K, L, X, B> >(
                Slice<Indices<B, P> >(idxs, K_),
                Slice<Indices<L, K, P, X> >(idxs, f12));
            //klxb,nmbx->klmn
            B_ -= Einsum<Indices<K, L, X, B>, Indices<N, M, B, X>, Indices<K, L, M, N> >(
                tmp,
                Slice<Indices<N, M, B, X> >(idxs, f12));

            //ap,klpy->klay
            tmp = Einsum<Indices<A, P>, Indices<K, L, P, Y>, Indices<K, L, A, Y> >(
                Slice<Indices<A, P> >(idxs, K_),
                Slice<Indices<K, L, P, Y> >(idxs, f12));
            //klay,mnay->klmn
            B_ -= Einsum<Indices<K, L, A, Y>, Indices<M, N, A, Y>, Indices<K, L, M, N> >(
                tmp,
                Slice<Indices<M, N, A, Y> >(idxs, f12));

            //yx,klax->klay
            tmp = Einsum<Indices<Y, X>, Indices<K, L, A, X>, Indices<K, L, A, Y> >(
                Slice<Indices<Y, X> >(idxs, K_),
                Slice<Indices<K, L, A, X> >(idxs, f12));
            //klay,mnay->klmn
            B_ -= Einsum<Indices<K, L, A, Y>, Indices<M, N, A, Y>, Indices<K, L, M, N> >(
                tmp,
                Slice<Indices<M, N, A, Y> >(idxs, f12));


            //by,lkyx->klxb
            tmp = Einsum<Indices<B, Y>, Indices<L, K, Y, X>, Indices<K, L, X, B> >(
                Slice<Indices<B, Y> >(idxs, K_),
                Slice<Indices<L, K, Y, X> >(idxs, f12));
            //klxb,nmbx->klmn
            B_ -= Einsum<Indices<K, L, X, B>, Indices<N, M, B, X>, Indices<K, L, M, N> >(
                tmp,
                Slice<Indices<N, M, B, X> >(idxs, f12));

            //ax,klxy->klay
            tmp = Einsum<Indices<A, X>, Indices<K, L, X, Y>, Indices<K, L, A, Y> >(
                Slice<Indices<A, X> >(idxs, K_),
                Slice<Indices<K, L, X, Y> >(idxs, f12));
            //klay,mnay->klmn
            B_ -= Einsum<Indices<K, L, A, Y>, Indices<M, N, A, Y>, Indices<K, L, M, N> >(
                tmp,
                Slice<Indices<M, N, A, Y> >(idxs, f12));

            //xp,klpy->klxy
            tmp = Einsum<Indices<X, P>, Indices<K, L, P, Y>, Indices<K, L, X, Y> >(
                Slice<Indices<X, P> >(idxs, K_),
                Slice<Indices<K, L, P, Y> >(idxs, f12));
            //klxy,mnxy->klmn
            B_ -= Einsum<Indices<K, L, X, Y>, Indices<M, N, X, Y>, Indices<K, L, M, N> >(
                tmp,
                Slice<Indices<M, N, X, Y> >(idxs, f12));

            //yp,lkpx->klxy
            tmp = Einsum<Indices<Y, P>, Indices<L, K, P, X>, Indices<K, L, X, Y> >(
                Slice<Indices<Y, P> >(idxs, K_),
                Slice<Indices<L, K, P, X> >(idxs, f12));
            //klxy,mnxy->klmn
            B_ -= Einsum<Indices<K, L, X, Y>, Indices<M, N, X, Y>, Indices<K, L, M, N> >(
                tmp,
                Slice<Indices<M, N, X, Y> >(idxs, f12));

            //yz,lkzx->klxy
            tmp = Einsum<Indices<Y, Z>, Indices<L, K, Z, X>, Indices<K, L, X, Y> >(
                Slice<Indices<Y, Z> >(idxs, K_),
                Slice<Indices<L, K, Z, X> >(idxs, f12));
            //klxy,mnxy->klmn
            B_ -= Einsum<Indices<K, L, X, Y>, Indices<M, N, X, Y>, Indices<K, L, M, N> >(
                tmp,
                Slice<Indices<M, N, X, Y> >(idxs, f12));

            //xz,klzy->klxy
            tmp = Einsum<Indices<X, Z>, Indices<K, L, Z, Y>, Indices<K, L, X, Y> >(
                Slice<Indices<X, Z> >(idxs, K_),
                Slice<Indices<K, L, Z, Y> >(idxs, f12));
            //klxy,mnxy->klmn
            B_ -= Einsum<Indices<K, L, X, Y>, Indices<M, N, X, Y>, Indices<K, L, M, N> >(
                tmp,
                Slice<Indices<M, N, X, Y> >(idxs, f12));


            //--F Squared Contributions--
            //kp,nmlp->klmn
            B_ += Einsum<Indices<K, P>, Indices<N, M, L, P>, Indices<K, L, M, N> >(
                Slice<Indices<K, P> >(idxs, F_K),
                Slice<Indices<N, M, L, P> >(idxs, f12Sq));
            //lp,mnkp->klmn
            B_ += Einsum<Indices<L, P>, Indices<M, N, K, P>, Indices<K, L, M, N> >(
                Slice<Indices<L, P> >(idxs, F_K),
                Slice<Indices<M, N, K, P> >(idxs, f12Sq));
            //kx,nmlx->klmn
            B_ += Einsum<Indices<K, X>, Indices<N, M, L, X>, Indices<K, L, M, N> >(
                Slice<Indices<K, X> >(idxs, F_K),
                Slice<Indices<N, M, L, X> >(idxs, f12Sq));
            //lx,mnkx->klmn
            B_ += Einsum<Indices<L, X>, Indices<M, N, K, X>, Indices<K, L, M, N> >(
                Slice<Indices<L, X> >(idxs, F_K),
                Slice<Indices<M, N, K, X> >(idxs, f12Sq));

            //-- C Tensor --
            //bx,klax->klab
            C_ += Einsum<Indices<B, X>, Indices<K, L, A, X>, Indices<K, L, A, B> >(
                Slice<Indices<B, X> >(idxs, F_),
                Slice<Indices<K, L, A, X> >(idxs, f12));
            //ax,lkbx->klab
            C_ += Einsum<Indices<A, X>, Indices<L, K, B, X>, Indices<K, L, A, B> >(
                Slice<Indices<A, X> >(idxs, F_),
                Slice<Indices<L, K, B, X> >(idxs, f12));


            //-- -FC contribution to B--
            //klab,mnab->klmn
            B_ -= Einsum<Indices<K, L, A, B>, Indices<M, N, A, B>, Indices<K, L, M, N> >(
                Slice<Indices<K, L, A, B> >(idxs, f12),
                C_);

            //--Symmetrise B--
            B_ = 0.5 * (B_ + Einsum<Indices<K, L, M, N>, Indices<M, N, K, L> >(B_));
        };

        TEITensor B_ = Slice<Indices<I, I, I, I> >(idxs, f12DelSq);
        TEITensor C_(f12Idxs.nOcc, f12Idxs.nOcc, f12Idxs.nVirt, f12Idxs.nVirt);
        C_.setZero();
        MakeB_C(B_, C_);

        //Calculate conventional MP2 energy
        double E_MP2 = 0.0;

        const auto &epsilon = refWfn.Ep_Alpha;

        //i and j run over the occupied orbitals
        for (Index i = 0; i < f12Idxs.nOcc; i++)
            for (Index j = 0; j < f12Idxs.nOcc; j++) {
                //a and b run over the virtual orbitals
                for (Index a = f12Idxs.nOcc; a < f12Idxs.nVirt + f12Idxs.nOcc; a++)
                    for (Index b = f12Idxs.nOcc; b < f12Idxs.nVirt + f12Idxs.nOcc; b++) {
                        //Integrals are in physicist's notation
                        double intIJAB = eri(i, j, a, b);
                        double intIJBA = eri(i, j, b, a);
                        double epsilonIJAB = epsilon[a] + epsilon[b] - epsilon[i] - epsilon[j];

                        //MP2 formula below already integrates spin out for a closed shell
                        //configuration
                        E_MP2 -= intIJAB * (2.0 * intIJAB - intIJBA) / epsilonIJAB;
                    }
            }


        double E_F12_X = 0.0;
        double E_F12_B = 0.0;
        double E_F12_V = 0.0;

        //Loop over pairs of electrons
        for (Index i = 0; i < f12Idxs.nOcc; i++)
            for (Index j = 0; j < f12Idxs.nOcc; j++) {
                E_F12_V += (5.0 / 8.0) * V_(i, j, i, j);
                E_F12_V -= (1.0 / 8.0) * V_(j, i, i, j);
                E_F12_X += (F_(i, i) + F_(j, j)) * (7.0 / 32.0) * X_(i, j, i, j);
                E_F12_X += (F_(i, i) + F_(j, j)) * (1.0 / 32.0) * X_(i, j, j, i);
                E_F12_B += (7.0 / 32.0) * B_(i, j, i, j);
                E_F12_B += (1.0 / 32.0) * B_(i, j, j, i);
            }


        //Calculate the F12 contribution to the MP2 energy
        double E_F12 = 2.0 * E_F12_V + E_F12_B - E_F12_X;

        if (F_pq_out)
            *F_pq_out = F_;
        if (V_ijxy_out)
            *V_ijxy_out = V_;
        if (X_xywz_out)
            *X_xywz_out = X_;
        if (B_xywz_out)
            *B_xywz_out = B_;
        if (C_xyab_out)
            *C_xyab_out = C_;

        return {E_MP2, E_F12};
    }
    //
    double CalculateMP2_F12_FixedEnergyElectronPair(long i, long j,
                                                    const Tensor2 &F, const TEITensor &V,
                                                    const TEITensor &X, const TEITensor &B) {
        double E_F12_V = (5.0 / 8.0) * V(i, j, i, j);
        E_F12_V -= (1.0 / 8.0) * V(j, i, i, j);

        double E_F12_X = (F(i, i) + F(j, j)) * (7.0 / 32.0) * X(i, j, i, j);
        E_F12_X += (F(i, i) + F(j, j)) * (1.0 / 32.0) * X(i, j, j, i);

        double E_F12_B = (7.0 / 32.0) * B(i, j, i, j);
        E_F12_B += (1.0 / 32.0) * B(i, j, j, i);

        return 2.0 * E_F12_V + E_F12_B - E_F12_X;
    }
    //
    //
}
