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
    //
    using namespace Eigen;
    //
    //
    //-- Indices --
    enum {
        //Occupied orbitals
        I, J, K, L, M, N,
        //OBS Virtual orbitals
        A, B, C, D,
        //All OBS orbitals
        P, Q, R, S, T, U, V, O,
        //CABS virtual orbitals
        X, Y, Z, W
    };
    //
    //
    void EvaluatePerturbative_CC_F12Intermediates(
        const HF_Wfn &refWfn,
        const Integrals::MolInts &mints,
        const Integrals::F12::F12Integrals &f12Ints,
        F12Idxs f12Idxs,
        Tensor2 *F_pq_out,
        TEITensor *V_pqrs_out,
        TEITensor *C_pqrs_out,
        TEITensor *B_pqrs_out,
        TEITensor *X_pqrs_out
    ) {
        TensorSliceIndices idxs = {
            //All Occupied orbital indices
            {I, {0, f12Idxs.nOcc}},
            {J, {0, f12Idxs.nOcc}},
            {K, {0, f12Idxs.nOcc}},
            {L, {0, f12Idxs.nOcc}},
            {M, {0, f12Idxs.nOcc}},
            {N, {0, f12Idxs.nOcc}},
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
            {V, {0, f12Idxs.nOBS}},
            {O, {0, f12Idxs.nOBS}},
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
        TEITensor V_pqrs = Slice<Indices<P, Q, R, S> >(idxs, f12Repulsion);

        //ijrs,pqrs->ijpq
        V_pqrs -= Einsum<Indices<P, Q, T, U>, Indices<R, S, T, U>, Indices<P, Q, R, S> >(
            Slice<Indices<P, Q, R, S> >(idxs, eri),
            Slice<Indices<P, Q, R, S> >(idxs, f12));

        //ijxm,pqxm->ijpq
        V_pqrs -= Einsum<Indices<P, Q, X, M>, Indices<R, S, X, M>, Indices<P, Q, R, S> >(
            Slice<Indices<P, Q, X, M> >(idxs, eri),
            Slice<Indices<P, Q, X, M> >(idxs, f12));

        //ijmx,pqmx->ijpq
        V_pqrs -= Einsum<Indices<P, Q, M, X>, Indices<R, S, M, X>, Indices<P, Q, R, S> >(
            Slice<Indices<P, Q, M, X> >(idxs, eri),
            Slice<Indices<R, S, M, X> >(idxs, f12));


        TEITensor X_pqrs = Slice<Indices<P, Q, R, S> >(idxs, f12Sq);
        //ijrs,klrs->ijkl
        X_pqrs -= Einsum<Indices<P, Q, T, U>, Indices<R, S, T, U>, Indices<P, Q, R, S> >(
            Slice<Indices<P, Q, R, S> >(idxs, f12),
            Slice<Indices<T, U, R, S> >(idxs, f12));

        //ijxm,klxm->ijkl
        X_pqrs -= Einsum<Indices<P, Q, X, M>, Indices<R, S, X, M>, Indices<P, Q, R, S> >(
            Slice<Indices<P, Q, X, M> >(idxs, f12),
            Slice<Indices<R, S, X, M> >(idxs, f12));

        //ijmx,klmx->ijkl
        X_pqrs -= Einsum<Indices<P, Q, M, X>, Indices<R, S, M, X>, Indices<P, Q, R, S> >(
            Slice<Indices<P, Q, M, X> >(idxs, f12),
            Slice<Indices<R, S, M, X> >(idxs, f12));

        //Full Fock Matrix
        Tensor<double, 2> F_pq = MatrixToTensor(TInts + VInts);
        //Coloumb Fock Terms
        //piqi->pq
        Slice<Indices<P, Q> >(idxs, F_pq) += 2.0 * Einsum<Indices<P, I, Q, I>, Indices<P, Q> >(
            Slice<Indices<P, I, Q, I> >(idxs, eri));
        //ipix->px
        Slice<Indices<P, X> >(idxs, F_pq) += 2.0 * Einsum<Indices<I, P, I, X>, Indices<P, X> >(
            Slice<Indices<I, P, I, X> >(idxs, eri));
        //ixiy->xy
        Slice<Indices<X, Y> >(idxs, F_pq) += 2.0 * Einsum<Indices<I, X, I, Y>, Indices<X, Y> >(
            Slice<Indices<I, X, I, Y> >(idxs, eri));

        Tensor<double, 2> K_(F_pq.dimensions());
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
        F_pq -= K_;

        //Fill in the remaining elements by symmetry
        Slice<Indices<X, P> >(idxs, F_pq) = Einsum<Indices<P, X>, Indices<X, P> >(Slice<Indices<P, X> >(idxs, F_pq));
        Slice<Indices<X, P> >(idxs, K_) = Einsum<Indices<P, X>, Indices<X, P> >(Slice<Indices<P, X> >(idxs, K_));

        //Fock matrix without the exchange terms
        Tensor<double, 2> F_K = F_pq + K_;

        //Makes the B and C intermediate tensors
        TEITensor B_pqrs = Slice<Indices<P, P, P, P> >(idxs, f12DelSq);
        //-- F12 Contributions --
        //pr,lkqr-> plkq -> klpq
        TEITensor tmp = Einsum<Indices<P, R>, Indices<U, T, Q, R>, Indices<T, U, P, Q> >(
            Slice<Indices<P, R> >(idxs, F_K),
            Slice<Indices<U, T, Q, R> >(idxs, f12));
        //mnpq,klpq-> mnkl -> klmn
        B_pqrs -= Einsum<Indices<V, O, P, Q>, Indices<T, U, P, Q>, Indices<T, U, V, O> >(
            Slice<Indices<V, O, P, Q> >(idxs, f12),
            tmp);

        //qr,klpr-> qklp -> klpq
        tmp = Einsum<Indices<Q, R>, Indices<T, U, P, R>, Indices<T, U, P, Q> >(
            Slice<Indices<Q, R> >(idxs, F_K),
            Slice<Indices<T, U, P, R> >(idxs, f12));
        //mnpq,klpq-> mnkl -> klmn
        B_pqrs -= Einsum<Indices<V, O, P, Q>, Indices<T, U, P, Q>, Indices<T, U, V, O> >(
            Slice<Indices<V, O, P, Q> >(idxs, f12),
            tmp);

        //px,lkqx-> plkq -> klpq
        tmp = Einsum<Indices<P, X>, Indices<U, T, Q, X>, Indices<T, U, P, Q> >(
            Slice<Indices<P, X> >(idxs, F_K),
            Slice<Indices<U, T, Q, X> >(idxs, f12));
        //mnpq,klpq-> mnkl -> klmn
        B_pqrs -= Einsum<Indices<V, O, P, Q>, Indices<T, U, P, Q>, Indices<T, U, V, O> >(
            Slice<Indices<V, O, P, Q> >(idxs, f12),
            tmp);

        //qx,klpx-> qklp -> klpq
        tmp = Einsum<Indices<Q, X>, Indices<T, U, P, X>, Indices<T, U, P, Q> >(
            Slice<Indices<Q, X> >(idxs, F_K),
            Slice<Indices<T, U, P, X> >(idxs, f12));
        //mnpq,klpq-> mnkl -> klmn
        B_pqrs -= Einsum<Indices<V, O, P, Q>, Indices<T, U, P, Q>, Indices<T, U, V, O> >(
            Slice<Indices<V, O, P, Q> >(idxs, f12),
            tmp);

        //yp,klip-> ykli ->  kliy
        tmp = Einsum<Indices<Y, P>, Indices<T, U, I, P>, Indices<T, U, I, Y> >(
            Slice<Indices<Y, P> >(idxs, F_K),
            Slice<Indices<T, U, I, P> >(idxs, f12));
        //mniy,kliy-> mnkl -> klmn
        B_pqrs -= Einsum<Indices<V, O, I, Y>, Indices<T, U, I, Y>, Indices<T, U, V, O> >(
            Slice<Indices<V, O, I, Y> >(idxs, f12),
            tmp);

        //xp,lkjp-> xlkj -> klxj
        tmp = Einsum<Indices<X, P>, Indices<U, T, J, P>, Indices<T, U, X, J> >(
            Slice<Indices<X, P> >(idxs, F_K),
            Slice<Indices<U, T, J, P> >(idxs, f12));
        //nmjx,klxj-> nmkl ->  klmn
        B_pqrs -= Einsum<Indices<O, V, J, X>, Indices<T, U, X, J>, Indices<T, U, V, O> >(
            Slice<Indices<O, V, J, X> >(idxs, f12),
            tmp);

        //ip,klpy-> ikly -> kliy
        tmp = Einsum<Indices<I, P>, Indices<T, U, P, Y>, Indices<T, U, I, Y> >(
            Slice<Indices<I, P> >(idxs, F_K),
            Slice<Indices<T, U, P, Y> >(idxs, f12));
        //mniy,kliy-> mnkl -> klmn
        B_pqrs -= Einsum<Indices<V, O, I, Y>, Indices<T, U, I, Y>, Indices<T, U, V, O> >(
            Slice<Indices<V, O, I, Y> >(idxs, f12),
            tmp);

        //yx,klix->kliy
        tmp = Einsum<Indices<Y, X>, Indices<T, U, I, X>, Indices<T, U, I, Y> >(
            Slice<Indices<Y, X> >(idxs, F_K),
            Slice<Indices<T, U, I, X> >(idxs, f12));
        //mniy,kliy->klmn
        B_pqrs -= Einsum<Indices<V, O, I, Y>, Indices<T, U, I, Y>, Indices<T, U, V, O> >(
            Slice<Indices<V, O, I, Y> >(idxs, f12),
            tmp);

        //xy,lkjy->klxj
        tmp = Einsum<Indices<X, Y>, Indices<U, T, J, Y>, Indices<T, U, X, J> >(
            Slice<Indices<X, Y> >(idxs, F_K),
            Slice<Indices<U, T, J, Y> >(idxs, f12));
        //nmjx,klxj->klmn
        B_pqrs -= Einsum<Indices<O, V, J, X>, Indices<T, U, X, J>, Indices<T, U, V, O> >(
            Slice<Indices<O, V, J, X> >(idxs, f12),
            tmp);

        //jp,lkpx->klxj
        tmp = Einsum<Indices<J, P>, Indices<U, T, P, X>, Indices<T, U, X, J> >(
            Slice<Indices<J, P> >(idxs, F_K),
            Slice<Indices<U, T, P, X> >(idxs, f12));
        //nmjx,klxj->klmn
        B_pqrs -= Einsum<Indices<O, V, J, X>, Indices<T, U, X, J>, Indices<T, U, V, O> >(
            Slice<Indices<O, V, J, X> >(idxs, f12),
            tmp);

        //ix,klxy->kliy
        tmp = Einsum<Indices<I, X>, Indices<T, U, X, Y>, Indices<T, U, I, Y> >(
            Slice<Indices<I, X> >(idxs, F_K),
            Slice<Indices<T, U, X, Y> >(idxs, f12));
        //mniy,kliy->klmn
        B_pqrs -= Einsum<Indices<V, O, I, Y>, Indices<T, U, I, Y>, Indices<T, U, V, O> >(
            Slice<Indices<V, O, I, Y> >(idxs, f12),
            tmp);

        //jy,lkyx->klxj
        tmp = Einsum<Indices<J, Y>, Indices<U, T, Y, X>, Indices<T, U, X, J> >(
            Slice<Indices<J, Y> >(idxs, F_K),
            Slice<Indices<U, T, Y, X> >(idxs, f12));
        //nmjx,klxj->klmn
        B_pqrs -= Einsum<Indices<O, V, J, X>, Indices<T, U, X, J>, Indices<T, U, V, O> >(
            Slice<Indices<O, V, J, X> >(idxs, f12),
            tmp);

        //-- Y contributions --
        //xp,lkbp->klxb
        tmp = Einsum<Indices<X, P>, Indices<U, T, B, P>, Indices<T, U, X, B> >(
            Slice<Indices<X, P> >(idxs, K_),
            Slice<Indices<U, T, B, P> >(idxs, f12));
        //klxb,nmbx->klmn
        B_pqrs -= Einsum<Indices<T, U, X, B>, Indices<O, V, B, X>, Indices<T, U, V, O> >(
            tmp,
            Slice<Indices<O, V, B, X> >(idxs, f12));

        //yp,klap->klay
        tmp = Einsum<Indices<Y, P>, Indices<T, U, A, P>, Indices<T, U, A, Y> >(
            Slice<Indices<Y, P> >(idxs, K_),
            Slice<Indices<T, U, A, P> >(idxs, f12));
        //klay,mnay->klmn
        B_pqrs -= Einsum<Indices<T, U, A, Y>, Indices<V, O, A, Y>, Indices<T, U, V, O> >(
            tmp,
            Slice<Indices<V, O, A, Y> >(idxs, f12));

        //xy,lkby->klxb
        tmp = Einsum<Indices<X, Y>, Indices<U, T, B, Y>, Indices<T, U, X, B> >(
            Slice<Indices<X, Y> >(idxs, K_),
            Slice<Indices<U, T, B, Y> >(idxs, f12));
        //klxb,nmbx->klmn
        B_pqrs -= Einsum<Indices<T, U, X, B>, Indices<O, V, B, X>, Indices<T, U, V, O> >(
            tmp,
            Slice<Indices<O, V, B, X> >(idxs, f12));

        //bp,lkpx->klxb
        tmp = Einsum<Indices<B, P>, Indices<U, T, P, X>, Indices<T, U, X, B> >(
            Slice<Indices<B, P> >(idxs, K_),
            Slice<Indices<U, T, P, X> >(idxs, f12));
        //klxb,nmbx->klmn
        B_pqrs -= Einsum<Indices<T, U, X, B>, Indices<O, V, B, X>, Indices<T, U, V, O> >(
            tmp,
            Slice<Indices<O, V, B, X> >(idxs, f12));

        //ap,klpy->klay
        tmp = Einsum<Indices<A, P>, Indices<T, U, P, Y>, Indices<T, U, A, Y> >(
            Slice<Indices<A, P> >(idxs, K_),
            Slice<Indices<T, U, P, Y> >(idxs, f12));
        //klay,mnay->klmn
        B_pqrs -= Einsum<Indices<T, U, A, Y>, Indices<V, O, A, Y>, Indices<T, U, V, O> >(
            tmp,
            Slice<Indices<V, O, A, Y> >(idxs, f12));

        //yx,klax->klay
        tmp = Einsum<Indices<Y, X>, Indices<T, U, A, X>, Indices<T, U, A, Y> >(
            Slice<Indices<Y, X> >(idxs, K_),
            Slice<Indices<T, U, A, X> >(idxs, f12));
        //klay,mnay->klmn
        B_pqrs -= Einsum<Indices<T, U, A, Y>, Indices<V, O, A, Y>, Indices<T, U, V, O> >(
            tmp,
            Slice<Indices<V, O, A, Y> >(idxs, f12));


        //by,lkyx->klxb
        tmp = Einsum<Indices<B, Y>, Indices<U, T, Y, X>, Indices<T, U, X, B> >(
            Slice<Indices<B, Y> >(idxs, K_),
            Slice<Indices<U, T, Y, X> >(idxs, f12));
        //klxb,nmbx->klmn
        B_pqrs -= Einsum<Indices<T, U, X, B>, Indices<O, V, B, X>, Indices<T, U, V, O> >(
            tmp,
            Slice<Indices<O, V, B, X> >(idxs, f12));

        //ax,klxy->klay
        tmp = Einsum<Indices<A, X>, Indices<T, U, X, Y>, Indices<T, U, A, Y> >(
            Slice<Indices<A, X> >(idxs, K_),
            Slice<Indices<T, U, X, Y> >(idxs, f12));
        //klay,mnay->klmn
        B_pqrs -= Einsum<Indices<T, U, A, Y>, Indices<V, O, A, Y>, Indices<T, U, V, O> >(
            tmp,
            Slice<Indices<V, O, A, Y> >(idxs, f12));

        //xp,klpy->klxy
        tmp = Einsum<Indices<X, P>, Indices<T, U, P, Y>, Indices<T, U, X, Y> >(
            Slice<Indices<X, P> >(idxs, K_),
            Slice<Indices<T, U, P, Y> >(idxs, f12));
        //klxy,mnxy->klmn
        B_pqrs -= Einsum<Indices<T, U, X, Y>, Indices<V, O, X, Y>, Indices<T, U, V, O> >(
            tmp,
            Slice<Indices<V, O, X, Y> >(idxs, f12));

        //yp,lkpx->klxy
        tmp = Einsum<Indices<Y, P>, Indices<U, T, P, X>, Indices<T, U, X, Y> >(
            Slice<Indices<Y, P> >(idxs, K_),
            Slice<Indices<U, T, P, X> >(idxs, f12));
        //klxy,mnxy->klmn
        B_pqrs -= Einsum<Indices<T, U, X, Y>, Indices<V, O, X, Y>, Indices<T, U, V, O> >(
            tmp,
            Slice<Indices<V, O, X, Y> >(idxs, f12));

        //yz,lkzx->klxy
        tmp = Einsum<Indices<Y, Z>, Indices<U, T, Z, X>, Indices<T, U, X, Y> >(
            Slice<Indices<Y, Z> >(idxs, K_),
            Slice<Indices<U, T, Z, X> >(idxs, f12));
        //klxy,mnxy->klmn
        B_pqrs -= Einsum<Indices<T, U, X, Y>, Indices<V, O, X, Y>, Indices<T, U, V, O> >(
            tmp,
            Slice<Indices<V, O, X, Y> >(idxs, f12));

        //xz,klzy->klxy
        tmp = Einsum<Indices<X, Z>, Indices<T, U, Z, Y>, Indices<T, U, X, Y> >(
            Slice<Indices<X, Z> >(idxs, K_),
            Slice<Indices<T, U, Z, Y> >(idxs, f12));
        //klxy,mnxy->klmn
        B_pqrs -= Einsum<Indices<T, U, X, Y>, Indices<V, O, X, Y>, Indices<T, U, V, O> >(
            tmp,
            Slice<Indices<V, O, X, Y> >(idxs, f12));


        //--F Squared Contributions--
        //kp,nmlp->klmn
        B_pqrs += Einsum<Indices<T, P>, Indices<O, V, U, P>, Indices<T, U, V, O> >(
            Slice<Indices<T, P> >(idxs, F_K),
            Slice<Indices<O, V, U, P> >(idxs, f12Sq));
        //lp,mnkp->klmn
        B_pqrs += Einsum<Indices<U, P>, Indices<V, O, T, P>, Indices<T, U, V, O> >(
            Slice<Indices<U, P> >(idxs, F_K),
            Slice<Indices<V, O, T, P> >(idxs, f12Sq));
        //kx,nmlx->klmn
        B_pqrs += Einsum<Indices<T, X>, Indices<O, V, U, X>, Indices<T, U, V, O> >(
            Slice<Indices<T, X> >(idxs, F_K),
            Slice<Indices<O, V, U, X> >(idxs, f12Sq));
        //lx,mnkx->klmn
        B_pqrs += Einsum<Indices<U, X>, Indices<V, O, T, X>, Indices<T, U, V, O> >(
            Slice<Indices<U, X> >(idxs, F_K),
            Slice<Indices<V, O, T, X> >(idxs, f12Sq));

        //-- C Tensor --
        //sx,pqrx-pqrs
        TEITensor C_pqrs = Einsum<Indices<S, X>, Indices<P, Q, R, X>, Indices<P, Q, R, S> >(
            Slice<Indices<S, X> >(idxs, F_pq),
            Slice<Indices<P, Q, R, X> >(idxs, f12));
        //rx,qpsx->pqrs
        C_pqrs += Einsum<Indices<R, X>, Indices<Q, P, S, X>, Indices<P, Q, R, S> >(
            Slice<Indices<R, X> >(idxs, F_pq),
            Slice<Indices<Q, P, S, X> >(idxs, f12));


        //-- -FC contribution to B--
        //pqtu,rstu->pqrs
        B_pqrs -= Einsum<Indices<P, Q, T, U>, Indices<R, S, T, U>, Indices<P, Q, R, S> >(
            Slice<Indices<P, Q, R, S> >(idxs, f12),
            C_pqrs);

        //--Symmetrise B--
        B_pqrs = 0.5 * (B_pqrs + Einsum<Indices<P, Q, R, S>, Indices<R, S, P, Q> >(B_pqrs));

        if (F_pq_out)
            *F_pq_out = F_pq;
        if (V_pqrs_out)
            *V_pqrs_out = V_pqrs;
        if (C_pqrs_out)
            *C_pqrs_out = C_pqrs;
        if (B_pqrs_out)
            *B_pqrs_out = B_pqrs;
        if (X_pqrs_out)
            *X_pqrs_out = X_pqrs;
    }
    //
    //
}
