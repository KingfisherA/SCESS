#include<Quantum/QuantumPCH.h>
#include "F12.h"

#include "Tensor/TensorOperations.h"
#include "Tensor/TensorSlice.h"
#include "Tensor/TensorEinsum.h"
#include "Quantum/Integrals/Integrals.h"
#include "Quantum/Integrals/IntegralsF12.h"
#include "Quantum/Integrals/IntegralTransforms.h"
#include "Quantum/HF/HartreeFock.h"

namespace SCESS::F12 {
    using namespace Eigen;

    // Tensor4 EvaluateVIntermediate(const TEITensor &f12Repulsion,
    //                             const TEITensor &f12,
    //                             const TEITensor &eri) {
    //     //Make the V intermediate
    //     //iiii
    //     TEITensor V_ = Slice<Indices<I, I, I, I> >(idxs, f12Repulsion);
    //
    //     //ijrs,klrs->ijkl
    //     V_ -= Einsum<Indices<I, J, R, S>, Indices<K, L, R, S>, Indices<I, J, K, L> >(
    //         Slice<Indices<I, J, R, S> >(idxs, eri),
    //         Slice<Indices<K, L, R, S> >(idxs, f12));
    //
    //     //ijxm,klxm->ijkl
    //     V_ -= Einsum<Indices<I, J, X, M>, Indices<K, L, X, M>, Indices<I, J, K, L> >(
    //         Slice<Indices<I, J, X, M> >(idxs, eri),
    //         Slice<Indices<K, L, X, M> >(idxs, f12));
    //
    //     //ijmx,klmx->ijkl
    //     V_ -= Einsum<Indices<I, J, M, X>, Indices<K, L, M, X>, Indices<I, J, K, L> >(
    //         Slice<Indices<I, J, M, X> >(idxs, eri),
    //         Slice<Indices<K, L, M, X> >(idxs, f12));
    //
    //     return V_;
    // }

    MP2::MP2F12Result MP2::MP2_F12_Optimise(const UHFWavefunction &refWfn,
                                   const Integrals::MolInts &mints,
                                   const Integrals::F12::F12Integrals &f12Ints) {
        // * We assume that the reference wavefunction is actually a closed shell RHF hence alternating
        //      alpha and beta spins (fix this in the future to support open shell :( )
        // * We need | Psi > = |ijk...> |alpha beta alpha beta ... > (alpha and beta are stored alternating)
        // * All calculations are done in SPIN ORBITALS so we need to integrate spin out the spin!
        // * All EVEN indices are spin alpha and ODD are spin beta

        //-- Indices --
        enum {
            //Occupied orbitals
            I, J, K, L, M, N, O,
            //OBS Virtual orbitals
            A, B, C, D,
            //All OBS orbitals
            P, Q, R, S,
            //CABS virtual orbitals
            X, Y, Z, W
        };

        //Number of orbitals in the ordinary basis set. Also the index for the beginning of the CABS
        constexpr Index nOBS = 2;
        //Number of orbitals in the complementary auxiliary set
        constexpr Index nCABS = 10;
        //Total number of basis functions in the combined OBS and CABS basis sets
        constexpr Index nbf = nOBS + nCABS;
        //Number of spatial molecular orbitals equals the number of OBS used in normal Hartree-Fock
        constexpr Index nMO = nOBS;
        //Only one occupied spatial MO in hydrogen
        constexpr Index nVirt = nMO - 1;
        //Two spin orbitals are occupied in hydrogen and one spatial orbital
        const Index nOcc = 1;
        // std::distance(refWfn.OccupationAlpha.begin(),
        //                              std::ranges::find(refWfn.OccupationAlpha, 0));;
        TensorSliceIndices idxs = {
            //All Occupied orbital indices
            {I, {0, nOcc}},
            {J, {0, nOcc}},
            {K, {0, nOcc}},
            {L, {0, nOcc}},
            {M, {0, nOcc}},
            {N, {0, nOcc}},
            {O, {0, nOcc}},
            //All virtual orbital indices
            {A, {nOcc, nVirt}},
            {B, {nOcc, nVirt}},
            {C, {nOcc, nVirt}},
            {D, {nOcc, nVirt}},
            //All OBS indices
            {P, {0, nOBS}},
            {Q, {0, nOBS}},
            {R, {0, nOBS}},
            {S, {0, nOBS}},
            //CABS indices
            {X, {nOBS, nCABS}},
            {Y, {nOBS, nCABS}},
            {Z, {nOBS, nCABS}},
            {W, {nOBS, nCABS}}
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
        //, &idxs, &nOcc, &nOBS, &nCABS]
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

            //-- Y contributions --
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
        TEITensor C_(nOcc, nOcc, nVirt, nVirt);
        C_.setZero();
        MakeB_C(B_, C_);

        //Calculate conventional MP2 energy
        double E_MP2 = 0.0;

        const auto &epsilon = refWfn.EnergiesAlpha;

        //i and j run over the occupied orbitals
        for (Index i = 0; i < nOcc; i++)
            for (Index j = 0; j < nOcc; j++) {
                //a and b run over the virtual orbitals
                for (Index a = nOcc; a < nVirt + nOcc; a++)
                    for (Index b = nOcc; b < nVirt + nOcc; b++) {
                        //Integrals are in physicist's notation
                        double intIJAB = eri(i, j, a, b);
                        double intIJBA = eri(i, j, b, a);
                        double epsilonIJAB = epsilon[a] + epsilon[b] - epsilon[i] - epsilon[j];

                        //MP2 formula below already integrates spin out for a closed shell
                        //configuration
                        E_MP2 -= intIJAB * (2.0 * intIJAB - intIJBA) / epsilonIJAB;
                    }
            }


        //Makes the V tilde intermediate integrating out spin for a closed shell reference
        //&nOcc,&nVirt,
        auto MakeVTildeSpin = [&eri, &V_, &C_, &F_](MatrixXd &V_AlphaAlpha,
                                                    MatrixXd &V_AlphaBeta,
                                                    Index i, Index j) {
            TEITensor eri_oovv = eri.Slice4({0, 0, nOcc, nOcc},
                                            {nOcc, nOcc, nVirt, nVirt});
            Index col = 0;
            for (Index k = 0; k < nOcc; k++)
                for (Index l = 0; l < nOcc; l++, col++) {
                    double VAlphaBeta_kl = V_(k, l, i, j);
                    double VAlphaAlpha_kl = V_(k, l, i, j) - V_(l, k, i, j);
                    for (Index a = 0; a < nVirt; a++)
                        for (Index b = 0; b < nVirt; b++) {
                            double D_ijab = F_(nOcc + a, nOcc + a) + F_(nOcc + b, nOcc + b)
                                            - F_(i, i) - F_(j, j);
                            double C_klab = C_(k, l, a, b);
                            double C_klba = C_(k, l, b, a);
                            double eri_ijab = eri_oovv(i, j, a, b);
                            double eri_ijba = eri_oovv(i, j, b, a);

                            VAlphaBeta_kl -= C_klab * eri_ijab / D_ijab;

                            //If there are two alpha or two beta spin functions
                            //exchange terms need to be considered
                            if (k < l && a < b) {
                                C_klab = C_klab - C_klba;
                                eri_ijab = eri_ijab - eri_ijba;
                                VAlphaAlpha_kl -= C_klab * eri_ijab / D_ijab;
                            }
                        }

                    V_AlphaBeta(0, col) = VAlphaBeta_kl;
                    V_AlphaAlpha(0, col) = VAlphaAlpha_kl;
                }

            V_AlphaBeta = V_AlphaBeta.transpose();
            V_AlphaAlpha = V_AlphaAlpha.transpose();
        };

        //Makes the B tilde intermediate integrating out spin for a closed shell reference
        //&nOcc
        auto MakeBTildeSpin = [&F_, &C_, &X_, &B_](MatrixXd &B_AlphaAlpha,
                                                          MatrixXd &B_AlphaBeta,
                                                          Index i, Index j) {
            double epsilon_ij = F_(i, i) + F_(j, j);
            Index row = 0, col = 0;
            for (Index k = 0; k < nOcc; k++)
                for (Index l = 0; l < nOcc; l++, row++) {
                    for (Index m = 0; m < nOcc; m++)
                        for (Index n = 0; n < nOcc; n++, col++) {
                            double B_AlphaBeta_klmn = B_(k, l, m, n)
                                                      - epsilon_ij * X_(k, l, m, n);
                            double B_AlphaAlpha_klmn = 0.0;

                            if (k < l && m < n)
                                B_AlphaAlpha_klmn =
                                        B_(k, l, m, n)
                                        - B_(k, l, n, m)
                                        - epsilon_ij * X_(k, l, m, n)
                                        + epsilon_ij * X_(k, l, n, m);

                            for (Index a = 0; a < nVirt; a++)
                                for (Index b = 0; b < nVirt; b++) {
                                    double D_ijab = F_(nOcc + a, nOcc + a) + F_(nOcc + b, nOcc + b) - F_(
                                                        i, i) - F_(j, j);
                                    double C_klab = C_(k, l, a, b);
                                    double C_klba = C_(k, l, b, a);
                                    double C_mnab = C_(m, n, a, b);
                                    double C_mnba = C_(m, n, b, a);

                                    B_AlphaBeta_klmn -= C_klab * C_mnab / D_ijab;

                                    if (a < b) {
                                        //Antisymmetrise matrix elements
                                        C_klab = C_klab - C_klba;
                                        C_mnab = C_mnab - C_mnba;
                                        B_AlphaAlpha_klmn -= C_klab * C_mnab / D_ijab;
                                    }
                                }

                            B_AlphaBeta(row, col) = B_AlphaBeta_klmn;
                            B_AlphaAlpha(row, col) = B_AlphaAlpha_klmn;
                        }
                }
        };

        //Calculate the F12 contribution to the MP2 energy
        double E_F12 = 0.0;

        //Loop over pairs of electrons
        for (Index i = 0; i < nOcc; i++)
            for (Index j = 0; j < nOcc; j++) {
                MatrixXd V_AlphaAlpha = MatrixXd::Zero(1, nOcc * nOcc);
                MatrixXd V_AlphaBeta = MatrixXd::Zero(1, nOcc * nOcc);
                MatrixXd B_AlphaAlpha = MatrixXd::Zero(nOcc * nOcc, nOcc * nOcc * nVirt * nVirt);
                MatrixXd B_AlphaBeta = MatrixXd::Zero(nOcc * nOcc, nOcc * nOcc * nVirt * nVirt);

                MakeVTildeSpin(V_AlphaAlpha, V_AlphaBeta, i, j);
                MakeBTildeSpin(B_AlphaAlpha, B_AlphaBeta, i, j);

                MatrixXd B_AlphaBetaInv = B_AlphaBeta.inverse();

                double AlphaBetaPairEnergy = -(V_AlphaBeta.transpose() * B_AlphaBetaInv * V_AlphaBeta)(0, 0);
                double AlphaAlphaPairEnergy = 0.0;

                //Given our wavefunction is alpha beta alpha beta ect. if i < j then we
                //have two alpha or two beta spin functions so we need to have an exchange
                //term
                if (i < j) {
                    MatrixXd BAlphaAlphaInv = B_AlphaAlpha.inverse();
                    AlphaAlphaPairEnergy = -(V_AlphaAlpha.transpose() * BAlphaAlphaInv * V_AlphaAlpha)(0, 0);
                }

                E_F12 += AlphaBetaPairEnergy + 2.0 * AlphaAlphaPairEnergy;
            }


        return {E_MP2, E_F12};
    }
}
