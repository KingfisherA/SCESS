#include<Quantum/QuantumPCH.h>
#include"F12.h"

#include"Tensor/TensorOperations.h"
#include"Tensor/TensorSlice.h"
#include "Tensor/TensorEinsum.h"

namespace SCESS::F12::QMC {
    using namespace Eigen;

    Tensor4 EvaluateW(const Tensor4 &f12, const Tensor2 &S) {
        constexpr Index nOcc = 1;
        constexpr Index nVirt = 1;
        constexpr Index nCABS = 10;

        enum { alpha, beta, I, J, M, N, A, B };

        TensorSliceIndices idxs = {
            {alpha, {nOcc, nVirt + nCABS}},
            {beta, {nOcc, nVirt + nCABS}},
            //
            {I, {0, nOcc}},
            {J, {0, nOcc}},
            {M, {0, nOcc}},
            {N, {0, nOcc}},
            //
            {A, {nOcc, nVirt}},
            {B, {nOcc, nVirt}}
        };

        Tensor4 W = Slice<Indices<alpha, beta, I, J> >(idxs, f12);
        //alpha m,beta m ij -> alpha beta ij
        W -= Einsum<Indices<alpha, M>, Indices<beta, M, I, J>, Indices<alpha, beta, I, J> >(
            Slice<Indices<alpha, M> >(idxs, S),
            Slice<Indices<beta, M, I, J> >(idxs, f12));
        // .Contract1(
        //     f12.Slice4({beta.first, 0, i.first, j.first},
        //         {beta.second, nOcc, i.second, j.second}), ({{1,1}}));


        //alpha m ij, beta m -> alpha ij beta -> alpha beta ij
        W -= Einsum<Indices<alpha, M, I, J>, Indices<beta, M>, Indices<alpha, beta, I, J> >(
            Slice<Indices<alpha, M, I, J> >(idxs, f12),
            Slice<Indices<beta, M> >(idxs, S));
        // W -= f12.Slice4({alpha.first, 0, i.first, j.first},
        //                 {alpha.second, nOcc, i.second, j.second}).Contract1(
        //             S.Slice2({beta.first, 0},{beta.second, nOcc}), ({{1,1}}))
        //         .Shuffle4({0, 3, 1, 2});

        //beta n, alpha m, mnij -> alpha beta ij
        W += Einsum<Indices<alpha, M>, Indices<beta, M, I, J>, Indices<alpha, beta, I, J> >(
            Slice<Indices<alpha, M> >(idxs, S),
            Einsum<Indices<beta, N>, Indices<M, N, I, J>, Indices<beta, M, I, J> >(
                Slice<Indices<beta, N> >(idxs, S),
                Slice<Indices<M, N, I, J> >(idxs, f12)));
        // W += S.Slice2({alpha.first, 0}, {alpha.second, nOcc})
        //         .Contract1(
        //             S.Slice2({alpha.first, 0}, {alpha.second, nOcc})
        //             .Contract1(
        //                 f12.Slice4({0,0,0,0}, {nOcc, nOcc, nOcc, nOcc}),
        //                 ({{1,0}})),
        //             ({{1,1}}));

        //beta b, alpha a, abij -> alpha beta ij
        W += Einsum<Indices<alpha, A>, Indices<beta, A, I, J>, Indices<alpha, beta, I, J> >(
            Slice<Indices<alpha, A> >(idxs, S),
            Einsum<Indices<beta, B>, Indices<A, B, I, J>, Indices<beta, A, I, J> >(
                Slice<Indices<beta, B> >(idxs, S),
                Slice<Indices<A, B, I, J> >(idxs, f12)));
        // W += S.Slice2({beta.first, nOcc}, {beta.second, nVirt})
        //         .Contract1(
        //             S.Slice2({alpha.first, nOcc}, {alpha.second,nVirt})
        //             .Contract1(
        //                 f12.Slice4({nOcc,nOcc,0,0}, {nVirt, nVirt, nOcc, nOcc}),
        //                 ({{1,0}})),
        //             ({{1,1}}));

        //alpha beta ij -> alpha beta ji
        Tensor4 W_tmp = Einsum<Indices<alpha, beta, I, J>, Indices<alpha, beta, J, I> >(W);

        W = (3.0 / 8.0) * W;
        W += (1.0 / 8.0) * W_tmp;

        return W;
    }
}
