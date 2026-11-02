#include<Quantum/Q_PCH.h>
#include"F12.h"
//
#include<Eigen/Core>
#include<unsupported/Eigen/CXX11/Tensor>
//
#include "Tensor/TensorOperations.h"
#include "Tensor/TensorSlice.h"
#include "Tensor/TensorEinsum.h"
#include "Utilities/Utilities.h"

//
//
namespace SCESS::F12 {
    using namespace Eigen;
    //
    //
    double EvaluateCABSSinglesCorrection(const Tensor2 &F,
                                           F12Idxs f12Idxs,
                                           Tensor2 *t1Amps) {
        enum {
            //All OBS orbitals
            P, Q,
            //Occupied orbitals
            I, J, K,
            //Virtual orbitals
            A, B,
            //Virtual CABS orbitals
            X, Y
        };
        TensorSliceIndices idxs = {
            //All Occupied orbital indices
            {I, {0, f12Idxs.nOcc}},
            {J, {0, f12Idxs.nOcc}},
            {K, {0, f12Idxs.nOcc}},
            //All virtual orbital indices
            {A, {f12Idxs.nOcc, f12Idxs.nVirt}},
            {B, {f12Idxs.nOcc, f12Idxs.nVirt}},
            //All OBS indices
            {P, {0, f12Idxs.nOBS}},
            {Q, {0, f12Idxs.nOBS}},
            //CABS indices
            {X, {f12Idxs.nOBS, f12Idxs.nCABS}},
            {Y, {f12Idxs.nOBS, f12Idxs.nCABS}}
        };


        //We need to assume canonical orbitals to solve equations iteratively
        Tensor2 Dai(f12Idxs.nVirt, f12Idxs.nOcc);
        Dai.setZero();

        for (Eigen::Index i = 0; i < f12Idxs.nOcc; i++)
            for (Eigen::Index a = 0; a < f12Idxs.nVirt; a++)
                Dai(a, i) = F(i, i) - F(a + f12Idxs.nOcc, a + f12Idxs.nOcc);


        Tensor2 Dxi(f12Idxs.nCABS, f12Idxs.nOcc);
        Dxi.setZero();

        for (Eigen::Index x = 0; x < f12Idxs.nCABS; x++)
            for (Eigen::Index i = 0; i < f12Idxs.nOcc; i++)
                Dxi(x, i) = F(i, i) - F(x + f12Idxs.nOBS, x + f12Idxs.nOBS);


        //Inital guess
        Tensor2 t_xi = Slice<Indices<X, I> >(idxs, F);
        Tensor2 t_ai(f12Idxs.nVirt, f12Idxs.nOcc);
        t_ai.setZero();

        Tensor2 t1Last(t_xi.dimensions());
        t1Last.setZero();

        constexpr uint32_t MAX_ITERATIONS = 100;
        double delta = std::numeric_limits<double>::infinity();

        uint32_t iteration = 0;
        for (; iteration < MAX_ITERATIONS; iteration++) {
            //ax, xi -> ai
            Tensor2 rhs_ai = Einsum<Indices<A, X>, Indices<X, I>, Indices<A, I> >(
                Slice<Indices<A, X> >(idxs, F),
                t_xi);

            Tensor2 rhs_xi = Slice<Indices<X, I> >(idxs, F);

            //xa,ai->xi
            rhs_xi += Einsum<Indices<X, A>, Indices<A, I>, Indices<X, I> >(
                Slice<Indices<X, A> >(idxs, F),
                t_ai);

            t_ai = rhs_ai / Dai;
            t_xi = rhs_xi / Dxi;

            //Largest change in any amplitude, so that changes of opposite sign cannot cancel
            Tensor2 deltaTensor = t_xi - t1Last;
            Tensor<double, 0> deltaTensor2 = deltaTensor.abs().maximum();
            delta = deltaTensor2(0);

            if (delta < 1E-8)
                break;

            t1Last = t_xi;
        }

        if (delta >= 1E-8)
            std::cerr << "Warning: CABS singles amplitudes failed to converge after " << MAX_ITERATIONS
                    << " iterations (max change " << delta << ")\n";

        //Times two for pair of electrons
        //ix, xi ->
        Tensor<double, 0> ECABSSingles = 2.0 * Einsum<Indices<I, X>, Indices<X, I>, Indices<> >(
                                             Slice<Indices<I, X> >(idxs, F),
                                             t_xi);
        if (t1Amps)
            *t1Amps = t_xi;

        return ECABSSingles(0);
    }
    //
    //
}
