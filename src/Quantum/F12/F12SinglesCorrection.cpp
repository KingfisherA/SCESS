#include<Quantum/QuantumPCH.h>
#include"F12.h"

#include<Eigen/Core>
#include<unsupported/Eigen/CXX11/Tensor>

#include "Tensor/TensorOperations.h"
#include "Tensor/TensorSlice.h"
#include "Tensor/TensorEinsum.h"
#include "Tensor/TensorEinsum.h"
#include "Utilities/Utilities.h"

namespace SCESS::F12 {
    using namespace Eigen;

    double EvaluateSinglesCorrectionEnergy(const Tensor2 &F, Tensor2 *t1Amps) {
        //Number of orbitals in the ordinary basis set. Also the index for the beginning of the CABS
        constexpr Index nOBS = 2;
        //Number of orbitals in the complementary auxiliary set
        constexpr Index nCABS = 10;
        //Total number of basis functions in the combined OBS and CABS basis sets
        constexpr Index nbf = nOBS + nCABS;
        constexpr Index nVirt = 1;
        //Two spin orbitals are occupied in hydrogen and one spatial orbital
        const Index nOcc = 1;

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
            {I, {0, nOcc}},
            {J, {0, nOcc}},
            {K, {0, nOcc}},
            //All virtual orbital indices
            {A, {nOcc, nVirt}},
            {B, {nOcc, nVirt}},
            //All OBS indices
            {P, {0, nOBS}},
            {Q, {0, nOBS}},
            //CABS indices
            {X, {nOBS, nCABS}},
            {Y, {nOBS, nCABS}}
        };


        //We need to assume canonical orbitals to solve equations iteratively
        Tensor2 Dai(nVirt, nOcc);
        Dai.setZero();

        for (Eigen::Index i = 0; i < nOcc; i++)
            for (Eigen::Index a = 0; a < nVirt; a++)
                Dai(a, i) = F(i, i) - F(a + nOcc, a + nOcc);


        Tensor2 Dxi(nCABS, nOcc);
        Dxi.setZero();

        for (Eigen::Index x = 0; x < nCABS; x++)
            for (Eigen::Index i = 0; i < nOcc; i++)
                Dxi(x, i) = F(i, i) - F(x + nOBS, x + nOBS);


        //Inital guess
        Tensor2 t_xi = Slice<Indices<X, I> >(idxs, F);
        Tensor2 t_ai(nVirt, nOcc);
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

            Tensor2 deltaTensor = t_xi - t1Last;
            Tensor<double, 0> deltaTensor2 = deltaTensor.sum();
            delta = abs(deltaTensor2(0));

            if (delta < 1E-8)
                break;

            t1Last = t_xi;
        }

        //Times two for pair of electrons
        //ix, xi ->
        Tensor<double, 0> ECABSSingles = 2.0 * Einsum<Indices<I, X>, Indices<X, I>, Indices<> >(
                                             Slice<Indices<I, X> >(idxs, F),
                                              t_xi);
        if (t1Amps)
            *t1Amps = t_xi;

        return ECABSSingles(0);
    }
}
