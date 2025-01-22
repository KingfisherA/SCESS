#include<Quantum/QuantumPCH.h>
#include "F12.h"

#include "Quantum/TensorOperations.h"
#include "Quantum/Integrals/Integrals.h"
#include "Quantum/Integrals/IntegralsF12.h"
#include "Quantum/Integrals/IntegralTransforms.h"
#include "Quantum/HF/HartreeFock.h"

namespace SCESS::F12 {
    using namespace Eigen;

    MP2::MP2F12Result MP2::MP2_F12(const UHFWavefunction &refWfn,
                                   const Integrals::MolInts &mints,
                                   const Integrals::F12::F12Integrals &f12Ints) {
        //We assume that the reference wavefunction is actually a closed shell RHF hence alternating alpha and beta spins (fix this in the future :( )
        // We need | Psi > = |ijk...> |alpha beta alpha beta ... > (alpha and beta are stored alternating)
        //All calculations are done in SPIN ORBITALS so we need to integrate spin out on the fly. All EVEN indices are spin alpha and ODD are spin beta

        const MatrixXd &TInts = mints.T;
        const MatrixXd &VInts = mints.V;

        // //In chemist's notation
        // TEITensor eri = Integrals::Transforms::ConvertChemistsToPhysicists(mints.eri);
        // TEITensor f12 = Integrals::Transforms::ConvertChemistsToPhysicists(f12Ints.f12);
        // TEITensor f12Sq = Integrals::Transforms::ConvertChemistsToPhysicists(f12Ints.f12Sq);
        // //F12 double commutator integrals?
        // TEITensor f12DelSq = Integrals::Transforms::ConvertChemistsToPhysicists(f12Ints.f12DelSq);
        // //Same as f12g12 integrals
        // TEITensor f12Repulsion = Integrals::Transforms::ConvertChemistsToPhysicists(f12Ints.f12Repulsion);

        const TEITensor &eri = mints.eri;
        const TEITensor &f12 = f12Ints.f12;
        const TEITensor &f12Sq = f12Ints.f12Sq;
        const TEITensor &f12DelSq = f12Ints.f12DelSq;
        const TEITensor &f12Repulsion = f12Ints.f12Repulsion;

        //nmo = Number of Molecular Orbitals
        //nalpha = Number of alpha electrons = Number of beta electrons in closed-shell RHF = ndocc
        //nfrzc = Number of frozen core electrons
        //naocc = Number of active occupied electrons
        //C = All wfn coefficients
        //Co= Occupied wfn coefficients
        //Cv = Virtual wfn coefficients
        //epsilon are the HF orbital eigenvalues

        //-- Indices --
        //i,j,k,l,o,... = Occupied orbitals
        //a,b,c,d,... = Virtual orbitals
        //p,q,r,s... = All OBS
        //x,y,w,z,... = CABS orbitals

        //!! Need to integrate out alpha and beta spin !!


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
        const Index nOcc = std::distance(refWfn.OccupationAlpha.begin(),
                                         std::ranges::find(refWfn.OccupationAlpha, 0));;

        //Make the V intermediate
        TEITensor V = f12Repulsion.Slice4({0, 0, 0, 0}, {nOcc, nOcc, nOcc, nOcc});
        //ijrs,klrs->ijkl
        V -= eri.Slice4({0, 0, 0, 0}, {nOcc, nOcc, nOBS, nOBS})
                .Contract2(
                    f12.Slice4({0,0,0,0}, {nOcc, nOcc, nOBS, nOBS}),
                    ({
                        {2,2},
                        {3,3}
                        }));
        //ijxm,klxm->ijkl
        V -= eri.Slice4({0, 0, nOBS, 0}, {nOcc, nOcc, nCABS, nOcc})
                .Contract2(
                    f12.Slice4({0, 0, nOBS, 0}, {nOcc, nOcc, nCABS, nOcc}),
                    ({
                        {2,2},
                        {3,3}
                        }));;
        //ijmx,klmx->ijkl
        V -= eri.Slice4({0, 0, 0, nOBS}, {nOcc, nOcc, nOcc, nCABS})
                .Contract2(
                    f12.Slice4({0, 0, 0, nOBS}, {nOcc, nOcc, nOcc, nCABS}),
                    ({
                        {2,2},
                        {3,3}
                        }));;


        TEITensor X = f12Sq.Slice4({0, 0, 0, 0}, {nOcc, nOcc, nOcc, nOcc});
        //ijrs,klrs->ijkl
        X -= f12.Slice4({0, 0, 0, 0}, {nOcc, nOcc, nOBS, nOBS})
                .Contract2(
                    f12.Slice4({0,0,0,0}, {nOcc, nOcc, nOBS, nOBS}),
                    ({
                        {2,2},
                        {3,3}
                        }));
        //ijxm,klxm->ijkl
        X -= f12.Slice4({0, 0, nOBS, 0}, {nOcc, nOcc, nCABS, nOcc})
                .Contract2(
                    f12.Slice4({0, 0, nOBS, 0}, {nOcc, nOcc, nCABS, nOcc}),
                    ({
                        {2,2},
                        {3,3}
                        }));;
        //ijmx,klmx->ijkl
        X -= f12.Slice4({0, 0, 0, nOBS}, {nOcc, nOcc, nOcc, nCABS})
                .Contract2(
                    f12.Slice4({0, 0, 0, nOBS}, {nOcc, nOcc, nOcc, nCABS}),
                    ({
                        {2,2},
                        {3,3}
                        }));

        //Full Fock Matrix
        Tensor<double, 2> F = MatrixToTensor(TInts + VInts);
        //Coloumb Fock Terms
        //piqi->pq
        F.Slice2({0, 0}, {nOBS, nOBS}) += 2.0 * eri.Slice4({0, 0, 0, 0}, {nOBS, nOcc, nOBS, nOcc}).Reduce2({1, 3});
        //ipix->px
        F.Slice2({0, nOBS}, {nOBS, nCABS}) += 2.0 * eri.Slice4({0, 0, 0, nOBS}, {nOcc, nOBS, nOcc, nCABS}).Reduce2({
            0, 2
        });
        //ixiy->xy
        F.Slice2({nOBS, nOBS}, {nCABS, nCABS}) += 2.0 * eri.Slice4({0, nOBS, 0, nOBS}, {nOcc, nCABS, nOcc, nCABS}).
                Reduce2({0, 2});

        Tensor<double, 2> K(F.dimensions());
        K.setZero();
        //Exchange Fock Terms
        //iipq->pq
        K.Slice2({0, 0}, {nOBS, nOBS}) = eri.Slice4({0, 0, 0, 0}, {nOcc, nOcc, nOBS, nOBS}).Reduce2({0, 1});
        //iipx->px
        K.Slice2({0, nOBS}, {nOBS, nCABS}) = eri.Slice4({0, 0, 0, nOBS}, {nOcc, nOcc, nOBS, nCABS}).Reduce2({0, 1});
        //iixy->xy
        K.Slice2({nOBS, nOBS}, {nCABS, nCABS}) = eri.Slice4({0, 0, nOBS, nOBS}, {nOcc, nOcc, nCABS, nCABS}).Reduce2({
            0, 1
        });
        F -= K;

        //Fill in the remaining elements by symmetry
        //Shuffle2(1,0)=transpose()
        F.Slice2({nOBS, 0}, {nCABS, nOBS}) = F.Slice2({0, nOBS}, {nOBS, nCABS}).Shuffle2({1, 0});
        K.Slice2({nOBS, 0}, {nCABS, nOBS}) = K.Slice2({0, nOBS}, {nOBS, nCABS}).Shuffle2({1, 0});

        //Fock matrix without the exchange terms
        Tensor<double, 2> F_K = F + K;

        //Makes the B and C intermediate tensors
        auto MakeB_C = [&f12Sq, &f12, &F, &F_K, &K, &nOcc, &nOBS, &nCABS]
        (TEITensor &B, TEITensor &C)-> void {
            //-- F12 Contributions --
            //pr,lkqr-> plkq -> klpq
            TEITensor tmp = F_K.Slice2({0, 0}, {nOBS, nOBS})
                    .Contract1(
                        f12.Slice4({0,0,0,0},{nOcc, nOcc, nOBS, nOBS}),
                        ({{1, 3}}))
                    .Shuffle4({2, 1, 0, 3});
            //mnpq,klpq-> mnkl -> klmn
            B -= f12.Slice4({0, 0, 0, 0}, {nOcc, nOcc, nOBS, nOBS})
                    .Contract2(tmp, ({{2,2},{3,3}}))
                    .Shuffle4({2, 3, 0, 1});

            //qr,klpr-> qklp -> klpq
            tmp = F_K.Slice2({0, 0}, {nOBS, nOBS})
                    .Contract1(
                        f12.Slice4({0,0,0,0},{nOcc, nOcc, nOBS, nOBS}),
                        ({{1, 3}}))
                    .Shuffle4({1, 2, 3, 0});
            //mnpq,klpq-> mnkl -> klmn
            B -= f12.Slice4({0, 0, 0, 0}, {nOcc, nOcc, nOBS, nOBS})
                    .Contract2(tmp, ({{2,2},{3,3}}))
                    .Shuffle4({2, 3, 0, 1});

            //px,lkqx-> plkq -> klpq
            tmp = F_K.Slice2({0, nOBS}, {nOBS, nCABS})
                    .Contract1(
                        f12.Slice4({0,0,0,nOBS},{nOcc, nOcc, nOBS, nCABS}),
                        ({{1, 3}}))
                    .Shuffle4({2, 1, 0, 3});
            //mnpq,klpq-> mnkl -> klmn
            B -= f12.Slice4({0, 0, 0, 0}, {nOcc, nOcc, nOBS, nOBS})
                    .Contract2(tmp, ({{2,2},{3,3}}))
                    .Shuffle4({2, 3, 0, 1});

            //qx,klpx-> qklp -> klpq
            tmp = F_K.Slice2({0, nOBS}, {nOBS, nCABS})
                    .Contract1(
                        f12.Slice4({0,0,0,nOBS},{nOcc, nOcc, nOBS, nCABS}),
                        ({{1, 3}}))
                    .Shuffle4({1, 2, 3, 0});
            //mnpq,klpq-> mnkl -> klmn
            B -= f12.Slice4({0, 0, 0, 0}, {nOcc, nOcc, nOBS, nOBS})
                    .Contract2(tmp, ({{2,2},{3,3}}))
                    .Shuffle4({2, 3, 0, 1});

            //yp,klip-> ykli ->  kliy
            tmp = F_K.Slice2({nOBS, 0}, {nCABS, nOBS})
                    .Contract1(
                        f12.Slice4({0,0,0,0},{nOcc, nOcc, nOcc, nOBS}),
                        ({{1, 3}}))
                    .Shuffle4({1, 2, 3, 0});
            //mniy,kliy-> mnkl -> klmn
            B -= f12.Slice4({0, 0, 0, nOBS}, {nOcc, nOcc, nOcc, nCABS})
                    .Contract2(tmp, ({{2,2},{3,3}}))
                    .Shuffle4({2, 3, 0, 1});

            //xp,lkjp-> xlkj -> klxj
            tmp = F_K.Slice2({nOBS, 0}, {nCABS, nOBS})
                    .Contract1(
                        f12.Slice4({0,0,0,0},{nOcc, nOcc, nOcc, nOBS}),
                        ({{1, 3}}))
                    .Shuffle4({2, 1, 0, 3});
            //nmjx,klxj-> nmkl ->  klmn
            B -= f12.Slice4({0, 0, 0, nOBS}, {nOcc, nOcc, nOcc, nCABS})
                    .Contract2(tmp, ({{2,3},{3,2}}))
                    .Shuffle4({2, 3, 1, 0});

            //ip,klpy-> ikly -> kliy
            tmp = F_K.Slice2({0, 0,}, {nOcc, nOBS})
                    .Contract1(
                        f12.Slice4({0,0,0,nOBS},{nOcc, nOcc, nOBS, nCABS}),
                        ({{1,2}}))
                    .Shuffle4({1, 2, 0, 3});
            //mniy,kliy-> mnkl -> klmn
            B -= f12.Slice4({0, 0, 0, nOBS}, {nOcc, nOcc, nOcc, nCABS})
                    .Contract2(tmp, ({{2,2},{3,3}}))
                    .Shuffle4({2, 3, 0, 1});

            //yx,klix->kliy
            tmp = F_K.Slice2({nOBS, nOBS}, {nCABS, nCABS})
                    .Contract1(
                        f12.Slice4({0, 0, 0, nOBS}, {nOcc, nOcc, nOcc, nCABS}),
                        ({{1, 3}}))
                    .Shuffle4({1, 2, 3, 0});
            //mniy,kliy->klmn
            B -= f12.Slice4({0, 0, 0, nOBS}, {nOcc, nOcc, nOcc, nCABS})
                    .Contract2(tmp,
                               ({{2, 2}, {3, 3}}))
                    .Shuffle4({3, 4, 0, 1});

            //xy,lkjy->klxj
            tmp = F_K.Slice2({nOBS, nOBS}, {nCABS, nCABS})
                    .Contract1(
                        f12.Slice4({0, 0, 0, nOBS}, {nOcc, nOcc, nOcc, nCABS}),
                        ({{1, 3}}))
                    .Shuffle4({2, 1, 0, 3});
            //nmjx,klxj->klmn
            B -= f12.Slice4({0, 0, 0, nOBS}, {nOcc, nOcc, nOcc, nCABS})
                    .Contract2(tmp,
                               ({{2, 3}, {3, 2}}))
                    .Shuffle4({3, 4, 1, 0});

            //jp,lkpx->klxj
            tmp = F_K.Slice2({0, 0}, {nOcc, nOBS})
                    .Contract1(
                        f12.Slice4({0, 0, 0, nOBS}, {nOcc, nOcc, nOBS, nCABS}),
                        ({{1, 2}}))
                    .Shuffle4({2, 1, 3, 0});
            //nmjx,klxj->klmn
            B -= f12.Slice4({0, 0, 0, nOBS}, {nOcc, nOcc, nOcc, nCABS})
                    .Contract2(tmp,
                               ({{2, 3}, {3, 2}}))
                    .Shuffle4({3, 4, 1, 0});

            //ix,klxy->kliy
            tmp = F_K.Slice2({0, nOBS}, {nOcc, nCABS})
                    .Contract1(
                        f12.Slice4({0, 0, nOBS, nOBS}, {nOcc, nOcc, nCABS, nCABS}),
                        ({{1, 2}}))
                    .Shuffle4({1, 2, 0, 3});
            //mniy,kliy->klmn
            B -= f12.Slice4({0, 0, 0, nOBS}, {nOcc, nOcc, nOcc, nCABS})
                    .Contract2(tmp,
                               ({{2, 2}, {3, 3}}))
                    .Shuffle4({3, 4, 0, 1});

            //jy,lkyx->klxj
            tmp = F_K.Slice2({0, nOBS}, {nOcc, nCABS})
                    .Contract1(
                        f12.Slice4({0, 0, nOBS, nOBS}, {nOcc, nOcc, nCABS, nCABS}),
                        ({{1, 2}}))
                    .Shuffle4({2, 1, 3, 0});
            //nmjx,klxj->klmn
            B -= f12.Slice4({0, 0, 0, nOBS}, {nOcc, nOcc, nOcc, nCABS})
                    .Contract2(tmp,
                               ({{2, 3}, {3, 2}}))
                    .Shuffle4({3, 4, 1, 0});

            //-- Y contributions --
            //xp,lkbp->klxb
            //xp,lkbp->klxb
            tmp = K.Slice2({nOBS, 0}, {nCABS, nOBS})
                    .Contract1(
                        f12.Slice4({0, 0, nOcc, 0}, {nOcc, nOcc, nVirt, nOBS}),
                        ({{1, 3}}))
                    .Shuffle4({2, 1, 0, 3});

            //klxb,nmbx->klmn
            B -= f12.Slice4({0, 0, nOcc, nOBS}, {nOcc, nOcc, nVirt, nCABS})
                    .Contract2(tmp,
                               ({{2, 3}, {3, 2}}))
                    .Shuffle4({0, 1, 4, 3});

            //yp,klap->klay
            tmp = K.Slice2({nOBS, 0}, {nCABS, nOBS})
                    .Contract1(
                        f12.Slice4({0, 0, nOcc, 0}, {nOcc, nOcc, nVirt, nOBS}),
                        ({{1, 3}}))
                    .Shuffle4({1, 2, 3, 0});

            //klay,mnay->klmn
            B -= f12.Slice4({0, 0, nOcc, nOBS}, {nOcc, nOcc, nVirt, nCABS})
                    .Contract2(tmp,
                               ({{2, 2}, {3, 3}}))
                    .Shuffle4({0, 1, 3, 4});

            //xy,lkby->klxb
            tmp = K.Slice2({nOBS, nOBS}, {nCABS, nCABS})
                    .Contract1(
                        f12.Slice4({0, 0, nOcc, nOBS}, {nOcc, nOcc, nVirt, nCABS}),
                        ({{1, 3}}))
                    .Shuffle4({2, 1, 0, 3});

            //klxb,nmbx->klmn
            B -= f12.Slice4({0, 0, nOcc, nOBS}, {nOcc, nOcc, nVirt, nCABS})
                    .Contract2(tmp,
                               ({{2, 3}, {3, 2}}))
                    .Shuffle4({0, 1, 4, 3});

            //bp,lkpx->klxb
            tmp = K.Slice2({nOcc, 0}, {nVirt, nOBS})
                    .Contract1(
                        f12.Slice4({0, 0, 0, nOBS}, {nOcc, nOcc, nOBS, nCABS}),
                        ({{1, 2}}))
                    .Shuffle4({2, 1, 3, 0});

            //klxb,nmbx->klmn
            B -= f12.Slice4({0, 0, nOcc, nOBS}, {nOcc, nOcc, nVirt, nCABS})
                    .Contract2(tmp,
                               ({{2, 3}, {3, 2}}))
                    .Shuffle4({0, 1, 4, 3});

            //ap,klpy->klay
            tmp = K.Slice2({nOcc, 0}, {nVirt, nOBS})
                    .Contract1(
                        f12.Slice4({0, 0, 0, nOBS}, {nOcc, nOcc, nOBS, nCABS}),
                        ({{1, 2}}))
                    .Shuffle4({1, 2, 0, 3});

            //klay,mnay->klmn
            B -= f12.Slice4({0, 0, nOcc, nOBS}, {nOcc, nOcc, nVirt, nCABS})
                    .Contract2(tmp,
                               ({{2, 2}, {3, 3}}))
                    .Shuffle4({0, 1, 3, 4});

            //yx,klax->klay
            tmp = K.Slice2({nOBS, nOBS}, {nCABS, nCABS})
                    .Contract1(
                        f12.Slice4({0, 0, nOcc, nOBS}, {nOcc, nOcc, nVirt, nCABS}),
                        ({{1, 3}}))
                    .Shuffle4({1, 2, 3, 0});

            //klay,mnay->klmn
            B -= f12.Slice4({0, 0, nOcc, nOBS}, {nOcc, nOcc, nVirt, nCABS})
                    .Contract2(tmp,
                               ({{2, 2}, {3, 3}}))
                    .Shuffle4({0, 1, 3, 4});

            //by,lkyx->klxb
            tmp = K.Slice2({nOcc, nOBS}, {nVirt, nCABS})
                    .Contract1(
                        f12.Slice4({0, 0, nOBS, nOBS}, {nOcc, nOcc, nCABS, nCABS}),
                        ({{1, 2}}))
                    .Shuffle4({2, 1, 3, 0});

            //klxb,nmbx->klmn
            B -= f12.Slice4({0, 0, nOcc, nOBS}, {nOcc, nOcc, nVirt, nCABS})
                    .Contract2(tmp,
                               ({{2, 3}, {3, 2}}))
                    .Shuffle4({0, 1, 4, 3});

            //ax,klxy->klay
            tmp = K.Slice2({nOcc, nOBS}, {nVirt, nCABS})
                    .Contract1(
                        f12.Slice4({0, 0, nOBS, nOBS}, {nOcc, nOcc, nCABS, nCABS}),
                        ({{1, 2}}))
                    .Shuffle4({1, 2, 0, 3});

            //klay,mnay->klmn
            B -= f12.Slice4({0, 0, nOcc, nOBS}, {nOcc, nOcc, nVirt, nCABS})
                    .Contract2(tmp,
                               ({{2, 2}, {3, 3}}))
                    .Shuffle4({0, 1, 3, 4});

            //xp,klpy->klxy
            tmp = K.Slice2({nOBS, 0}, {nCABS, nOBS})
                    .Contract1(
                        f12.Slice4({0, 0, 0, nOBS}, {nOcc, nOcc, nOBS, nCABS}),
                        ({{1, 2}}))
                    .Shuffle4({1, 2, 0, 3});

            //klxy,mnxy->klmn
            B -= f12.Slice4({0, 0, nOBS, nOBS}, {nOcc, nOcc, nCABS, nCABS})
                    .Contract2(tmp,
                               ({{2, 2}, {3, 3}}))
                    .Shuffle4({0, 1, 3, 4});

            //yp,lkpx->klxy
            tmp = K.Slice2({nOBS, 0}, {nCABS, nOBS})
                    .Contract1(
                        f12.Slice4({0, 0, 0, nOBS}, {nOcc, nOcc, nOBS, nCABS}),
                        ({{1, 2}}))
                    .Shuffle4({2, 1, 3, 0});

            //klxy,mnxy->klmn
            B -= f12.Slice4({0, 0, nOBS, nOBS}, {nOcc, nOcc, nCABS, nCABS})
                    .Contract2(tmp,
                               ({{2, 2}, {3, 3}}))
                    .Shuffle4({0, 1, 3, 4});

            //yz,lkzx->klxy
            tmp = K.Slice2({nOBS, nOBS}, {nCABS, nCABS})
                    .Contract1(
                        f12.Slice4({0, 0, nOBS, nOBS}, {nOcc, nOcc, nCABS, nCABS}),
                        ({{1, 2}}))
                    .Shuffle4({2, 1, 3, 0});

            //klxy,mnxy->klmn
            B -= f12.Slice4({0, 0, nOBS, nOBS}, {nOcc, nOcc, nCABS, nCABS})
                    .Contract2(tmp,
                               ({{2, 2}, {3, 3}}))
                    .Shuffle4({0, 1, 3, 4});

            //xz,klzy->klxy
            tmp = K.Slice2({nOBS, nOBS}, {nCABS, nCABS})
                    .Contract1(
                        f12.Slice4({0, 0, nOBS, nOBS}, {nOcc, nOcc, nCABS, nCABS}),
                        ({{1, 2}}))
                    .Shuffle4({1, 2, 0, 3});

            //klxy,mnxy->klmn
            B -= f12.Slice4({0, 0, nOBS, nOBS}, {nOcc, nOcc, nCABS, nCABS})
                    .Contract2(tmp,
                               ({{2, 2}, {3, 3}}))
                    .Shuffle4({0, 1, 3, 4});


            //--F Squared Contributions--
            //kp,nmlp->klmn
            B += F_K.Slice2({0, 0}, {nOcc, nOBS})
                    .Contract1(
                        f12Sq.Slice4({0, 0, 0, 0}, {nOcc, nOcc, nOcc, nOBS}),
                        ({{1, 3}}))
                    .Shuffle4({0, 3, 2, 1});
            //lp,mnkp->klmn
            B += F_K.Slice2({0, 0}, {nOcc, nOBS})
                    .Contract1(
                        f12Sq.Slice4({0, 0, 0, 0}, {nOcc, nOcc, nOcc, nOBS}),
                        ({{1, 3}}))
                    .Shuffle4({3, 0, 1, 2});
            //kx,nmlx->klmn
            B += F_K.Slice2({0, nOBS}, {nOcc, nCABS})
                    .Contract1(
                        f12Sq.Slice4({0, 0, 0, nOBS}, {nOcc, nOcc, nOcc, nCABS}),
                        ({{1, 3}}))
                    .Shuffle4({0, 3, 2, 1});
            //lx,mnkx->klmn
            B += F_K.Slice2({0, nOBS}, {nOcc, nCABS})
                    .Contract1(
                        f12Sq.Slice4({0, 0, 0, nOBS}, {nOcc, nOcc, nOcc, nCABS}),
                        ({{1, 3}}))
                    .Shuffle4({3, 0, 1, 2});

            //-- C Tensor --
            //bx,klax->klab
            C += F.Slice2({nOcc, nOBS}, {nVirt, nCABS})
                    .Contract1(
                        f12.Slice4({0, 0, nOcc, nOBS}, {nOcc, nOcc, nVirt, nCABS}),
                        ({{1, 3}}))
                    .Shuffle4({1, 2, 3, 0});
            //ax,lkbx->klab
            C += F.Slice2({nOcc, nOBS}, {nVirt, nCABS})
                    .Contract1(
                        f12.Slice4({0, 0, nOcc, nOBS}, {nOcc, nOcc, nVirt, nCABS}),
                        ({{1, 3}}))
                    .Shuffle4({2, 1, 0, 3});


            //-- -FC contribution to B--
            //klab,mnab->klmn
            B -= f12.Slice4({0, 0, nOcc, nOcc}, {nOcc, nOcc, nVirt, nVirt})
                    .Contract2(C,
                               ({{2, 2}, {3, 3}}))
                    .Shuffle4({0, 1, 3, 4});

            //--Symmetrise B--
            B = 0.5 * (B + B.Shuffle4({2, 3, 0, 1}));
        };

        TEITensor B = f12DelSq.Slice4({0, 0, 0, 0}, {nOcc, nOcc, nOcc, nOcc});
        TEITensor C(nOcc, nOcc, nVirt, nVirt);
        C.setZero();
        MakeB_C(B, C);

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
        auto MakeVTildeSpin = [&eri, &nOcc,&nVirt, &V, &C, &F](MatrixXd &V_AlphaAlpha,
                                                               MatrixXd &V_AlphaBeta,
                                                               Index i, Index j) {
            TEITensor eri_oovv = eri.Slice4({0, 0, nOcc, nOcc},
                                            {nOcc, nOcc, nVirt, nVirt});
            Index col = 0;
            for (Index k = 0; k < nOcc; k++)
                for (Index l = 0; l < nOcc; l++, col++) {
                    double VAlphaBeta_kl = V(k, l, i, j);
                    double VAlphaAlpha_kl = V(k, l, i, j) - V(l, k, i, j);
                    for (Index a = 0; a < nVirt; a++)
                        for (Index b = 0; b < nVirt; b++) {
                            double D_ijab = F(nOcc + a, nOcc + a) + F(nOcc + b, nOcc + b)
                                            - F(i, i) - F(j, j);
                            double C_klab = C(k, l, a, b);
                            double C_klba = C(k, l, b, a);
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
        auto MakeBTildeSpin = [&F, &C, &X, &B, &nOcc](MatrixXd &B_AlphaAlpha,
                                                      MatrixXd &B_AlphaBeta,
                                                      Index i, Index j) {
            double epsilon_ij = F(i, i) + F(j, j);
            Index row = 0, col = 0;
            for (Index k = 0; k < nOcc; k++)
                for (Index l = 0; l < nOcc; l++, row++) {
                    for (Index m = 0; m < nOcc; m++)
                        for (Index n = 0; n < nOcc; n++, col++) {
                            double B_AlphaBeta_klmn = B(k, l, m, n)
                                                      - epsilon_ij * X(k, l, m, n);
                            double B_AlphaAlpha_klmn = 0.0;

                            if (k < l && m < n)
                                B_AlphaAlpha_klmn =
                                        B(k, l, m, n)
                                        - B(k, l, n, m)
                                        - epsilon_ij * X(k, l, m, n)
                                        + epsilon_ij * X(k, l, n, m);

                            for (Index a = 0; a < nVirt; a++)
                                for (Index b = 0; b < nVirt; b++) {
                                    double D_ijab = F(nOcc + a, nOcc + a) + F(nOcc + b, nOcc + b) - F(
                                                        i, i) - F(j, j);
                                    double C_klab = C(k, l, a, b);
                                    double C_klba = C(k, l, b, a);
                                    double C_mnab = C(m, n, a, b);
                                    double C_mnba = C(m, n, b, a);

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
