#include<pch.h>
#include"IntegralTransforms.h"

#include "IntegralsF12.h"
#include "Quantum/HF/HartreeFock.h"
#include "Quantum/Basis/BasisUtils.h"
#include "Utilities/Utilities.h"

namespace SCESS::Integrals {
    using namespace Eigen;

    TEITensor Transforms::TransformTEI(const TEITensor &tei, const Eigen::MatrixXd &C) {
        //Note that the tei tensors are in chemists’ notation

        //tei(p,q,r,s)= C(mu,p)*C(nu,q)*tei(mu, nu, kap, lam)*C(kap,r)*C(lam, s)

        TEITensor transformed1(tei.dimensions());
        transformed1.setZero();
        TEITensor transformed2(tei.dimensions());
        transformed2.setZero();
        TEITensor transformed3(tei.dimensions());
        transformed3.setZero();
        TEITensor transformed4(tei.dimensions());
        transformed4.setZero();

        Index dim = tei.dimension(0);

        for (Index s = 0; s < dim; ++s) {
            for (Index mu = 0; mu < dim; ++mu)
                for (Index nu = 0; nu < dim; ++nu)
                    for (Index kap = 0; kap < dim; ++kap)
                        for (Index lam = 0; lam < dim; ++lam)
                            //From <mu nu | kap lam> to <mu nu | kap s>
                            //Einsum("mnkl,ls->mnks")
                            transformed1(mu, nu, kap, s) += C(lam, s) * tei(mu, nu, kap, lam);

            for (Index r = 0; r < dim; ++r) {
                for (Index mu = 0; mu < dim; ++mu)
                    for (Index nu = 0; nu < dim; ++nu)
                        for (Index kap = 0; kap < dim; ++kap)
                            //From <mu nu | kap s> to <mu nu | r s>
                            //Einsum("mnks,kr->mnrs")
                            transformed2(mu, nu, r, s) += C(kap, r) * transformed1(mu, nu, kap, s);

                for (Index q = 0; q < dim; ++q) {
                    for (Index mu = 0; mu < dim; ++mu)
                        for (Index nu = 0; nu < dim; ++nu)
                            //From <mu nu | r s> to <mu q | r s>
                            //Einsum("mnrs,nq->mqrs")
                            transformed3(mu, q, r, s) += C(nu, q) * transformed2(mu, nu, r, s);

                    for (Index p = 0; p < dim; ++p) {
                        for (Index mu = 0; mu < dim; ++mu)
                            //From <mu q | r s> to <p q | r s>
                            //Einsum("mqrs,mp->pqrs")
                            transformed4(p, q, r, s) += C(mu, p) * transformed3(mu, q, r, s);
                    }
                }
            }
        }

        return transformed4;
    }

    void Transforms::TransformF12Basis(const BasisUtils::F12BasisSet &basisSet,
                                       ConfigIntegrals &mints,
                                       F12::F12Integrals &f12Ints,
                                       uRef<UHFWavefunction> *hfWfnOut,
                                       uRef<UHFWavefunction> *cabsHFWfnOut) {
        using namespace Eigen;

        const Index nHF = basisSet.GetNoHFFuncs();
        const Index nCABS = basisSet.GetNoCABSFuncs();
        const Index N = basisSet.GetNumOfFuncs();

        constexpr double Vnn = 0.74298642461902142;

        assert(N == mints.S.rows());

        MatrixXd &S = mints.S;
        MatrixXd &T = mints.T;
        MatrixXd &V = mints.V;
        TEITensor &eri = mints.eri;

        uRef<UHFWavefunction> hfWfn; {
            ConfigIntegrals hfInts;
            hfInts.S = S.block(0, 0, nHF, nHF);
            hfInts.T = T.block(0, 0, nHF, nHF);
            hfInts.V = V.block(0, 0, nHF, nHF);
            hfInts.eri = eri.slice(Eigen::array<Index, 4>{0, 0, 0, 0},
                                   Eigen::array<Index, 4>{nHF, nHF, nHF, nHF});

            hfWfn = HartreeFock::UHF(hfInts, 1, 1, Vnn,
                                     std::numeric_limits<double>::epsilon());
        }

        //Orthoganalise the OBS AOs
        // {
        //     SelfAdjointEigenSolver<MatrixXd> eigensolver(S.block(0, 0, nHF, nHF));
        //     const VectorXd &SEigenValues = eigensolver.eigenvalues().real();
        //     const MatrixXd &SEigenVects = eigensolver.eigenvectors().real();
        //     MatrixXd SValsInvSqrt = SEigenValues.cwiseInverse().cwiseSqrt().asDiagonal();
        //     MatrixXd X = SEigenVects * SValsInvSqrt * SEigenVects.transpose();
        //     X = (X + X.transpose()) * 0.5;
        //
        //     auto TransformMatrix = [nHF, nCABS](MatrixXd &M, const MatrixXd& X) {
        //         MatrixXd tmp = X.transpose() * M.block(0, 0, nHF, nHF) * X;
        //         M.block(0, 0, nHF, nHF) = tmp;
        //
        //         tmp = M.block(nHF, 0, nCABS, nHF) * X;
        //         M.block(nHF, 0, nCABS, nHF) = tmp;
        //         M.block(0, nHF, nHF, nCABS) = tmp.transpose();
        //     };
        //
        //     // auto TransformTensor = [nHF, nCABS](TEITensor &tei, const MatrixXd& X) {
        //     //     tei.slice(Eigen::array<Index, 4>{0, 0, 0, 0},
        //     //               Eigen::array<Index, 4>{nHF, nHF, nHF, nHF}) = TransformTEI(tei, X);
        //     // };
        //
        //
        //     TransformMatrix(S, X);
        //     TransformMatrix(T, X);
        //     TransformMatrix(V, X);
        //     //TransformTensor(eri, X);
        //
        //
        //     // for (Index i = 0; i < N; i++)
        //     //     for (Index j = i; j < N; j++) {
        //     //         if (i + 1 > nHF && j + 1 > nHF) {
        //     //             S_tmp(i, j) = S_tmp(j, i) = S(i, j);
        //     //             T_tmp(i, j) = T_tmp(j, i) = T(i, j);
        //     //             V_tmp(i, j) = V_tmp(j, i) = V(i, j);
        //     //             continue;
        //     //         }
        //     //
        //     //         if (i < nHF && j < nHF) {
        //     //             for (Index mu = 0; mu < X.rows(); mu++)
        //     //                 for (Index nu = 0; nu < X.cols(); nu++) {
        //     //                     S_tmp(i, j) = S_tmp(j, i)
        //     //                                   += X(mu, i) * X(nu, j) * S(mu, nu);
        //     //                     T_tmp(i, j) = T_tmp(j, i)
        //     //                                   += X(mu, i) * X(nu, j) * T(mu, nu);
        //     //                     V_tmp(i, j) = V_tmp(j, i)
        //     //                                   += X(mu, i) * X(nu, j) * V(mu, nu);
        //     //                 }
        //     //             continue;
        //     //         }
        //     //
        //     //
        //     //         //k is the HF index, l is the CABS index always
        //     //         Index k = i, l = j;
        //     //         if (l < nHF)
        //     //             std::swap(k, l);
        //     //
        //     //         for (Index mu = 0; mu < X.rows(); mu++) {
        //     //             S_tmp(i, j) = S_tmp(j, i)
        //     //                           += X(mu, k) * S(mu, l);
        //     //             T_tmp(i, j) = T_tmp(j, i)
        //     //                           += X(mu, k) * T(mu, l);
        //     //             V_tmp(i, j) = V_tmp(j, i)
        //     //                           += X(mu, k) * V(mu, l);
        //     //         }
        //     //     }
        // }

        *hfWfnOut = std::move(hfWfn);

        std::cout << "S\n";
        Utilities::PrintMatrix(S);

        //Inverse of HF block of matrix
        MatrixXd SInv = S.block(0, 0, nHF, nHF).inverse();

        //Projection Matrix
        MatrixXd P = MatrixXd::Zero(N, N);

        for (Index i = 0; i < N; ++i) {
            P(i, i) = 1.0;
            for (Index j = 0; j < nHF; ++j) {
                double c = 0.0;
                for (Index mu = 0; mu < nHF; mu++)
                    c -= SInv(j, mu) * S(mu, i);
                P(j, i) = c;
            }
        }


        //
        {
            MatrixXd S_tmp = P.transpose() * S * P;
            MatrixXd T_tmp = P.transpose() * T * P;
            MatrixXd V_tmp = P.transpose() * V * P;
            S = S_tmp;
            T = T_tmp;
            V = V_tmp;

            std::cout << "STilde\n";
            Utilities::PrintMatrix(S);


            eri = TransformTEI(eri, P);
            f12Ints.f12 = TransformTEI(f12Ints.f12, P);
            f12Ints.f12Sq = TransformTEI(f12Ints.f12Sq, P);
            f12Ints.f12DelSq = TransformTEI(f12Ints.f12DelSq, P);
            f12Ints.f12Repulsion = TransformTEI(f12Ints.f12Repulsion, P);
        }

        ConfigIntegrals cabsConfigInts;
        cabsConfigInts.S = S.block(nHF, nHF, nCABS, nCABS);
        cabsConfigInts.T = T.block(nHF, nHF, nCABS, nCABS);
        cabsConfigInts.V = V.block(nHF, nHF, nCABS, nCABS);
        cabsConfigInts.eri = eri.slice(Eigen::array<Index, 4>{nHF, nHF, nHF, nHF},
                                       Eigen::array<Index, 4>{nCABS, nCABS, nCABS, nCABS});

        uRef<UHFWavefunction> cabsHFWfn = HartreeFock::UHF(cabsConfigInts, 1, 1, Vnn,
                                                           std::numeric_limits<double>::epsilon());


        MatrixXd SCabsTilde = cabsHFWfn->CAlpha.transpose() * cabsHFWfn->FockAlpha * cabsHFWfn->CAlpha;

        std::cout << "CABS F\n";
        Utilities::PrintMatrix(SCabsTilde);

        const MatrixXd &C = cabsHFWfn->CAlpha;
        f12Ints.eri = ConvertChemistsToPhysicists(TransformTEI(
            cabsConfigInts.eri, C));
        f12Ints.f12 = ConvertChemistsToPhysicists(TransformTEI(
            f12Ints.f12.slice(
                Eigen::array<Index, 4>{nHF, nHF, nHF, nHF},
                Eigen::array<Index, 4>{
                    nCABS, nCABS, nCABS, nCABS
                }), C));
        f12Ints.f12Sq = ConvertChemistsToPhysicists(TransformTEI(
            f12Ints.f12Sq.slice(
                Eigen::array<Index, 4>{nHF, nHF, nHF, nHF},
                Eigen::array<Index, 4>{
                    nCABS, nCABS, nCABS, nCABS
                }), C));
        f12Ints.f12DelSq = ConvertChemistsToPhysicists(TransformTEI(
            f12Ints.f12DelSq.slice(
                Eigen::array<Index, 4>{nHF, nHF, nHF, nHF},
                Eigen::array<Index, 4>{
                    nCABS, nCABS, nCABS, nCABS
                }), C));
        f12Ints.f12Repulsion = ConvertChemistsToPhysicists(TransformTEI(
            f12Ints.f12Repulsion.slice(Eigen::array<Index, 4>{nHF, nHF, nHF, nHF},
                                       Eigen::array<Index, 4>{
                                           nCABS, nCABS, nCABS, nCABS
                                       }), C));


        *cabsHFWfnOut = std::move(cabsHFWfn);
    }

    TEITensor Transforms::AntisymmetriseTEI(const TEITensor &tei) {
        TEITensor transformed(tei.dimensions());

        for (Index p = 0; p < transformed.dimension(0); ++p)
            for (Index q = 0; q < transformed.dimension(1); ++q)
                for (Index r = 0; r < transformed.dimension(2); ++r)
                    for (Index s = 0; s < transformed.dimension(3); ++s)
                        transformed(p, q, r, s) = tei(p, r, q, s) - tei(p, s, q, r);

        return transformed;
    }

    TEITensor Transforms::AntisymmetriseIntegrateOutSpinClosedShell(const TEITensor &moTEI) {
        Eigen::DSizes<Index, 4> dimensions = moTEI.dimensions();

        //Make the tensor twice as big to hold the spin-integrated integrals
        //Each integral in the moTEI actually represents two electrons
        for (Index &dim: dimensions)
            dim *= 2;

        TEITensor transformed(dimensions);

        for (Index p = 0; p < transformed.dimension(0); ++p)
            for (Index q = 0; q < transformed.dimension(1); ++q)
            for (Index q = 0; q < transformed.dimension(1); ++q)
                for (Index r = 0; r < transformed.dimension(2); ++r)
                    for (Index s = 0; s < transformed.dimension(3); ++s) {
                        //Only allow the integral to exist if the alpha is integrated with an alpha
                        //or beta with a beta (spin function orthogonality)
                        Index coloumbMask = (p % 2 == r % 2) * (q % 2 == s % 2);
                        //One if the spin functions commute. Zero otherwise
                        Index exchangeMask = (p % 2 == s % 2) * (q % 2 == r % 2);

                        //Convert from chemist's to physicists' notation and apply spin integration
                        transformed(p, q, r, s) =
                                moTEI(p / 2, r / 2, q / 2, s / 2) * static_cast<double>(coloumbMask)
                                - moTEI(p / 2, s / 2, q / 2, r / 2) * static_cast<double>(exchangeMask);
                    }

        return transformed;
    }

    TEITensor Transforms::ConvertChemistsToPhysicists(const TEITensor &tei) {
        TEITensor transformed(tei.dimensions());
        Index dim = tei.dimension(0);

        for (Index mu = 0; mu < dim; ++mu)
            for (Index nu = 0; nu < dim; ++nu)
                for (Index kap = 0; kap < dim; ++kap)
                    for (Index lam = 0; lam < dim; ++lam)
                        transformed(mu, nu, kap, lam) = tei(mu, kap, nu, lam);

        return transformed;
    }
}
