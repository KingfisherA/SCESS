#include<pch.h>
#include"IntegralTransforms.h"

#include "IntegralsF12.h"
#include "Quantum/HF/HartreeFock.h"
#include "Quantum/Basis/BasisSet.h"
#include "Utilities/Utilities.h"

namespace SCESS::Integrals {
    using namespace Eigen;

    TEITensor Transforms::TransformTEI(const TEITensor &tei, const Eigen::MatrixXd &C) {
        return TransformTEI(tei, C, C, C, C);
    }

    TEITensor Transforms::TransformTEI(const TEITensor &tei, const Eigen::MatrixXd &A, const Eigen::MatrixXd &B,
                                       const Eigen::MatrixXd &C, const Eigen::MatrixXd &D) {
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

        for (Index s = 0; s < tei.dimension(3); ++s) {
            for (Index mu = 0; mu < tei.dimension(0); ++mu)
                for (Index nu = 0; nu < tei.dimension(1); ++nu)
                    for (Index kap = 0; kap < tei.dimension(2); ++kap)
                        for (Index lam = 0; lam < tei.dimension(3); ++lam)
                            //From <mu nu | kap lam> to <mu nu | kap s>
                            //Einsum("mnkl,ls->mnks")
                            transformed1(mu, nu, kap, s) += D(lam, s) * tei(mu, nu, kap, lam);

            for (Index r = 0; r < tei.dimension(2); ++r) {
                for (Index mu = 0; mu < tei.dimension(0); ++mu)
                    for (Index nu = 0; nu < tei.dimension(1); ++nu)
                        for (Index kap = 0; kap < tei.dimension(2); ++kap)
                            //From <mu nu | kap s> to <mu nu | r s>
                            //Einsum("mnks,kr->mnrs")
                            transformed2(mu, nu, r, s) += C(kap, r) * transformed1(mu, nu, kap, s);

                for (Index q = 0; q < tei.dimension(1); ++q) {
                    for (Index mu = 0; mu < tei.dimension(0); ++mu)
                        for (Index nu = 0; nu < tei.dimension(1); ++nu)
                            //From <mu nu | r s> to <mu q | r s>
                            //Einsum("mnrs,nq->mqrs")
                            transformed3(mu, q, r, s) += B(nu, q) * transformed2(mu, nu, r, s);

                    for (Index p = 0; p < tei.dimension(0); ++p) {
                        for (Index mu = 0; mu < tei.dimension(0); ++mu)
                            //From <mu q | r s> to <p q | r s>
                            //Einsum("mqrs,mp->pqrs")
                            transformed4(p, q, r, s) += A(mu, p) * transformed3(mu, q, r, s);
                    }
                }
            }
        }

        return transformed4;
    }

    void Transforms::TransformF12Basis(const F12BasisSet &basisSet,
                                       ConfigInts &mints,
                                       F12::F12Integrals &f12Ints,
                                       uRef<UHFWavefunction> *hfWfnOut,
                                       double *cabsSinglesCorrection) {
        using namespace Eigen;

        const Index nOBS = (Index) basisSet.GetNoHFFuncs();
        const Index nCABS = (Index) basisSet.GetNoCABSFuncs();
        const Index N = (Index) basisSet.GetNumOfFuncs();

        constexpr double Vnn = 0.74298642461902142;

        assert(N == mints.S.rows());

        MatrixXd &S = mints.S;
        MatrixXd &T = mints.T;
        MatrixXd &V = mints.V;
        TEITensor &eri = mints.eri;

        uRef<UHFWavefunction> hfWfn; {
            ConfigInts hfInts;
            hfInts.S = S.block(0, 0, nOBS, nOBS);
            hfInts.T = T.block(0, 0, nOBS, nOBS);
            hfInts.V = V.block(0, 0, nOBS, nOBS);
            hfInts.eri = eri.slice(Eigen::array<Index, 4>{0, 0, 0, 0},
                                   Eigen::array<Index, 4>{nOBS, nOBS, nOBS, nOBS});

            hfWfn = HartreeFock::UHF(hfInts, 1, 1, Vnn,
                                     std::numeric_limits<double>::epsilon());
        }

        //Inverse of HF block of matrix
        MatrixXd SInv = S.block(0, 0, nOBS, nOBS).inverse();

        //Projector Matrix
        MatrixXd P = MatrixXd::Zero(N, N);

        for (Index i = 0; i < N; ++i) {
            P(i, i) = 1.0;
            for (Index mu = 0; mu < nOBS; ++mu) {
                double c = 0.0;
                for (Index nu = 0; nu < nOBS; nu++)
                    c -= SInv(mu, nu) * S(nu, i);
                P(mu, i) = c;
            }
        }

        //
        {
            MatrixXd S_tmp = P.transpose() * S * P;
            VectorXd normalisation = S_tmp.diagonal();

            //Normalise the projector
            for (Index i = 0; i < N; ++i)
                for (Index j = 0; j < N; ++j) {
                    P(i, j) /= sqrt(normalisation(j));
                }

            S_tmp = P.transpose() * S * P;

            S = S_tmp;
            T = P.transpose() * T * P;
            V = P.transpose() * V * P;

            eri = TransformTEI(eri, P);
            f12Ints.f12 = TransformTEI(f12Ints.f12, P);
            f12Ints.f12Sq = TransformTEI(f12Ints.f12Sq, P);
            f12Ints.f12DelSq = TransformTEI(f12Ints.f12DelSq, P);
            f12Ints.f12Repulsion = TransformTEI(f12Ints.f12Repulsion, P);
        }


        //Fock matrix built in CABS virtual orbitals
        MatrixXd F_CABS = MatrixXd::Zero(nCABS, nCABS);
        //HF Density matrix
        MatrixXd &D = hfWfn->DensityAlpha;
        for (Index i = nOBS; i < N; i++)
            for (Index j = nOBS; j < N; j++) {
                double g = 0.0;
                for (Index mu = 0; mu < nOBS; mu++)
                    for (Index nu = 0; nu < nOBS; nu++) {
                        double J = eri(i, j, mu, nu);
                        double K = eri(i, nu, mu, j);
                        g += D(mu, nu) * (2.0 * J - K);
                    }

                F_CABS(i - nOBS, j - nOBS) =
                        F_CABS(j - nOBS, i - nOBS) = T(i, j) + V(i, j) + g;
            }


        //Matrix to transform into canonical CABS orbitals from projected CABS orbitals
        MatrixXd C_CABS;
        //Matrix to transform from HF AOs to HF MOs
        MatrixXd &C_HF = hfWfn->CAlpha;
        //
        {
            auto S_CABS = S.block(nOBS, nOBS, nCABS, nCABS);
            SelfAdjointEigenSolver<MatrixXd> eigensolver(S_CABS);
            const VectorXd &SEigenValues = eigensolver.eigenvalues().real();
            const MatrixXd &SEigenVects = eigensolver.eigenvectors().real();
            MatrixXd SValsInvSqrt = SEigenValues.cwiseInverse().cwiseSqrt().asDiagonal();

            //Organisational matrix
            MatrixXd X = SEigenVects * SValsInvSqrt * SEigenVects.transpose();
            X = (X + X.transpose()) * 0.5;

            auto F_CABS_Prime = X.transpose() * F_CABS * X;
            SelfAdjointEigenSolver<MatrixXd> eigenSolver2(F_CABS_Prime);
            C_CABS = X * eigenSolver2.eigenvectors().real();

            if (cabsSinglesCorrection) {
                //Build combined Fock matrix
                MatrixXd F = MatrixXd::Zero(nOBS + nCABS, nOBS + nCABS);
                for (Index p = 0; p < F.rows(); p++)
                    for (Index q = p; q < F.cols(); q++) {
                        double g = 0.0;
                        for (Index mu = 0; mu < nOBS; mu++)
                            for (Index nu = 0; nu < nOBS; nu++) {
                                double J = eri(p, q, mu, nu);
                                double K = eri(p, nu, mu, q);
                                g += D(mu, nu) * (2.0 * J - K);
                            }

                        F(p, q) =
                                F(q, p) = T(p, q) + V(p, q) + g;
                    }

                GeneralizedSelfAdjointEigenSolver<MatrixXd> eigenSolver3(F, S);
                VectorXd epsilon = eigenSolver3.eigenvalues().real();
                MatrixXd F_Eigenvectors = eigenSolver3.eigenvectors().real();
                F = F_Eigenvectors.transpose() * F * F_Eigenvectors;

                (*cabsSinglesCorrection) = 0.0;
                for (Index x = nOBS; x < nCABS + nOBS; x++)
                    for (Index i = 0; i < nOBS; i++) {
                        double F_ix = F(i, x);
                        double D_xi = epsilon(x) - epsilon(i);
                        (*cabsSinglesCorrection) -= abs(F_ix * F_ix) / D_xi;
                    }
            }
        }

        array<uint8_t, 16> combinations = {0};

        for (uint8_t i = 0, x = 0; i < 2; ++i)
            for (uint8_t j = 0; j < 2; ++j)
                for (uint8_t k = 0; k < 2; ++k)
                    for (uint8_t l = 0; l < 2; ++l, ++x)
                        combinations[x] |= (i << 0) | (j << 1) | (k << 2) | (l << 3);

        auto TransformF12Tensor = [&C_CABS, &C_HF, nOBS, nCABS, &combinations](TEITensor &tei) {
            for (auto comb: combinations) {
                array<Index, 4> startIndices = {0};
                array<Index, 4> sizes = {0};

                for (uint8_t j = 0; j < 4; ++j)
                    if (comb & (0x1 << j)) {
                        startIndices[j] = nOBS;
                        sizes[j] = nCABS;
                    } else {
                        startIndices[j] = 0;
                        sizes[j] = nOBS;
                    }


                const MatrixXd &A = comb & (0x1 << 0) ? C_CABS : C_HF;
                const MatrixXd &B = comb & (0x1 << 1) ? C_CABS : C_HF;
                const MatrixXd &C = comb & (0x1 << 2) ? C_CABS : C_HF;
                const MatrixXd &D = comb & (0x1 << 3) ? C_CABS : C_HF;

                tei.slice(startIndices, sizes) = TransformTEI(tei.slice(startIndices, sizes), A, B, C, D);
            }
        };


        auto TransformF12Matrix = [&C_CABS, &C_HF, nOBS, nCABS](MatrixXd &mat) {
            MatrixXd hfBlock = mat.block(0, 0, nOBS, nOBS);
            MatrixXd cabsBlock = mat.block(nOBS, nOBS, nCABS, nCABS);
            MatrixXd hfCABSBlock = mat.block(0, nOBS, nOBS, nCABS);

            mat.block(0, 0, nOBS, nOBS) = C_HF.transpose() * hfBlock * C_HF;
            mat.block(nOBS, nOBS, nCABS, nCABS) = C_CABS.transpose() * cabsBlock * C_CABS;
            mat.block(0, nOBS, nOBS, nCABS) = C_HF.transpose() * hfCABSBlock * C_CABS;
            mat.block(nOBS, 0, nCABS, nOBS) = mat.block(0, nOBS, nOBS, nCABS).transpose();
        };

        TransformF12Matrix(mints.S);
        TransformF12Matrix(mints.T);
        TransformF12Matrix(mints.V);

        TransformF12Tensor(eri);
        f12Ints.eri = eri;

        TransformF12Tensor(f12Ints.f12);
        TransformF12Tensor(f12Ints.f12Sq);
        TransformF12Tensor(f12Ints.f12DelSq);
        TransformF12Tensor(f12Ints.f12Repulsion);

        *hfWfnOut = std::move(hfWfn);
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
