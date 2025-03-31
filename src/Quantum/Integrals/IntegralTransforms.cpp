#include<Quantum/Q_PCH.h>
#include"IntegralTransforms.h"
//
#include "IntegralsF12.h"
#include "Tensor/TensorOperations.h"
#include "Quantum/HF/HF.h"
#include "Quantum/Basis/BasisSet.h"

//
//
namespace SCESS::Integrals::Transforms {
    //
    //
    using namespace Eigen;
    //
    TEITensor TransformTEI(const TEITensor &tei, const Eigen::MatrixXd &C, bool convertChemToPhys) {
        return TransformTEI(tei, C, C, C, C, convertChemToPhys);
    }
    //
    TEITensor TransformTEI(const TEITensor &tei, const Eigen::MatrixXd &A, const Eigen::MatrixXd &B,
                           const Eigen::MatrixXd &C, const Eigen::MatrixXd &D, bool convertChemToPhys) {
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

        if (convertChemToPhys) {
            transformed3 = transformed4;
            //mu, nu, kap, lam -> mu, kap, nu, lam
            transformed4 = transformed3.Shuffle4({0, 2, 1, 3});
        }

        return transformed4;
    }
    //
    void TransformF12Basis(const F12BasisSet &basisSet,
                           const Eigen::MatrixXd &HF_Density,
                           const Eigen::MatrixXd &C_HF,
                           MolInts &mints,
                           F12::F12Integrals &f12Ints,
                           bool convertChemToPhys,
                           Eigen::MatrixXd *fockOut) {
        using namespace Eigen;

        const auto nOBS = (Index) basisSet.GetNoHFFuncs();
        const auto nCABS = (Index) basisSet.GetNoCABSFuncs();
        const auto N = (Index) basisSet.GetNumOfFuncs();

        assert(N == mints.S.rows());

        MatrixXd &S = mints.S;
        MatrixXd &T = mints.T;
        MatrixXd &V = mints.V;
        TEITensor &eri = mints.eri;

        //Inverse of HF block of matrix
        MatrixXd SInv = S.block(0, 0, nOBS, nOBS).inverse();

        //Projector Matrix
        MatrixXd P = MatrixXd::Zero(N, N);

        for (Index i = 0; i < N; ++i) {
            if (i >= nOBS)
                P(i, i) = 1.0;
            for (Index mu = 0; mu < nOBS; ++mu) {
                double c = 0.0;
                for (Index nu = 0; nu < nOBS; nu++)
                    c -= SInv(mu, nu) * S(nu, i);

                if (i != mu)
                    P(mu, i) = c;
            }
        }

        //Canonical orthogonalisation with normalisation
        auto S_tmp = P.transpose() * S * P;
        SelfAdjointEigenSolver<MatrixXd> orthogonalEig(S_tmp);
        MatrixXd U = orthogonalEig.eigenvectors().real();
        VectorXd s = orthogonalEig.eigenvalues().real();
        MatrixXd sInvSqrt = MatrixXd::Zero(s.rows(), s.rows());

        for (Index i = 0; i < s.rows(); ++i)
            if (abs(s(i)) > 1E-8)
                sInvSqrt(i, i) = 1.0 / std::sqrt(s(i));

        MatrixXd X = U * sInvSqrt;
        MatrixXd P_Prime = P * X;

        //Fock matrix built in CABS virtual orbitals
        MatrixXd F; {
            MatrixXd G = MatrixXd::Zero(N, N);
            //HF Density matrix
            const MatrixXd &D = HF_Density;
            for (Index i = 0; i < N; i++)
                for (Index j = 0; j < N; j++)
                    for (Index mu = 0; mu < nOBS; mu++)
                        for (Index nu = 0; nu < nOBS; nu++) {
                            double J = eri(i, j, mu, nu);
                            double K = eri(i, nu, mu, j);

                            G(i, j) += D(mu, nu) * (2.0 * J - K);
                        }

            F = T + V + G;
        }

        auto F_Prime = P_Prime.transpose() * F * P_Prime;
        SelfAdjointEigenSolver<MatrixXd> fockEig(F_Prime);

        //Matrix which transforms both the OBS and CABS into the projected canonicalised basis set
        //This matrix is different to MRCC's as the degenerate eigenvectors are different
        //but any observable property will be the same
        MatrixXd Q = P_Prime * fockEig.eigenvectors().real();
        Q.block(0, 0, nOBS, nOBS) = C_HF;

        auto TransformF12Tensor = [&Q, convertChemToPhys](TEITensor &tei) {
            tei = TransformTEI(tei, Q, convertChemToPhys);
        };

        auto TransformF12Matrix = [&Q](MatrixXd &m) {
            m = Q.transpose() * m * Q;
        };

        //Transform the configuration integrals
        TransformF12Matrix(F);
        TransformF12Matrix(S);
        TransformF12Matrix(T);
        TransformF12Matrix(V);
        TransformF12Tensor(eri);

        //Transform the F12 integrals
        TransformF12Tensor(f12Ints.f12);
        TransformF12Tensor(f12Ints.f12Sq);
        TransformF12Tensor(f12Ints.f12DelSq);
        TransformF12Tensor(f12Ints.f12Repulsion);

        if (fockOut)
            *fockOut = std::move(F);
    }
    //
    TEITensor AntisymmetriseTEI(const TEITensor &tei) {
        TEITensor transformed(tei.dimensions());

        for (Index p = 0; p < transformed.dimension(0); ++p)
            for (Index q = 0; q < transformed.dimension(1); ++q)
                for (Index r = 0; r < transformed.dimension(2); ++r)
                    for (Index s = 0; s < transformed.dimension(3); ++s)
                        transformed(p, q, r, s) = tei(p, r, q, s) - tei(p, s, q, r);

        return transformed;
    }
    //
    TEITensor AntisymmetriseIntegrateOutSpinClosedShell(const TEITensor &moTEI) {
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
    //
    TEITensor ConvertChemistsToPhysicists(const TEITensor &tei) {
        //mu, nu, kap, lam -> mu, kap, nu, lam
        return tei.Shuffle4({0, 2, 1, 3});

        // TEITensor transformed(tei.dimensions());
        // Index dim = tei.dimension(0);
        //
        // for (Index mu = 0; mu < dim; ++mu)
        //     for (Index nu = 0; nu < dim; ++nu)
        //         for (Index kap = 0; kap < dim; ++kap)
        //             for (Index lam = 0; lam < dim; ++lam)
        //                 transformed(mu, nu, kap, lam) = tei(mu, kap, nu, lam);

        //return transformed;
    }
    //
    TEITensor ConvertPhysicistsToChemists(const TEITensor &tei) {
        //mu, nu, kap, lam -> mu, kap, nu, lam
        return tei.Shuffle4({0, 2, 1, 3});
    }
    //
    //
}
