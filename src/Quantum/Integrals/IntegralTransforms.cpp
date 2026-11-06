#include<Quantum/Q_PCH.h>
#include"IntegralTransforms.h"
//
#include "Tensor/TensorEinsum.h"
#include "Tensor/TensorOperations.h"

//
//
namespace SCESS::Integrals::Transforms {
    //
    //
    using namespace Eigen;
    //
    MatrixXd CoulombMatrix(const TEITensor &eri, const MatrixXd &D) {
        const Index n = D.rows();
        //The column-major (mu nu | lam sig) is the n^2 x n^2 matrix M(mu + n nu, lam + n sig), so vec(J) = M vec(D)
        const Map<const MatrixXd> M(eri.data(), n * n, n * n);
        const Map<const VectorXd> vecD(D.data(), n * n);
        VectorXd vecJ = M * vecD;
        return Map<MatrixXd>(vecJ.data(), n, n);
    }
    //
    MatrixXd ExchangeMatrix(const TEITensor &eri, const MatrixXd &D) {
        enum { Mu, Nu, Lam, Sig };
        return TensorToMatrix(Evaluate(Einsum<Indices<Mu, Lam, Nu, Sig>, Indices<Lam, Sig>, Indices<Mu, Nu> >(
            eri, MatrixToTensor(D))));
    }
    //
    TEITensor TransformTEI(const TEITensor &tei, const Eigen::MatrixXd &C, bool convertChemToPhys) {
        return TransformTEI(tei, C, C, C, C, convertChemToPhys);
    }
    //
    TEITensor TransformTEI(const TEITensor &tei, const Eigen::MatrixXd &A, const Eigen::MatrixXd &B,
                           const Eigen::MatrixXd &C, const Eigen::MatrixXd &D, bool convertChemToPhys) {
        //Note that the tei tensors are in chemists’ notation
        //tei(p,q,r,s)= C(mu,p)*C(nu,q)*tei(mu, nu, kap, lam)*C(kap,r)*C(lam, s), done one index at a time
        enum { Mu, Nu, Kap, Lam, P, Q, R, S };

        //Each quarter transform is evaluated into a new tensor before the previous one is released, so at most two
        //four-index tensors are alive at once
        //From (mu nu | kap lam) to (mu nu | kap s)
        TEITensor transformed = Evaluate(Einsum<Indices<Mu, Nu, Kap, Lam>, Indices<Lam, S>, Indices<Mu, Nu, Kap, S> >(
            tei, MatrixToTensor(D)));
        //From (mu nu | kap s) to (mu nu | r s)
        transformed = Evaluate(Einsum<Indices<Mu, Nu, Kap, S>, Indices<Kap, R>, Indices<Mu, Nu, R, S> >(
            transformed, MatrixToTensor(C)));
        //From (mu nu | r s) to (mu q | r s)
        transformed = Evaluate(Einsum<Indices<Mu, Nu, R, S>, Indices<Nu, Q>, Indices<Mu, Q, R, S> >(
            transformed, MatrixToTensor(B)));
        //From (mu q | r s) to (p q | r s)
        transformed = Evaluate(Einsum<Indices<Mu, Q, R, S>, Indices<Mu, P>, Indices<P, Q, R, S> >(
            transformed, MatrixToTensor(A)));

        if (convertChemToPhys)
            //p q r s -> p r q s
            transformed = Evaluate(Einsum<Indices<P, Q, R, S>, Indices<P, R, Q, S> >(transformed));

        return transformed;
    }
    //
    MatrixXd BuildCABSCoefficients(const MatrixXd &S,
                                   const MatrixXd &F_AO,
                                   const MatrixXd &C_HF,
                                   const Index nOBS,
                                   const Index nCABS) {
        const Index N = nOBS + nCABS;
        assert(N == S.rows() && N == F_AO.rows());

        //Inverse of HF block of matrix
        const MatrixXd SInv = S.block(0, 0, nOBS, nOBS).inverse();

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

        //Canonical orthogonalisation with normalisation. The OBS directions have been projected out so the
        //projected overlap has (at least) nOBS null eigenvalues; only the remaining directions span the CABS
        const MatrixXd S_tmp = P.transpose() * S * P;
        SelfAdjointEigenSolver<MatrixXd> orthogonalEig(S_tmp);
        MatrixXd U = orthogonalEig.eigenvectors().real();
        VectorXd s = orthogonalEig.eigenvalues().real();

        std::vector<Index> cabsDirections;
        for (Index i = 0; i < s.rows(); ++i)
            if (abs(s(i)) > 1E-8)
                cabsDirections.push_back(i);

        if ((Index) cabsDirections.size() > nCABS)
            throw std::runtime_error("Failed to project the orbital basis set out of the CABS");

        MatrixXd X = MatrixXd::Zero(N, (Index) cabsDirections.size());
        for (Index k = 0; k < (Index) cabsDirections.size(); ++k)
            X.col(k) = U.col(cabsDirections[k]) / std::sqrt(s(cabsDirections[k]));

        const MatrixXd P_Prime = P * X;

        //Canonicalise the CABS with the Fock operator
        const MatrixXd F_Prime = P_Prime.transpose() * F_AO * P_Prime;
        SelfAdjointEigenSolver<MatrixXd> fockEig(F_Prime);

        //Matrix which transforms both the OBS and CABS into the projected canonicalised basis set
        //This matrix is different to MRCC's as the degenerate eigenvectors are different
        //but any observable property will be the same. Should the CABS be linearly dependent on the OBS, the
        //redundant CABS functions are left as zero columns
        MatrixXd Q = MatrixXd::Zero(N, N);
        Q.block(0, 0, nOBS, nOBS) = C_HF;
        Q.block(0, nOBS, N, (Index) cabsDirections.size()) = P_Prime * fockEig.eigenvectors().real();
        return Q;
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
