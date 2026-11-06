#include <Quantum/Q_PCH.h>
#include <libint2.hpp>
#include <atomic>
#include <exception>
//
#include "IntegralsF12.h"
#include "Quantum/Basis/BasisSet.h"
#include "Tensor/TensorEinsum.h"
#include "Tensor/TensorOperations.h"
#include "Dahlia/Tools/Threading.h"

//
//
namespace SCESS::Integrals::F12 {
    //
    //
    using namespace Eigen;
    using libint2::Shell;
    //
    //
    namespace {
        //Index of the first function of each shell, and the total number of functions, over the first nShells shells
        std::vector<Index> FirstFunctions(const std::vector<Shell> &shells, size_t nShells, Index *nFuncs = nullptr) {
            std::vector<Index> first(nShells);
            Index n = 0;
            for (size_t s = 0; s < nShells; s++) {
                first[s] = n;
                n += (Index) shells[s].size();
            }
            if (nFuncs)
                *nFuncs = n;
            return first;
        }
        //
        libint2::Engine MakeEngine(const std::vector<Shell> &shells, Operator2e op, double gamma) {
            using libint2::Operator;
            size_t maxNprim = 0;
            int maxL = 0;
            for (const Shell &shell: shells) {
                maxNprim = std::max(maxNprim, shell.nprim());
                for (const auto &c: shell.contr)
                    maxL = std::max(maxL, c.l);
            }
            constexpr double precision = std::numeric_limits<double>::epsilon();

            switch (op) {
                case Operator2e::Coulomb:
                    return {Operator::coulomb, maxNprim, maxL, 0, precision};
                case Operator2e::F12:
                    return {Operator::cgtg, maxNprim, maxL, 0, precision, CalculateSlaterContraction(gamma)};
                case Operator2e::F12Squared:
                    return {
                        Operator::cgtg, maxNprim, maxL, 0, precision,
                        CalculateSlaterSquaredContraction(CalculateSlaterContraction(gamma))
                    };
                case Operator2e::F12Coulomb:
                    return {Operator::cgtg_x_coulomb, maxNprim, maxL, 0, precision, CalculateSlaterContraction(gamma)};
                case Operator2e::DelF12Squared:
                    return {Operator::delcgtg2, maxNprim, maxL, 0, precision, CalculateSlaterContraction(gamma)};
            }
            throw std::invalid_argument("Unknown two-electron operator");
        }
        //
        //One engine per thread: libint2 engines are not thread-safe
        std::vector<libint2::Engine> MakeEngines(const std::vector<Shell> &shells, Operator2e op, double gamma) {
            std::vector<libint2::Engine> engines;
            engines.reserve(NumComputeThreads());
            for (int t = 0; t < NumComputeThreads(); t++)
                engines.push_back(MakeEngine(shells, op, gamma));
            return engines;
        }
    }
    //
    //
    TEITensor AOTwoElectronTensor(const std::vector<Shell> &shells, size_t nShells, Operator2e op, double gamma,
                                  const DH::Thread *thread) {
        Index nbf = 0;
        const std::vector<Index> first = FirstFunctions(shells, nShells, &nbf);

        TEITensor tei(nbf, nbf, nbf, nbf);
        // Screened shell quartets are skipped below so must already be zero
        tei.setZero();

        std::vector<libint2::Engine> engines = MakeEngines(shells, op, gamma);

        //Unique shell quartets s1 >= s2, s1 >= s3 >= s4 and (s1 s2) >= (s3 s4). Every quartet is written by the item
        //of its largest shell alone, so the threads never write to the same elements. The largest (most expensive)
        //items are handed out first.
        ParallelForDynamic(nShells, [&](size_t item, int threadIdx) {
            if (thread && thread->ShouldStop())
                return;
            const size_t s1 = nShells - 1 - item;
            libint2::Engine &engine = engines[threadIdx];
            const auto &results = engine.results();

            const Index n1 = (Index) shells[s1].size(), bf1 = first[s1];
            for (size_t s2 = 0; s2 <= s1; s2++) {
                const Index n2 = (Index) shells[s2].size(), bf2 = first[s2];
                for (size_t s3 = 0; s3 <= s1; s3++) {
                    const Index n3 = (Index) shells[s3].size(), bf3 = first[s3];
                    const size_t s4Max = s3 == s1 ? s2 : s3;
                    for (size_t s4 = 0; s4 <= s4Max; s4++) {
                        const Index n4 = (Index) shells[s4].size(), bf4 = first[s4];

                        engine.compute(shells[s1], shells[s2], shells[s3], shells[s4]);
                        const double *ints = results[0];
                        if (ints == nullptr)
                            continue;

                        for (Index f1 = 0, f1234 = 0; f1 < n1; f1++) {
                            const Index i = bf1 + f1;
                            for (Index f2 = 0; f2 < n2; f2++) {
                                const Index j = bf2 + f2;
                                for (Index f3 = 0; f3 < n3; f3++) {
                                    const Index k = bf3 + f3;
                                    for (Index f4 = 0; f4 < n4; f4++, f1234++) {
                                        const Index l = bf4 + f4;
                                        const double value = ints[f1234];
                                        tei(i, j, k, l) = tei(j, i, k, l) = tei(i, j, l, k) = tei(j, i, l, k) = value;
                                        tei(k, l, i, j) = tei(l, k, i, j) = tei(k, l, j, i) = tei(l, k, j, i) = value;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        });

        if (thread && thread->ShouldStop())
            return {};
        return tei;
    }
    //
    TEITensor OBSTwoElectronRepulsionTensor(const F12BasisSet &basisSet) {
        return AOTwoElectronTensor(basisSet.GetShells(), basisSet.GetNoHFShells(), Operator2e::Coulomb);
    }
    //
    TEITensor HalfTransformTwoElectron(const F12BasisSet &basisSet, Operator2e op,
                                       const MatrixXd &C1, AOSpace space2,
                                       const MatrixXd &C3, AOSpace space4,
                                       double gamma) {
        const std::vector<Shell> &shells = basisSet.GetShells();
        const size_t nShells = shells.size();
        Index N = 0;
        const std::vector<Index> first = FirstFunctions(shells, nShells, &N);

        if (C1.rows() != N || C3.rows() != N)
            throw std::invalid_argument("HalfTransformTwoElectron: the coefficients need one row per basis function");

        //Shells (and number of functions) spanned by an AO index
        auto Shells = [&](AOSpace space) { return space == AOSpace::OBS ? basisSet.GetNoHFShells() : nShells; };
        auto Funcs = [&](AOSpace space) { return space == AOSpace::OBS ? (Index) basisSet.GetNoHFFuncs() : N; };
        const size_t nShells2 = Shells(space2), nShells4 = Shells(space4);
        const Index N2 = Funcs(space2), N4 = Funcs(space4);

        //Shells with a non-zero coefficient
        auto ActiveShells = [&](const MatrixXd &C) {
            std::vector<size_t> active;
            for (size_t s = 0; s < nShells; s++)
                if (C.cols() > 0 && C.middleRows(first[s], (Index) shells[s].size()).cwiseAbs().maxCoeff() > 0.0)
                    active.push_back(s);
            return active;
        };
        const std::vector<size_t> active1 = ActiveShells(C1), active3 = ActiveShells(C3);

        //(mu nu|lam sig) = (lam sig|mu nu) so when the two electrons are transformed alike only the quartets with
        //s1 <= s3 are needed: H is then the half with s1 < s3, plus half the s1 = s3 part, plus its transpose
        const bool symmetric = space2 == space4 && C1.cols() == C3.cols() && C1 == C3;

        std::vector<std::pair<size_t, size_t> > pairs13;
        for (size_t a = 0; a < active1.size(); a++)
            for (size_t b = symmetric ? a : 0; b < active3.size(); b++)
                pairs13.emplace_back(active1[a], active3[b]);

        const Index P = C1.cols(), R = C3.cols();
        //Row-major copies so that the coefficients of one AO are contiguous
        const Matrix<double, Dynamic, Dynamic, RowMajor> C1r = C1, C3r = C3;

        TEITensor H(P, N2, R, N4);
        H.setZero();
        double *h = H.data();

        std::vector<libint2::Engine> engines = MakeEngines(shells, op, gamma);
        std::vector<std::vector<double> > buffers(NumComputeThreads());

        //Each item owns the block of H with nu in shell s2 and sig in shell s4, so the threads never write to the same
        //elements
        ParallelForDynamic(nShells2 * nShells4, [&](size_t item, int thread) {
            const size_t s2 = item / nShells4, s4 = item % nShells4;
            const Index n2 = (Index) shells[s2].size(), bf2 = first[s2];
            const Index n4 = (Index) shells[s4].size(), bf4 = first[s4];

            libint2::Engine &engine = engines[thread];
            const auto &results = engine.results();
            std::vector<double> &Y = buffers[thread];

            for (const auto &[s1, s3]: pairs13) {
                const Index n1 = (Index) shells[s1].size(), bf1 = first[s1];
                const Index n3 = (Index) shells[s3].size(), bf3 = first[s3];
                const double weight = symmetric && s1 == s3 ? 0.5 : 1.0;

                engine.compute(shells[s1], shells[s2], shells[s3], shells[s4]);
                const double *ints = results[0];
                if (ints == nullptr)
                    continue;

                //Y(p, f2, f3, f4) = sum_f1 C1(mu, p) (mu nu|lam sig)
                Y.assign(P * n2 * n3 * n4, 0.0);
                for (Index f1 = 0, f1234 = 0; f1 < n1; f1++) {
                    const double *c1 = C1r.data() + (bf1 + f1) * P;
                    for (Index f2 = 0; f2 < n2; f2++)
                        for (Index f3 = 0; f3 < n3; f3++)
                            for (Index f4 = 0; f4 < n4; f4++, f1234++) {
                                const double value = weight * ints[f1234];
                                if (value == 0.0)
                                    continue;
                                double *y = Y.data() + P * (f2 + n2 * (f3 + n3 * f4));
                                for (Index p = 0; p < P; p++)
                                    y[p] += c1[p] * value;
                            }
                }

                //H(p, nu, r, sig) += sum_f3 C3(lam, r) Y(p, f2, f3, f4)
                for (Index f4 = 0; f4 < n4; f4++)
                    for (Index f3 = 0; f3 < n3; f3++) {
                        const double *c3 = C3r.data() + (bf3 + f3) * R;
                        for (Index f2 = 0; f2 < n2; f2++) {
                            const double *y = Y.data() + P * (f2 + n2 * (f3 + n3 * f4));
                            for (Index r = 0; r < R; r++) {
                                const double c = c3[r];
                                if (c == 0.0)
                                    continue;
                                double *hOut = h + P * ((bf2 + f2) + N2 * (r + R * (bf4 + f4)));
                                for (Index p = 0; p < P; p++)
                                    hOut[p] += c * y[p];
                            }
                        }
                    }
            }
        });

        if (symmetric) {
            enum { p, nu, r, sig };
            //p nu r sig -> r sig p nu, evaluated into a new tensor as H would otherwise be read while overwritten
            return Evaluate(H + Einsum<Indices<p, nu, r, sig>, Indices<r, sig, p, nu> >(H));
        }

        return H;
    }
    //
    TEITensor FinishTransform(const TEITensor &halfTransformed, const MatrixXd &C2, const MatrixXd &C4) {
        enum { p, nu, q, r, sig, s };
        //p nu r sig, sig s -> p nu r s
        const TEITensor quarter = Evaluate(Einsum<Indices<p, nu, r, sig>, Indices<sig, s>, Indices<p, nu, r, s> >(
            halfTransformed, MatrixToTensor(C4)));
        //p nu r s, nu q -> p q r s
        return Evaluate(Einsum<Indices<p, nu, r, s>, Indices<nu, q>, Indices<p, q, r, s> >(
            quarter, MatrixToTensor(C2)));
    }
    //
    MatrixXd CoulombMatrix(const F12BasisSet &basisSet, const MatrixXd &D_OBS) {
        const std::vector<Shell> &shells = basisSet.GetShells();
        const size_t nShells = shells.size(), nOBSShells = basisSet.GetNoHFShells();
        Index N = 0;
        const std::vector<Index> first = FirstFunctions(shells, nShells, &N);

        MatrixXd J = MatrixXd::Zero(N, N);
        std::vector<libint2::Engine> engines = MakeEngines(shells, Operator2e::Coulomb, 1.0);

        //Each item (s1 >= s2) owns the J blocks (s1, s2) and (s2, s1)
        std::vector<std::pair<size_t, size_t> > pairs12;
        for (size_t s1 = 0; s1 < nShells; s1++)
            for (size_t s2 = 0; s2 <= s1; s2++)
                pairs12.emplace_back(s1, s2);

        ParallelForDynamic(pairs12.size(), [&](size_t item, int thread) {
            const auto [s1, s2] = pairs12[item];
            const Index n1 = (Index) shells[s1].size(), bf1 = first[s1];
            const Index n2 = (Index) shells[s2].size(), bf2 = first[s2];
            libint2::Engine &engine = engines[thread];
            const auto &results = engine.results();

            MatrixXd block = MatrixXd::Zero(n1, n2);
            //(mu nu|lam sig) = (mu nu|sig lam) and D is symmetric so s3 > s4 counts for both orders
            for (size_t s3 = 0; s3 < nOBSShells; s3++)
                for (size_t s4 = 0; s4 <= s3; s4++) {
                    const Index n3 = (Index) shells[s3].size(), bf3 = first[s3];
                    const Index n4 = (Index) shells[s4].size(), bf4 = first[s4];
                    const double factor = s3 == s4 ? 1.0 : 2.0;

                    engine.compute(shells[s1], shells[s2], shells[s3], shells[s4]);
                    const double *ints = results[0];
                    if (ints == nullptr)
                        continue;

                    for (Index f1 = 0, f1234 = 0; f1 < n1; f1++)
                        for (Index f2 = 0; f2 < n2; f2++) {
                            double sum = 0.0;
                            for (Index f3 = 0; f3 < n3; f3++)
                                for (Index f4 = 0; f4 < n4; f4++, f1234++)
                                    sum += ints[f1234] * D_OBS(bf3 + f3, bf4 + f4);
                            block(f1, f2) += factor * sum;
                        }
                }

            J.block(bf1, bf2, n1, n2) = block;
            J.block(bf2, bf1, n2, n1) = block.transpose();
        });

        return J;
    }
    //
    //
}
