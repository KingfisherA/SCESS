#include <Quantum/Q_PCH.h>
//
#include "Dahlia/Tools/Threading.h"
//
#include "Quantum/CC/CC.h"
#include "Quantum/HF/HF.h"
#include "Quantum/Integrals/IntegralTransforms.h"
#include "Tensor/TensorParallel.h"

//Closed-shell perturbative triples correction, as in PySCF's ccsd_t_slow for RHF references. For each unique triple
//of virtual orbitals a >= b >= c the connected (W) and disconnected (V) triples are built for every i, j, k in o^3
//work arrays, so the memory is that of the ovvv integrals rather than of the o^3 v^3 triples amplitudes.
namespace SCESS::CC {
    using namespace Eigen;
    //
    //
    double RCCSD_TriplesCorrection(const TEITensor &eri,
                                   const VectorXd &orbitalEnergies,
                                   const Index nOcc,
                                   const T1Tensor &t1,
                                   const T2Tensor &t2,
                                   const DH::Thread *thread) {
        const Index o = nOcc;
        const Index v = eri.dimension(0) - nOcc;
        const Index o2 = o * o, o3 = o * o * o;

        //Integrals and amplitudes laid out so that the contractions of each triple are matrix products over
        //contiguous blocks. eri is <pq|rs> = (pr|qs).
        //P(i, f, a, b) = (ia|bf), an o x v matrix for each a, b
        Tensor<double, 4> P(o, v, v, v);
        //R(i, j, m, a) = (ia|mj), an o^2 x o matrix for each a
        Tensor<double, 4> R(o, o, o, v);
        //S(i, j, a, b) = (ia|jb), an o x o matrix for each a, b
        Tensor<double, 4> S(o, o, v, v);
        //Q(f, j, k, c) = t_kj^cf, a v x o^2 matrix for each c
        Tensor<double, 4> Q(v, o, o, v);
        for (Index b = 0; b < v; b++)
            for (Index a = 0; a < v; a++)
                for (Index f = 0; f < v; f++)
                    for (Index i = 0; i < o; i++)
                        P(i, f, a, b) = eri(i, o + b, o + a, o + f);
        for (Index a = 0; a < v; a++)
            for (Index m = 0; m < o; m++)
                for (Index j = 0; j < o; j++)
                    for (Index i = 0; i < o; i++)
                        R(i, j, m, a) = eri(i, m, o + a, j);
        for (Index b = 0; b < v; b++)
            for (Index a = 0; a < v; a++)
                for (Index j = 0; j < o; j++)
                    for (Index i = 0; i < o; i++)
                        S(i, j, a, b) = eri(i, j, o + a, o + b);
        for (Index c = 0; c < v; c++)
            for (Index k = 0; k < o; k++)
                for (Index j = 0; j < o; j++)
                    for (Index f = 0; f < v; f++)
                        Q(f, j, k, c) = t2(k, j, c, f);

        using MatMap = Map<MatrixXd>;
        using ConstMatMap = Map<const MatrixXd>;

        //w_ijk^abc = sum_f (ia|bf) t_kj^cf - sum_m (ia|mj) t_mk^bc into out(i, j, k)
        auto GetW = [&](Index a, Index b, Index c, double *out) {
            MatMap w(out, o, o2);
            w.noalias() = ConstMatMap(P.data() + o * v * (a + v * b), o, v) *
                          ConstMatMap(Q.data() + v * o2 * c, v, o2);
            MatMap(out, o2, o).noalias() -= ConstMatMap(R.data() + o3 * a, o2, o) *
                                             ConstMatMap(t2.data() + o2 * (b + v * c), o, o);
        };
        //v_ijk^abc = (ia|jb) t_k^c / 2 into out(i, j, k). The f_kc term vanishes for canonical orbitals.
        auto GetV = [&](Index a, Index b, Index c, double *out) {
            const ConstMatMap Sab(S.data() + o2 * (a + v * b), o2, 1);
            for (Index k = 0; k < o; k++)
                MatMap(out + o2 * k, o2, 1) = (0.5 * t1(k, c)) * Sab;
        };

        //Sums the six simultaneous permutations of the pairs (ia), (jb), (kc) of get(...) into out:
        //X_ijk^abc + X_ikj^acb + X_jik^bac + X_jki^bca + X_kij^cab + X_kji^cba
        auto Symmetrise = [&](Index a, Index b, Index c, auto &&get, double *out, double *scratch) {
            auto idx = [o, o2](Index i, Index j, Index k) { return i + o * j + o2 * k; };
            get(a, b, c, out);
            get(a, c, b, scratch);
            for (Index k = 0; k < o; k++)
                for (Index j = 0; j < o; j++)
                    for (Index i = 0; i < o; i++)
                        out[idx(i, j, k)] += scratch[idx(i, k, j)];
            get(b, a, c, scratch);
            for (Index k = 0; k < o; k++)
                for (Index j = 0; j < o; j++)
                    for (Index i = 0; i < o; i++)
                        out[idx(i, j, k)] += scratch[idx(j, i, k)];
            get(b, c, a, scratch);
            for (Index k = 0; k < o; k++)
                for (Index j = 0; j < o; j++)
                    for (Index i = 0; i < o; i++)
                        out[idx(i, j, k)] += scratch[idx(j, k, i)];
            get(c, a, b, scratch);
            for (Index k = 0; k < o; k++)
                for (Index j = 0; j < o; j++)
                    for (Index i = 0; i < o; i++)
                        out[idx(i, j, k)] += scratch[idx(k, i, j)];
            get(c, b, a, scratch);
            for (Index k = 0; k < o; k++)
                for (Index j = 0; j < o; j++)
                    for (Index i = 0; i < o; i++)
                        out[idx(i, j, k)] += scratch[idx(k, j, i)];
        };

        //The unique virtual triples a >= b >= c (every triple costs the same)
        std::vector<std::array<Index, 3> > triples;
        triples.reserve(v * (v + 1) * (v + 2) / 6);
        for (Index a = 0; a < v; a++)
            for (Index b = 0; b <= a; b++)
                for (Index c = 0; c <= b; c++)
                    triples.push_back({a, b, c});

        const int nThreads = NumComputeThreads();
        //Per thread: W, V, Z and a scratch array, each o^3
        std::vector<VectorXd> buffers(nThreads, VectorXd(4 * o3));
        std::vector<double> energies(nThreads, 0.0);
        const VectorXd &e = orbitalEnergies;

        ParallelForDynamic(triples.size(), [&](size_t item, int threadIdx) {
            if (thread && thread->ShouldStop())
                return;
            const auto [a, b, c] = triples[item];
            double *W = buffers[threadIdx].data();
            double *V = W + o3;
            double *Z = V + o3;
            double *scratch = Z + o3;

            Symmetrise(a, b, c, GetW, W, scratch);
            Symmetrise(a, b, c, GetV, V, scratch);

            //Unique triples stand for every ordering of a, b, c
            double degeneracy = 1.0;
            if (a == c)
                degeneracy = 6.0;
            else if (a == b || b == c)
                degeneracy = 2.0;

            //Y = W + V (V already holds the half of (ia|jb) t_k^c), then Z = r3(Y) / D with
            //r3(Y) = 4 Y_ijk + Y_kij + Y_jki - 2 Y_kji - 2 Y_ikj - 2 Y_jik. As r3 commutes with the simultaneous
            //permutations this equals PySCF's sum over the 36 pairs of unsymmetrised w and z.
            auto idx = [o, o2](Index i, Index j, Index k) { return i + o * j + o2 * k; };
            for (Index n = 0; n < o3; n++)
                scratch[n] = W[n] + V[n];
            const double eVirt = e(o + a) + e(o + b) + e(o + c);
            double energy = 0.0;
            for (Index k = 0; k < o; k++)
                for (Index j = 0; j < o; j++)
                    for (Index i = 0; i < o; i++) {
                        const double r3 = 4.0 * scratch[idx(i, j, k)] + scratch[idx(k, i, j)] +
                                          scratch[idx(j, k, i)] - 2.0 * scratch[idx(k, j, i)] -
                                          2.0 * scratch[idx(i, k, j)] - 2.0 * scratch[idx(j, i, k)];
                        const double d3 = (e(i) + e(j) + e(k) - eVirt) * degeneracy;
                        energy += W[idx(i, j, k)] * r3 / d3;
                    }
            energies[threadIdx] += energy;
        });

        if (thread && thread->ShouldStop())
            return 0.0;

        double E_T = 0.0;
        for (double energy: energies)
            E_T += energy;
        return 2.0 * E_T;
    }
    //
    RCCSD_Wfn RCCSD_T(const HFWfn &rhfWfn, const Integrals::MolInts &mints,
                      const DH::Thread *thread, std::stringstream *console) {
        const Index nOcc = std::distance(rhfWfn.Occ_Alpha.begin(), std::ranges::find(rhfWfn.Occ_Alpha, 0));
        const TEITensor eri = Integrals::Transforms::TransformTEI(mints.eri, rhfWfn.C_Alpha, true);

        if (console)
            *console << "--- Closed-shell CCSD with perturbative triples (CCSD(T)) ---\n";

        RCCSD_Wfn wfn = SolveRCCSD_DIIS(eri, rhfWfn.Ep_Alpha, nOcc, nullptr, thread, console);
        if (!wfn.converged || (thread && thread->ShouldStop())) {
            if (console)
                *console << "Cannot calculate the triples correction as the CCSD amplitudes did not converge.\n";
            return wfn;
        }

        if (console)
            *console << "Calculating the triples correction to the energy...\n";
        wfn.E_T = RCCSD_TriplesCorrection(eri, rhfWfn.Ep_Alpha, nOcc, wfn.t1, wfn.t2, thread);

        if (console)
            *console << std::setprecision(12)
                    << "E(T) = " << wfn.E_T << " hartrees\n"
                    << "E(Total) = E(Hartree Fock) + E(CCSD) + E((T)) = "
                    << rhfWfn.E_HF + wfn.E_CCSD + wfn.E_T << " hartrees\n\n\n";
        return wfn;
    }
    //
    //
}
