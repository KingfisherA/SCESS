#include <Quantum/Q_PCH.h>
//
#include <Eigen/Core>
#include <unsupported/Eigen/CXX11/Tensor>
//
#include "Dahlia/Tools/Threading.h"
#include "Quantum/CC/CC.h"
#include "Tensor/TensorEinsum.h"
#include "Tensor/TensorOperations.h"
#include "Tensor/TensorSlice.h"

//
//
namespace SCESS::CC {
using namespace Eigen;
//
//
CCSD_Wfn CCSD_T(const HFWfn &rhfResult, const Integrals::MolInts &mints,
                const DH::Thread *thread, std::stringstream *console) {
  if (console)
    *console << "--- Coupled Cluster Singles & Doubles with perturbative "
                "Triples (CCSD(T)) ---\n";

  // Perform CCSD with DIIS acceleration and get t1 and t2 amplitudes
  CCSD_Wfn ccsdWfn = CCSD_DIIS(rhfResult, mints, thread, console);

  // if (ccsdWfn == nullptr) {
  //     if (console)
  //         *console << "Cannot calculate triples correction as the CCSD t
  //         amplitudes did not converge.\n";
  //
  //     return nullptr;
  // }

  if (console)
    *console << "Calculating triples correction to the energy ... \n";

  // Number of occupied spin orbitals
  Index nOcc = std::distance(rhfResult.Occ_Alpha.begin(),
                             std::ranges::find(rhfResult.Occ_Alpha, 0));
  // For two electrons in each orbital
  nOcc *= 2;
  // Number of virtual orbitals
  Index nVirt = rhfResult.Occ_Alpha.size() * 2 - nOcc;

  // Get the Fock matrix in canonical MO basis set
  FockTensor F = TransformFockToMO(rhfResult);

  if (thread)
    if (thread->ShouldStop())
      return {};

  // Calculate the triples perturbation correction energy
  double E_T = CalculateTriplesCorrection(ccsdWfn, F, ccsdWfn.mints, nOcc,
                                          nVirt, thread);

  if (console) {
    *console << "E(T) = " << E_T << "\n";
    *console << "E(Total) = E(Hartree Fock) + E(CCSD) + E((T)) = "
             << rhfResult.E_HF + ccsdWfn.E_CCSD + E_T << '\n';
  }

  ccsdWfn.E_T = E_T;

  return ccsdWfn;
}
//
// Indices used for calculating the perturbative triples correction
enum {
  // Occupied orbitals
  I,
  J,
  K,
  L,
  M,
  N,
  // Unoccupied orbitals
  A,
  B,
  C,
  D,
  E,
  F
};
//
// Computes the triple correction energy with the spin integrated out
// Implemented from PROF. T. DANIEL CRAWFORD, VIRGINIA TECH Lecture 3
// eri = Must have spin integrated out and be the MO transformed tensor
double CalculateTriplesCorrection(const CCSD_Wfn &ccsd, const FockTensor &Fock,
                                  const TEITensor &eri, Index nOcc, Index nVirt,
                                  const DH::Thread *thread) {
  using IdxP = Eigen::IndexPair<int>;

  TensorSliceIndices idxs = {// All occupied orbital indices
                             {I, {0, nOcc}},
                             {J, {0, nOcc}},
                             {K, {0, nOcc}},
                             {L, {0, nOcc}},
                             {M, {0, nOcc}},
                             {N, {0, nOcc}},
                             // All virtual orbital indices
                             {A, {nOcc, nVirt}},
                             {B, {nOcc, nVirt}},
                             {C, {nOcc, nVirt}},
                             {D, {nOcc, nVirt}},
                             {E, {nOcc, nVirt}},
                             {F, {nOcc, nVirt}}};

  // Construct D; The denominator tensor composed of the diagonal elements
  // of the Fock matrix
  auto MakeDijkabc = [&Fock, &idxs, nOcc, nVirt]() {
    Tensor<double, 6> Dijkabc(nOcc, nOcc, nOcc, nVirt, nVirt, nVirt);

    // Slices of the Fock matricies for occupied and virtual orbitals
    Tensor<double, 2> Focc = Slice<Indices<I, J>>(idxs, Fock);
    Tensor<double, 2> Fvirt = Slice<Indices<A, B>>(idxs, Fock);

    for (Eigen::Index i = 0; i < nOcc; i++)
      for (Eigen::Index j = 0; j < nOcc; j++)
        for (Eigen::Index k = 0; k < nOcc; k++)
          for (Eigen::Index a = 0; a < nVirt; a++)
            for (Eigen::Index b = 0; b < nVirt; b++)
              for (Eigen::Index c = 0; c < nVirt; c++) {
                Dijkabc(i, j, k, a, b, c) = Focc(i, i) + Focc(j, j) +
                                            Focc(k, k) - Fvirt(a, a) -
                                            Fvirt(b, b) - Fvirt(c, c);
              }

    return Dijkabc;
  };

  auto MakeDijab = [&Fock, &idxs, nOcc, nVirt]() {
    Tensor<double, 4> Dijab(nOcc, nOcc, nVirt, nVirt);

    // Slices of the Fock matricies for occupied and virtual orbitals
    Tensor<double, 2> Focc = Slice<Indices<I, J>>(idxs, Fock);
    Tensor<double, 2> Fvirt = Slice<Indices<A, B>>(idxs, Fock);

    for (Eigen::Index i = 0; i < nOcc; i++)
      for (Eigen::Index j = 0; j < nOcc; j++)
        for (Eigen::Index a = 0; a < nVirt; a++)
          for (Eigen::Index b = 0; b < nVirt; b++) {
            Dijab(i, j, a, b) =
                Focc(i, i) + Focc(j, j) - Fvirt(a, a) - Fvirt(b, b);
          }

    return Dijab;
  };

  // Compute P(i/jk)P(a/bc) =
  // P(i/jk) = ijk - jik - kji
  // P(a/bc) = abc - bac - cba
  // (ijk - jik - kji) * (abc - bac - cba)
  // ijkabc - ijkbac - ijkcba
  // - jikabc + jikbac + jikcba
  // - kjiabc + kjibac + kjicba
  auto PerformPermute_ijk_abc =
      [&idxs](Tensor<double, 6> &tensor_ijkabc) -> void {
    Tensor<double, 6> result =
        // ijkabc - ijkbac - ijkcba
        tensor_ijkabc -
        Einsum<Indices<I, J, K, A, B, C>, Indices<I, J, K, B, A, C>>(
            tensor_ijkabc) -
        Einsum<Indices<I, J, K, A, B, C>, Indices<I, J, K, C, B, A>>(
            tensor_ijkabc)
        // - jikabc + jikbac + jikcba
        - Einsum<Indices<I, J, K, A, B, C>, Indices<J, I, K, A, B, C>>(
              tensor_ijkabc) +
        Einsum<Indices<I, J, K, A, B, C>, Indices<J, I, K, B, A, C>>(
            tensor_ijkabc) +
        Einsum<Indices<I, J, K, A, B, C>, Indices<J, I, K, C, B, A>>(
            tensor_ijkabc)
        // - kjiabc + kjibac + kjicba
        - Einsum<Indices<I, J, K, A, B, C>, Indices<K, J, I, A, B, C>>(
              tensor_ijkabc) +
        Einsum<Indices<I, J, K, A, B, C>, Indices<K, J, I, B, A, C>>(
            tensor_ijkabc) +
        Einsum<Indices<I, J, K, A, B, C>, Indices<K, J, I, C, B, A>>(
            tensor_ijkabc);
    tensor_ijkabc = result;
  };

  // Third order corrected double amplitudes
  Tensor<double, 4> t_ijab_3;
  // Second order correted triple amplitudes
  Tensor<double, 6> t_ijkabc;
  {
    {
      // g_{bc}^{dk} t_{ad}^{ij}
      Tensor<double, 6> lhs =
          Einsum<Indices<D, K, B, C>, Indices<I, J, A, D>,
                 Indices<I, J, K, A, B, C>, double, 4, double, 4>(
              Slice<Indices<D, K, B, C>>(idxs, eri), ccsd.t2);
      PerformPermute_ijk_abc(lhs);

      t_ijkabc = lhs;
      auto &rhs = lhs;

      // g_{lc}^{jk} t_{ab}^{il}
      rhs = Einsum<Indices<J, K, L, C>, Indices<I, L, A, B>,
                   Indices<I, J, K, A, B, C>>(
          Slice<Indices<J, K, L, C>>(idxs, eri), ccsd.t2);
      PerformPermute_ijk_abc(rhs);
      t_ijkabc -= rhs;
    }

    // Divide through by the denominator tensor
    t_ijkabc /= MakeDijkabc();

    Tensor<double, 4> lhs =
        Einsum<Indices<E, F, B, M>, Indices<I, J, M, A, E, F>,
               Indices<I, J, A, B>, double, 4, double, 6>(
            Slice<Indices<E, F, B, M>>(idxs, eri), t_ijkabc);

    // P(ab)ijab = ijab - ijba
    lhs = (lhs - Einsum<Indices<I, J, A, B>, Indices<I, J, B, A>>(lhs)).eval();

    t_ijab_3 = 0.5 * lhs;

    auto &rhs = lhs;
    rhs = Einsum<Indices<J, E, M, N>, Indices<I, M, N, A, B, E>,
                 Indices<I, J, A, B>, double, 4, double, 6>(
        Slice<Indices<J, B, M, N>>(idxs, eri), t_ijkabc);
    // P(ij)ijab = ijab - jiab
    rhs = (rhs - Einsum<Indices<I, J, A, B>, Indices<J, I, A, B>>(rhs)).eval();

    t_ijab_3 -= 0.5 * rhs;
    t_ijab_3 /= MakeDijab();
  }

  Tensor<double, 0> E_T_4 =
      Einsum<Indices<I, J, A, B>, Indices<I, J, A, B>, Indices<>, double, 4,
             double, 4>(Slice<Indices<I, J, A, B>>(idxs, eri), t_ijab_3);

  // Fourth order energy correction using the third order \hat{T}_2 (using
  // second order wavefunction)
  E_T_4(0) *= 1.0 / 4.0;

  // Calculate the fifth order energy correction using the fourth order
  // \hat{T}_1
  // E^{(5)} = < HF | \hat{H}^{(5)} | HF > = < HF | (\hat{V}^{(1)}
  // \hat{T}_1^{(4)}  )_c

  auto g_jkbc_with_t_ijkabc =
      Einsum<Indices<B, C, J, K>, Indices<I, J, K, A, B, C>, Indices<I, A>>(
          Slice<Indices<B, C, J, K>>(idxs, eri), t_ijkabc);
  Tensor<double, 0> E_ST_5 = Einsum<Indices<I, A>, Indices<I, A>, Indices<>>(
      g_jkbc_with_t_ijkabc, ccsd.t1);

  E_ST_5(0) *= 1.0 / 4.0;

  return E_T_4(0) + E_ST_5(0);
}
//
//
} // namespace SCESS::CC
