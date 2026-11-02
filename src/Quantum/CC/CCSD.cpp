#include <Quantum/Q_PCH.h>
#include <iostream>
//
#include "Dahlia/Tools/Threading.h"
//
#include <cfloat>

#include "Quantum/CC/CC.h"
#include "Quantum/HF/HF.h"
#include "Quantum/Integrals/IntegralTransforms.h"

//
//
namespace SCESS::CC {
using namespace Eigen;
//
//
template <typename T>
using MatrixType = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
//
CCSD_Wfn CCSD(const HFWfn &rhfWfn, const Integrals::MolInts &mints,
              const DH::Thread *thread, std::stringstream *console) {
  // The delta change in energy which must be achieved for convergence
  constexpr double CC_CONVERGENCE = 1.0E-6;
  // Maximum coupled cluster iterations
  constexpr size_t MAX_ITERATIONS = 50;

  // Transform the two electron repulsion integrals into the MO basis set
  // for MP2 calculation. Integrals are in chemist's notation
  TEITensor eri =
      Integrals::Transforms::TransformTEI(mints.eri, rhfWfn.C_Alpha);
  // Evaluate MP2 energy
  double mp2Energy = HF::MP2_RHF_ClosedShell(rhfWfn, eri, thread, console);

  if (console) {
    *console << std::setprecision(12);
    *console << "--- Coupled Cluster Singles & Doubles (CCSD) ---\n";
  }

  // Antisymmetrise the two-electron integrals and integrate out the alpha and
  // beta spin functions for a closed-shell wavefunction. Also converts from
  // chemist's to physicits' notation
  eri = Integrals::Transforms::AntisymmetriseIntegrateOutSpinClosedShell(eri);

  // Number of occupied spin orbitals
  Index nOcc = std::distance(rhfWfn.Occ_Alpha.begin(),
                             std::ranges::find(rhfWfn.Occ_Alpha, 0));
  // For two electrons in each orbital
  nOcc *= 2;
  // Number of virtual spin orbitals
  Index nVirt = rhfWfn.Occ_Alpha.size() * 2 - nOcc;

  // Fock matrix in canonical MO basis set
  FockTensor F = TransformFockToMO(rhfWfn);

  // Denominator arrays
  Tensor<double, 2> Dia(nOcc, nVirt);
  Tensor<double, 4> Dijab(nOcc, nOcc, nVirt, nVirt);
  CreateDenominatorTensors(F, nOcc, nVirt, Dia, Dijab);

  // The first excitation amplitudes tensor
  // Dimensions of o,v
  T1Tensor t1Amps(nOcc, nVirt);
  t1Amps.setZero();

  // The second excitation amplitudes tensor
  // Dimensions of o,o,v,v
  T2Tensor t2Amps = CreateInitialT2Tensor(nOcc, nVirt, eri, Dijab);

  double E_CC_last = DBL_MAX, E_CC = 0.0;
  size_t ccIteration = 0;
  bool converged = false;

  for (; ccIteration < MAX_ITERATIONS; ccIteration++) {
    if (thread)
      if (thread->ShouldStop())
        return {};

    // Solve the CC equations
    IterateCCSD(eri, F, nOcc, nVirt, Dia, Dijab, t1Amps, t2Amps, thread);

    // Calculate CC correlation energy
    E_CC = CalculateEnergyCCSD(F, nOcc, nVirt, eri, t1Amps, t2Amps);

    double dE = abs(E_CC - E_CC_last);
    if (console)
      *console << "   CCSD(" << ccIteration << "): E(corr) = " << E_CC
               << "      dE = " << dE << '\n';

    if (dE < CC_CONVERGENCE) {
      converged = true;
      if (console)
        *console << "---- CCSD procedure converged! ----\n";
      break;
    }

    E_CC_last = E_CC;
  }

  if (!converged) {
    if (console)
      *console << "!!! CCSD procedure failed to converge !!!\n\n\n";
    return {};
  }

  if (console) {
    *console << "E(CCSD correlation) = " << E_CC << " hartrees\n";
    *console << "E(Total) = E(Hartree Fock) + E(CCSD) = " << rhfWfn.E_HF + E_CC
             << " hartrees\n\n\n";
  }

  return {E_CC, 0.0, t1Amps, t2Amps};
}
//
FockTensor TransformFockToMO(const HFWfn &rhfWfn) {
  // Transform the Fock matrix into the MO basis set
  MatrixXd FockMat =
      rhfWfn.C_Alpha.transpose() * rhfWfn.F_Alpha * rhfWfn.C_Alpha;

  // Assuming canonical HF orbitals only the diagonals are needed
  // Hence, we copy over only the diagonal elements.
  // Also, even indices are alpha spin and odd indices are beta spin
  FockTensor FockTensor(FockMat.rows() * 2, FockMat.rows() * 2);
  FockTensor.setZero();
  for (Index i = 0; i < FockTensor.dimension(0); i++)
    // Assume canonical HF orbitals
    FockTensor(i, i) = FockMat(i / 2, i / 2);

  return FockTensor;
}
//
void CreateDenominatorTensors(const FockTensor &F, const Index nOcc,
                              const Index nVirt, Eigen::Tensor<double, 2> &Dia,
                              Eigen::Tensor<double, 4> &Dijab) {
  // Slices of the Fock matricies for occupied and virtual orbitals
  Tensor<double, 2> Focc =
      F.slice(array<Index, 2>{0, 0}, array<Index, 2>{nOcc, nOcc});
  Tensor<double, 2> Fvirt =
      F.slice(array<Index, 2>{nOcc, nOcc}, array<Index, 2>{nVirt, nVirt});

  for (Index i = 0; i < nOcc; i++)
    for (Index a = 0; a < nVirt; a++)
      Dia(i, a) = Focc(i, i) - Fvirt(a, a);

  for (Index i = 0; i < nOcc; i++)
    for (Index j = 0; j < nOcc; j++)
      for (Index a = 0; a < nVirt; a++)
        for (Index b = 0; b < nVirt; b++)
          Dijab(i, j, a, b) =
              Focc(i, i) + Focc(j, j) - Fvirt(a, a) - Fvirt(b, b);
}
//
T2Tensor CreateInitialT2Tensor(const Index nOcc, const Index nVirt,
                               const TEITensor &moTEI,
                               const Eigen::Tensor<double, 4> &Dijab) {
  T2Tensor MOijab = moTEI.slice(array<Index, 4>{0, 0, nOcc, nOcc},
                                array<Index, 4>{nOcc, nOcc, nVirt, nVirt});
  return MOijab / Dijab;
}
//
void IterateCCSD(const TEITensor &moTEI, const FockTensor &F, const Index nOcc,
                 const Index nVirt, const Eigen::Tensor<double, 2> &Dia,
                 const Eigen::Tensor<double, 4> &Dijab, T1Tensor &t1,
                 T2Tensor &t2, const DH::Thread *thread) {
  // Intermediates
  auto MakeTildeTau = [&t1, &t2]() -> TauTensor {
    Tensor<double, 4> tildeTau = t2;

    // ia,jb-> iajb -> ijab
    Tensor<double, 4> tmp = 0.5 * t1.contract(t1, array<IndexPair<int>, 0>{})
                                      .shuffle(array<int, 4>{0, 2, 1, 3});
    tildeTau += tmp;
    tildeTau -= tmp.shuffle(array<int, 4>{0, 1, 3, 2});

    return tildeTau;
  };
  //
  auto MakeTau = [&t1, &t2]() -> TauTensor {
    TauTensor tau = t2;

    // ia,jb-> iajb -> ijab
    TauTensor tmp = t1.contract(t1, array<IndexPair<int>, 0>{})
                        .shuffle(array<int, 4>{0, 2, 1, 3});
    tau += tmp;
    tau -= tmp.shuffle(array<int, 4>{0, 1, 3, 2});

    return tau;
  };
  //
  auto MakeFae = [&moTEI, &F, &nOcc, &nVirt,
                  &t1](const TauTensor &tildeTau) -> FTensor {
    FTensor Fae(nVirt, nVirt);
    Fae.setZero();

    // me,ma-> ea -> ae
    Fae -=
        0.5 * F.slice(array<Index, 2>{0, nOcc}, array<Index, 2>{nOcc, nVirt})
                  .contract(t1, array<IndexPair<int>, 1>{IndexPair<int>{0, 0}})
                  .shuffle(array<int, 2>{1, 0});

    // MO[o,v,v,v]
    // mf,mafe -> ae
    Fae += t1.contract(
        moTEI.slice(array<Index, 4>{0, nOcc, nOcc, nOcc},
                    array<Index, 4>{nOcc, nVirt, nVirt, nVirt}),
        array<IndexPair<int>, 2>{IndexPair<int>{0, 0}, IndexPair<int>{1, 2}});

    // MO[o, o, v, v]
    // mnaf,mnef->ae
    Fae -= 0.5 * tildeTau.contract(
                     moTEI.slice(array<Index, 4>{0, 0, nOcc, nOcc},
                                 array<Index, 4>{nOcc, nOcc, nVirt, nVirt}),
                     array<IndexPair<int>, 3>{IndexPair<int>{0, 0},
                                              IndexPair<int>{1, 1},
                                              IndexPair<int>{3, 3}});

    return Fae;
  };
  //
  auto MakeFmi = [&moTEI, &F, &nOcc, &nVirt,
                  &t1](const TauTensor &tildeTau) -> FTensor {
    FTensor Fmi(nOcc, nOcc);
    Fmi.setZero();

    // F[o, v]
    // ie,me-> im -> mi
    Fmi +=
        0.5 * F.slice(array<Index, 2>{0, nOcc}, array<Index, 2>{nOcc, nVirt})
                  .contract(t1, array<IndexPair<int>, 1>{IndexPair<int>{1, 1}})
                  .shuffle(array<int, 2>{1, 0});

    // MO[o, o, o, v]
    // ne,mnie->mi
    Fmi += t1.contract(
        moTEI.slice(array<Index, 4>{0, 0, 0, nOcc},
                    array<Index, 4>{nOcc, nOcc, nOcc, nVirt}),
        array<IndexPair<int>, 2>{IndexPair<int>{0, 1}, IndexPair<int>{1, 3}});

    // MO[o, o, v, v]
    // inef,mnef-> im -> mi
    Fmi += 0.5 *
           tildeTau
               .contract(moTEI.slice(array<Index, 4>{0, 0, nOcc, nOcc},
                                     array<Index, 4>{nOcc, nOcc, nVirt, nVirt}),
                         array<IndexPair<int>, 3>{IndexPair<int>{1, 1},
                                                  IndexPair<int>{2, 2},
                                                  IndexPair<int>{3, 3}})
               .shuffle(array<int, 2>{1, 0});

    return Fmi;
  };
  //
  auto MakeFme = [&moTEI, &F, &nOcc, &nVirt, &t1]() -> FTensor {
    FTensor Fme =
        F.slice(array<Index, 2>{0, nOcc}, array<Index, 2>{nOcc, nVirt});
    // nf,mnef->me
    Fme += t1.contract(
        moTEI.slice(array<Index, 4>{0, 0, nOcc, nOcc},
                    array<Index, 4>{nOcc, nOcc, nVirt, nVirt}),
        array<IndexPair<int>, 2>{IndexPair<int>{0, 1}, IndexPair<int>{1, 3}});

    return Fme;
  };
  //
  auto MakeWmnij = [&moTEI, &nOcc, &nVirt,
                    &t1](const TauTensor &tau) -> WTensor {
    WTensor Wmnij = moTEI.slice(array<Index, 4>{0, 0, 0, 0},
                                array<Index, 4>{nOcc, nOcc, nOcc, nOcc});

    // je,mnie-> jmni -> mnij
    WTensor slicedMO = moTEI.slice(array<Index, 4>{0, 0, 0, nOcc},
                                   array<Index, 4>{nOcc, nOcc, nOcc, nVirt});
    WTensor Pij =
        t1.contract(slicedMO, array<IndexPair<int>, 1>{IndexPair<int>{1, 3}})
            .shuffle(array<int, 4>{1, 2, 3, 0});
    Wmnij += Pij;
    Wmnij -= Pij.shuffle(array<int, 4>{0, 1, 3, 2});

    // ijef,mnef-> ijmn -> mnij
    Wmnij +=
        0.25 *
        tau.contract(moTEI.slice(array<Index, 4>{0, 0, nOcc, nOcc},
                                 array<Index, 4>{nOcc, nOcc, nVirt, nVirt}),
                     array<IndexPair<int>, 2>{IndexPair<int>{2, 2},
                                              IndexPair<int>{3, 3}})
            .shuffle(array<int, 4>{2, 3, 0, 1});

    return Wmnij;
  };
  //
  auto MakeWabef = [&moTEI, &nOcc, &nVirt,
                    &t1](const TauTensor &tau) -> WTensor {
    WTensor Wabef = moTEI.slice(array<Index, 4>{nOcc, nOcc, nOcc, nOcc},
                                array<Index, 4>{nVirt, nVirt, nVirt, nVirt});

    // mb,amef-> baef -> abef
    WTensor Pab =
        t1.contract(moTEI.slice(array<Index, 4>{nOcc, 0, nOcc, nOcc},
                                array<Index, 4>{nVirt, nOcc, nVirt, nVirt}),
                    array<IndexPair<int>, 1>{IndexPair<int>{0, 1}})
            .shuffle(array<int, 4>{1, 0, 2, 3});

    Wabef -= Pab;
    Wabef += Pab.shuffle(array<int, 4>{1, 0, 2, 3});

    // MO[o, o, v, v]
    // mnab,mnef->abef
    Wabef += 0.25 * tau.contract(
                        moTEI.slice(array<Index, 4>{0, 0, nOcc, nOcc},
                                    array<Index, 4>{nOcc, nOcc, nVirt, nVirt}),
                        array<IndexPair<int>, 2>{IndexPair<int>{0, 0},
                                                 IndexPair<int>{1, 1}});

    return Wabef;
  };
  //
  auto MakeWmbej = [&moTEI, &nOcc, &nVirt, &t1, &t2]() -> WTensor {
    WTensor Wmbej = moTEI.slice(array<Index, 4>{0, nOcc, nOcc, 0},
                                array<Index, 4>{nOcc, nVirt, nVirt, nOcc});

    // jf,mbef-> jmbe -> mbej
    Wmbej +=
        t1.contract(moTEI.slice(array<Index, 4>{0, nOcc, nOcc, nOcc},
                                array<Index, 4>{nOcc, nVirt, nVirt, nVirt}),
                    array<IndexPair<int>, 1>{IndexPair<int>{1, 3}})
            .shuffle(array<int, 4>{1, 2, 3, 0});

    // nb,mnej-> bmej -> mbej
    Wmbej -= t1.contract(moTEI.slice(array<Index, 4>{0, 0, nOcc, 0},
                                     array<Index, 4>{nOcc, nOcc, nVirt, nOcc}),
                         array<IndexPair<int>, 1>{IndexPair<int>{0, 1}})
                 .shuffle(array<int, 4>{1, 0, 2, 3});

    // jf,nb-> jfnb -> jnfb
    WTensor tmp = (0.5 * t2) + t1.contract(t1, array<IndexPair<int>, 0>{})
                                   .shuffle(array<int, 4>{0, 2, 1, 3});

    // MO[o, o, v, v]
    // jnfb,mnef-> jbme -> mbej
    Wmbej -=
        tmp.contract(moTEI.slice(array<Index, 4>{0, 0, nOcc, nOcc},
                                 array<Index, 4>{nOcc, nOcc, nVirt, nVirt}),
                     array<IndexPair<int>, 2>{IndexPair<int>{1, 1},
                                              IndexPair<int>{2, 3}})
            .shuffle(array<int, 4>{2, 1, 3, 0});

    return Wmbej;
  };

  if (thread)
    if (thread->ShouldStop())
      return;

  TauTensor tildeTau = MakeTildeTau(), tau = MakeTau();

  FTensor Fae = MakeFae(tildeTau), Fmi = MakeFmi(tildeTau), Fme = MakeFme();

  WTensor Wmnij = MakeWmnij(tau), Wabef = MakeWabef(tau), Wmbej = MakeWmbej();

  if (thread)
    if (thread->ShouldStop())
      return;

  // Make the t1 RHS
  // Dimensions of nOcc and nVirt
  T1Tensor t1RHS =
      F.slice(array<Index, 2>{0, nOcc}, array<Index, 2>{nOcc, nVirt});

  // ie,ae -> ia
  t1RHS += t1.contract(Fae, array<IndexPair<int>, 1>{IndexPair<int>{1, 1}});
  // ma,mi -> ai -> ia
  t1RHS -= t1.contract(Fmi, array<IndexPair<int>, 1>{IndexPair<int>{0, 0}})
               .shuffle(array<int, 2>{1, 0});

  // imae,me -> ia
  t1RHS += t2.contract(Fme, array<IndexPair<int>, 2>{IndexPair<int>{1, 0},
                                                     IndexPair<int>{3, 1}});

  // nf,naif -> ai -> ia
  t1RHS -= t1.contract(moTEI.slice(array<Index, 4>{0, nOcc, 0, nOcc},
                                   array<Index, 4>{nOcc, nVirt, nOcc, nVirt}),
                       array<IndexPair<int>, 2>{IndexPair<int>{0, 0},
                                                IndexPair<int>{1, 3}})
               .shuffle(array<int, 2>{1, 0});

  // imef, maef -> ia
  t1RHS -=
      0.5 * t2.contract(moTEI.slice(array<Index, 4>{0, nOcc, nOcc, nOcc},
                                    array<Index, 4>{nOcc, nVirt, nVirt, nVirt}),
                        array<IndexPair<int>, 3>{IndexPair<int>{1, 0},
                                                 IndexPair<int>{2, 2},
                                                 IndexPair<int>{3, 3}});

  // mnae, nmei -> ai -> ia
  t1RHS -=
      0.5 * t2.contract(moTEI.slice(array<Index, 4>{0, 0, nOcc, 0},
                                    array<Index, 4>{nOcc, nOcc, nVirt, nOcc}),
                        array<IndexPair<int>, 3>{IndexPair<int>{0, 1},
                                                 IndexPair<int>{1, 0},
                                                 IndexPair<int>{3, 2}})
                .shuffle(array<int, 2>{1, 0});

  if (thread)
    if (thread->ShouldStop())
      return;

  // Make the t2 RHS
  T2Tensor t2RHS = moTEI.slice(array<Index, 4>{0, 0, nOcc, nOcc},
                               array<Index, 4>{nOcc, nOcc, nVirt, nVirt});
  // Fae - 1/2 Sum_m(tmb * Fme)
  // mb,me->be
  T1Tensor tmp =
      Fae - 0.5 * Tensor<double, 2>(t1.contract(
                      Fme, array<IndexPair<int>, 1>{IndexPair<int>{0, 0}}));

  // Do ijae,be->ijab = 0123,01->0123
  T2Tensor Pab =
      t2.contract(tmp, array<IndexPair<int>, 1>{IndexPair<int>{3, 1}});
  t2RHS += Pab;
  t2RHS -= Pab.shuffle(array<int, 4>{0, 1, 3, 2});

  // me,je->mj
  tmp = Fmi +
        0.5 * Fme.contract(t1, array<IndexPair<int>, 1>{IndexPair<int>{1, 1}});
  // Do imab,mj->iabj->ijab
  T2Tensor Pij =
      t2.contract(tmp, array<IndexPair<int>, 1>{IndexPair<int>{1, 0}})
          .shuffle(array<int, 4>{0, 3, 1, 2});
  t2RHS -= Pij;
  t2RHS += Pij.shuffle(array<int, 4>{1, 0, 2, 3});

  // mnij,mnab->ijab
  t2RHS +=
      0.5 * Wmnij.contract(tau, array<IndexPair<int>, 2>{IndexPair<int>{0, 0},
                                                         IndexPair<int>{1, 1}});
  // ijef,abef->ijab
  t2RHS +=
      0.5 * tau.contract(Wabef, array<IndexPair<int>, 2>{IndexPair<int>{2, 2},
                                                         IndexPair<int>{3, 3}});

  // mbej,ie -> mbji,ma->bjia -> ijab
  T2Tensor tmp2 =
      Tensor<double, 4>(moTEI.slice(array<Index, 4>{0, nOcc, nOcc, 0},
                                    array<Index, 4>{nOcc, nVirt, nVirt, nOcc}))
          .contract(t1, array<IndexPair<int>, 1>{IndexPair<int>{2, 1}})
          .contract(t1, array<IndexPair<int>, 1>{IndexPair<int>{0, 0}})
          .shuffle(array<int, 4>{2, 1, 3, 0});

  // imae,mbej->iabj->ijab
  T2Tensor Pijab =
      t2.contract(Wmbej, array<IndexPair<int>, 2>{IndexPair<int>{1, 0},
                                                  IndexPair<int>{3, 2}})
          .shuffle(array<int, 4>{0, 3, 1, 2});
  Pijab -= tmp2;
  //
  t2RHS += Pijab;
  t2RHS -= Pijab.shuffle(array<int, 4>{0, 1, 3, 2});
  t2RHS -= Pijab.shuffle(array<int, 4>{1, 0, 2, 3});
  t2RHS += Pijab.shuffle(array<int, 4>{1, 0, 3, 2});

  // abej,ie -> abji -> ijab
  Pij =
      Tensor<double, 4>(moTEI.slice(array<Index, 4>{nOcc, nOcc, nOcc, 0},
                                    array<Index, 4>{nVirt, nVirt, nVirt, nOcc}))
          .contract(t1, array<IndexPair<int>, 1>{IndexPair<int>{2, 1}})
          .shuffle(array<int, 4>{3, 2, 0, 1});

  t2RHS += Pij;
  t2RHS -= Pij.shuffle(array<int, 4>{1, 0, 2, 3});

  // mbij,ma -> bija -> ijab
  Pab = Tensor<double, 4>(moTEI.slice(array<Index, 4>{0, nOcc, 0, 0},
                                      array<Index, 4>{nOcc, nVirt, nOcc, nOcc}))
            .contract(t1, array<IndexPair<int>, 1>{IndexPair<int>{0, 0}})
            .shuffle(array<int, 4>{1, 2, 3, 0});

  t2RHS -= Pab;
  t2RHS += Pab.shuffle(array<int, 4>{0, 1, 3, 2});

  t1 = t1RHS / Dia;
  t2 = t2RHS / Dijab;
}
//
double CalculateEnergyCCSD(const FockTensor &F, const Index nOcc,
                           const Index nVirt, const TEITensor &moTEI,
                           const T1Tensor &t1, const T2Tensor &t2) {
  Tensor<double, 0> correlationEnergyTensor =
      F.slice(array<Index, 2>{0, nOcc}, array<Index, 2>{nOcc, nVirt})
          .contract(t1, array<IndexPair<int>, 2>{IndexPair<int>{0, 0},
                                                 IndexPair<int>{1, 1}});
  Tensor<double, 4> moSlice =
      moTEI.slice(array<Index, 4>{0, 0, nOcc, nOcc},
                  array<Index, 4>{nOcc, nOcc, nVirt, nVirt});
  correlationEnergyTensor +=
      0.25 *
      moSlice.contract(t2, array<IndexPair<int>, 4>{
                               IndexPair<int>{0, 0}, IndexPair<int>{1, 1},
                               IndexPair<int>{2, 2}, IndexPair<int>{3, 3}});
  correlationEnergyTensor += 0.5 * moSlice
                                       .contract(t1,
                                                 array<IndexPair<int>, 2>{
                                                     IndexPair<int>{0, 0},
                                                     IndexPair<int>{2, 1},
                                                 })
                                       .contract(t1, array<IndexPair<int>, 2>{
                                                         IndexPair<int>{0, 0},
                                                         IndexPair<int>{1, 1},
                                                     });

  return correlationEnergyTensor(0);
}

//
} // namespace SCESS::CC
