// Adapted from https://github.com/mkodrycka/mp2f12.git

#include "F12.h"
#include <Quantum/Q_PCH.h>
//
#include "Quantum/HF/HF.h"
#include "Quantum/Integrals/Integrals.h"
#include "Quantum/Integrals/IntegralsF12.h"

//
//
namespace SCESS::F12 {
using namespace Eigen;
//
//
MP2F12Result MP2_F12_Optimise(const HF_Wfn &refWfn,
                              const Integrals::MolInts &mints,
                              const Integrals::F12::F12Integrals &f12Ints,
                              F12Idxs f12Idxs, Tensor2 *F_out, TEITensor *V_out,
                              TEITensor *X_out, TEITensor *B_out) {
  // * We assume that the reference wavefunction is actually a closed shell RHF
  // * The integrals have been canonicalised and transformed into the CABS/OBS
  //   MO basis set and are in physicists' notation
  // * Spin is integrated out: each spatial pair ij has an opposite-spin
  //   component, and for i < j also a same-spin component
  const TEITensor &eri = mints.eri;
  const MP2F12Intermediates ints =
      EvaluateMP2F12Intermediates(mints, f12Ints, f12Idxs);
  const Tensor2 &F_ = ints.F;
  const TEITensor &V_ = ints.V;
  const TEITensor &X_ = ints.X;
  const TEITensor &B_ = ints.B;
  const TEITensor &C_ = ints.C;

  // Calculate conventional MP2 energy
  const double E_MP2 = EvaluateMP2Energy(eri, refWfn.Ep_Alpha, f12Idxs);

  const Index nOcc = f12Idxs.nOcc;
  const Index nVirt = f12Idxs.nVirt;

  // Geminal pairs kl: every ordered pair for opposite spins, but only k < l
  // for same spins as the antisymmetrised geminal |kl> - |lk> vanishes for
  // k = l and is redundant for k > l
  std::vector<std::pair<Index, Index>> oppositeSpinPairs, sameSpinPairs;
  for (Index k = 0; k < nOcc; k++)
    for (Index l = 0; l < nOcc; l++) {
      oppositeSpinPairs.emplace_back(k, l);
      if (k < l)
        sameSpinPairs.emplace_back(k, l);
    }

  // Minimise the Hylleraas functional for the geminal amplitudes of the
  // electron pair ij, returning the pair energy -V~^T B~^-1 V~ where
  // V~_kl = V_ij^kl - sum_ab C_kl^ab <ij|ab> / D_ij^ab and
  // B~_kl,mn = B_kl^mn - (e_i + e_j) X_kl^mn - sum_ab C_kl^ab C_mn^ab / D_ij^ab
  auto PairEnergy = [&](Index i, Index j, bool sameSpin) -> double {
    const auto &pairs = sameSpin ? sameSpinPairs : oppositeSpinPairs;
    const Index nPairs = (Index)pairs.size();
    const double epsilon_ij = F_(i, i) + F_(j, j);

    // Same-spin quantities are antisymmetrised in the geminal and virtual
    // indices
    auto V_kl = [&](Index k, Index l) {
      return sameSpin ? V_(i, j, k, l) - V_(i, j, l, k) : V_(i, j, k, l);
    };
    auto B_klmn = [&](Index k, Index l, Index m, Index n) {
      return sameSpin ? B_(k, l, m, n) - B_(k, l, n, m) -
                            epsilon_ij * (X_(k, l, m, n) - X_(k, l, n, m))
                      : B_(k, l, m, n) - epsilon_ij * X_(k, l, m, n);
    };
    auto C_klab = [&](Index k, Index l, Index a, Index b) {
      return sameSpin ? C_(k, l, a, b) - C_(k, l, b, a) : C_(k, l, a, b);
    };
    auto eri_ijab = [&](Index a, Index b) {
      return sameSpin ? eri(i, j, nOcc + a, nOcc + b) -
                            eri(i, j, nOcc + b, nOcc + a)
                      : eri(i, j, nOcc + a, nOcc + b);
    };

    VectorXd V_Tilde(nPairs);
    MatrixXd B_Tilde(nPairs, nPairs);
    for (Index p = 0; p < nPairs; p++) {
      const auto [k, l] = pairs[p];
      V_Tilde(p) = V_kl(k, l);
      for (Index q = 0; q < nPairs; q++) {
        const auto [m, n] = pairs[q];
        B_Tilde(p, q) = B_klmn(k, l, m, n);
      }
    }

    VectorXd C_ab(nPairs);
    for (Index a = 0; a < nVirt; a++)
      for (Index b = sameSpin ? a + 1 : 0; b < nVirt; b++) {
        const double D_ijab = F_(nOcc + a, nOcc + a) +
                              F_(nOcc + b, nOcc + b) - epsilon_ij;
        const double g_ijab = eri_ijab(a, b);

        for (Index p = 0; p < nPairs; p++)
          C_ab(p) = C_klab(pairs[p].first, pairs[p].second, a, b);

        V_Tilde -= C_ab * (g_ijab / D_ijab);
        B_Tilde -= C_ab * C_ab.transpose() / D_ijab;
      }

    return -V_Tilde.dot(B_Tilde.fullPivLu().solve(V_Tilde));
  };

  // Calculate the F12 contribution to the MP2 energy
  double E_F12 = 0.0;

  // Loop over pairs of electrons
  for (Index i = 0; i < nOcc; i++)
    for (Index j = 0; j < nOcc; j++) {
      E_F12 += PairEnergy(i, j, false);

      // For i < j there are both an alpha-alpha and a beta-beta pair
      if (i < j)
        E_F12 += 2.0 * PairEnergy(i, j, true);
    }

  if (F_out)
    *F_out = F_;
  if (V_out)
    *V_out = V_;
  if (X_out)
    *X_out = X_;
  if (B_out)
    *B_out = B_;

  return {E_MP2, E_F12};
}
//
//
} // namespace SCESS::F12
