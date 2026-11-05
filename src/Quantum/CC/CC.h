#pragma once
#include <Eigen/Core>
#include <unsupported/Eigen/CXX11/Tensor>
//
#include "Quantum/HF/HF.h"
#include "Quantum/Integrals/Integrals.h"

//
//
namespace DH {
class Thread;
}

namespace SCESS::HF {
struct HF_Wfn;
}

namespace SCESS::Integrals {
struct MolInts;
}

namespace SCESS::CC {
//
using HFWfn = HF::HF_Wfn;
//
//
struct CCSD_Wfn {
public:
  double E_CCSD = 0.0;
  double E_T = 0.0;
  Eigen::Tensor<double, 2> t1;
  Eigen::Tensor<double, 4> t2;
  Eigen::Tensor<double, 4> mints;
};
//
//
using T1Tensor = Eigen::Tensor<double, 2>;
using T2Tensor = Eigen::Tensor<double, 4>;
using FTensor = Eigen::Tensor<double, 2>;
using WTensor = Eigen::Tensor<double, 4>;
using TEITensor = Eigen::Tensor<double, 4>;
using TauTensor = Eigen::Tensor<double, 4>;
using FockTensor = Eigen::Tensor<double, 2>;
//
/**
 * An implementation to solve the coupled cluster equations
 * to single and double excitations for a closed shell molecule starting from a
 * RHF determinant
 * @param rhfResult The RHF determinant
 */
CCSD_Wfn CCSD(const HFWfn &rhfResult, const Integrals::MolInts &mints,
              const DH::Thread *thread = nullptr,
              std::stringstream *console = nullptr);
//
/**
 * An implementation to solve the coupled cluster equations
 * to single and double excitations for a closed shell molecule starting from a
 * RHF determinant with DIIS acceleration.
 * @param rhfWfn The RHF determinant
 */
CCSD_Wfn CCSD_DIIS(const HFWfn &rhfWfn, const Integrals::MolInts &mints,
                   const DH::Thread *thread = nullptr,
                   std::stringstream *console = nullptr);
//
CCSD_Wfn CCSD_T(const HFWfn &rhfWfn, const Integrals::MolInts &mints,
                const DH::Thread *thread = nullptr,
                std::stringstream *console = nullptr);
//
double CalculateTriplesCorrection(const CCSD_Wfn &ccsd, const FockTensor &F,
                                  const TEITensor &eri, Eigen::Index nOcc,
                                  Eigen::Index nVirt,
                                  const DH::Thread *thread = nullptr);
//
FockTensor TransformFockToMO(const HFWfn &rhfWfn);
void CreateDenominatorTensors(const FockTensor &F, const Eigen::Index nOcc,
                              const Eigen::Index nVirt,
                              Eigen::Tensor<double, 2> &Dia,
                              Eigen::Tensor<double, 4> &Dijab);
T2Tensor CreateInitialT2Tensor(const Eigen::Index nOcc,
                               const Eigen::Index nVirt, const TEITensor &moTEI,
                               const Eigen::Tensor<double, 4> &Dijab);
void IterateCCSD(const TEITensor &moTEI, const FockTensor &F,
                 const Eigen::Index nOcc, const Eigen::Index nVirt,
                 const Eigen::Tensor<double, 2> &Dia,
                 const Eigen::Tensor<double, 4> &Dijab, T1Tensor &t1,
                 T2Tensor &t2, const DH::Thread *thread = nullptr);
//
double CalculateEnergyCCSD(const FockTensor &F, const Eigen::Index nOcc,
                           const Eigen::Index nVirt, const TEITensor &moTEI,
                           const T1Tensor &t1, const T2Tensor &t2);
//
//
//-- Closed-shell (spin-adapted, spatial orbital) CCSD --
//
// The amplitudes are those of the spatial orbitals, t1 indexed (i, a) and t2
// indexed (i, j, a, b) with t_ij^ab = t_ji^ba (the opposite-spin amplitudes of
// the spin-orbital formulation). The integrals and amplitudes are 16 times
// smaller than in the spin-orbital CCSD above.
struct RCCSD_Wfn {
  double E_CCSD = 0.0;
  // Energy of the coupling between the conventional doubles and the explicitly
  // correlated geminals (CCSD-F12 methods only, already included in E_CCSD)
  double E_F12_Coupling = 0.0;
  bool converged = false;
  T1Tensor t1;
  T2Tensor t2;
};
//
// Couples explicitly correlated geminals into the closed-shell CCSD equations.
// Both are indexed (i, j, a, b) over the occupied and virtual spatial orbitals:
//   residual: added to the driving <ij|ab> term of the doubles residual
//   energy:   the energy contribution sum_ijab (2 e_ij^ab - e_ij^ba) t_ij^ab
struct RCCSD_F12Coupling {
  Eigen::Tensor<double, 4> residual;
  Eigen::Tensor<double, 4> energy;
};
//
/**
 * Solves the closed-shell CCSD equations (Hirata et al., J. Chem. Phys. 120,
 * 2581 (2004), as in PySCF's rccsd) with DIIS acceleration in canonical
 * orbitals.
 * @param eri <pq|rs> over the spatial molecular orbitals (physicists' notation),
 *            ordered occupied then virtual
 * @param orbitalEnergies Canonical orbital energies of the same orbitals
 * @param nOcc Number of doubly occupied orbitals
 * @param f12 Optional coupling of the doubles to F12 geminals
 */
RCCSD_Wfn SolveRCCSD_DIIS(const TEITensor &eri,
                          const Eigen::VectorXd &orbitalEnergies,
                          Eigen::Index nOcc,
                          const RCCSD_F12Coupling *f12 = nullptr,
                          const DH::Thread *thread = nullptr,
                          std::stringstream *console = nullptr);
//
// Closed-shell CCSD with DIIS from a RHF determinant using SolveRCCSD_DIIS
RCCSD_Wfn RCCSD_DIIS(const HFWfn &rhfWfn, const Integrals::MolInts &mints,
                     const DH::Thread *thread = nullptr,
                     std::stringstream *console = nullptr);
} // namespace SCESS::CC
