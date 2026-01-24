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
} // namespace SCESS::CC
