#pragma once
#include "Quantum/HF/HartreeFock.h"
#include <unsupported/Eigen/CXX11/Tensor>

//#include "Tensor/TensorOperations.h"

namespace SCESS {
    namespace Integrals {
        struct MolInts;
    };

    namespace Integrals::F12 {
        struct F12Integrals;
    }

    //using TensorSliceIndices = std::unordered_map<int, std::pair<Idx, Idx> >;

    struct RHFWavefunction;
    struct UHFWavefunction;

    namespace F12 {
        using Tensor2 = Eigen::Tensor<double, 2>;
        using Tensor4 = Eigen::Tensor<double, 4>;

        double EvaluateSinglesCorrectionEnergy(const Tensor2 &F,
                                               Tensor2 *t1Amps = nullptr);

        namespace CCSD {
            double CCSD_F12(const UHFWavefunction &refWfn, Integrals::F12::F12Integrals &f12Ints);
        }

        namespace MP2 {
            struct MP2F12Result {
                double E_MP2;
                double E_F12_MP2;
            };

            MP2F12Result MP2_F12_Optimise(const UHFWavefunction &refWfn,
                                          const Integrals::MolInts &configIntegrals,
                                          const Integrals::F12::F12Integrals &f12Ints);

            MP2F12Result MP2_F12_Fixed(const UHFWavefunction &refWfn,
                                       const Integrals::MolInts &configIntegrals,
                                       const Integrals::F12::F12Integrals &f12Ints,
                                       Tensor2*   F_out = nullptr,
                                       TEITensor* V_out = nullptr,
                                       TEITensor* X_out = nullptr,
                                       TEITensor* B_out = nullptr);

            double CalculateMP2_F12_FixedEnergyElectronPair(long electronI, long electronJ,
                                                const Tensor2 &F,
                                                const TEITensor &V,
                                                const TEITensor &X,
                                                const TEITensor &B);
        }

        namespace QMC {
            /**
             * Generate the W intermediate for CCQMC F12
             * @param f12Int The <pq | f(12) | rs> integrals
             * @return The W tensors
             */
            Tensor4 EvaluateW(const Tensor4 &f12Int);

            double CCQMC_F12(const UHFWavefunction &refWfn, const Integrals::MolInts &mints,
                             const Integrals::F12::F12Integrals &f12Ints);
        }
    }
}
