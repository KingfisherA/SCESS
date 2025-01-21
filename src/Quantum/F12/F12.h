#pragma once
#include "Quantum/HF/HartreeFock.h"
#include <unsupported/Eigen/CXX11/Tensor>

namespace SCESS {
    namespace Integrals {
        struct MolInts;
    };

    namespace Integrals::F12 {
        struct F12Integrals;
    }

    struct RHFWavefunction;
    struct UHFWavefunction;

    namespace F12 {
        using Tensor4 = Eigen::Tensor<double, 4>;

        namespace CCSD {
            double CCSD_F12(const UHFWavefunction &refWfn, Integrals::F12::F12Integrals &f12Ints);
        }

        namespace MP2 {

            struct MP2F12Result {
                double E_MP2;
                double E_F12_MP2;
            };

            MP2F12Result MP2_F12(const UHFWavefunction &refWfn, const Integrals::MolInts& configIntegrals, const Integrals::F12::F12Integrals &f12Ints);
        }
    }
}
