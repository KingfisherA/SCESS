#pragma once
#include "Quantum/HF/HartreeFock.h"
#include<unsupported/Eigen/CXX11/Tensor>

namespace SCESS {
    namespace Integrals::F12 {
        struct F12Integrals;
    }

    struct RHFWavefunction;
    struct UHFWavefunction;

    namespace F12 {
        using Tensor4 = Eigen::Tensor<double, 4>;

        double CCSD_F12(const UHFWavefunction& refWfn, Integrals::F12::F12Integrals& f12Ints);
    }
}
