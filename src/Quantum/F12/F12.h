#pragma once
#include "Quantum/HF/HartreeFock.h"
#include<unsupported/Eigen/CXX11/Tensor>

namespace SCESS {
    struct RHFWavefunction;

    class F12 {
    public:
        using Tensor4 = Eigen::Tensor<double, 4>;

        static double F12_MP2(const Molecule& mol, const RHFWavefunction &rhfWfn);

        static Tensor4 MakeBxyvwTensor();
        static Tensor4 MakeXxyvwTensor();
        static Tensor4 MakeCabxyTensor();
        static Tensor4 MakeVijxyTensor();
    };
}
