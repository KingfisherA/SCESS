#include "Tests.h"
#include "Tensor/TensorOperations.h"
#include "Tensor/TensorSlice.h"
#include "Tensor/TensorEinsum.h"
#include<iostream>

namespace SCESS::Tests {
    using namespace Eigen;

    template<typename Arr> requires requires(Arr arr)
    {
        { arr.size() };
    }
    void PrintArray(const Arr &arr) {
        std::cout << " [";
        for (int i = 0; i < arr.size(); i++) {
            std::cout << arr[i] << ' ';
        }
        std::cout << "] \n";
    }

    void TestTensorEinsum() {
        enum { alpha, beta, I, J, K, L, M, N, O, P, Q, R, S, T, U, V };

        constexpr std::array<Index, 4> dims{5, 5, 5, 5};
        Tensor<double, 4> t1(dims);
        Tensor<double, 4> t2(dims);
        t1.setRandom();
        t2.setRandom();

        auto op = Einsum<Indices<alpha, M>, Indices<beta, M, I, J>, Indices<alpha, beta, I, J> >(t1,t2);

    }
}
