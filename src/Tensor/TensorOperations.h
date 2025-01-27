#pragma once

#include <array>
#include <unsupported/Eigen/CXX11/Tensor>

#define Slice2 slice<std::array<Eigen::Index, 2>, std::array<Eigen::Index, 2>>
#define Slice4 slice<std::array<Eigen::Index, 4>, std::array<Eigen::Index, 4>>

#define Shuffle2 shuffle<std::array<Eigen::Index, 2>>
#define Shuffle4 shuffle<std::array<Eigen::Index, 4>>

#define Contract1(other, dimensions) contract(other, std::array<Eigen::IndexPair<Eigen::Index>,1>dimensions)
#define Contract2(other, dimensions) contract(other, std::array<Eigen::IndexPair<Eigen::Index>,2>dimensions)
#define Contract3(other, dimensions) contract(other, std::array<Eigen::IndexPair<Eigen::Index>,2>dimensions)
#define Contract4(other, dimensions) contract(other, std::array<Eigen::IndexPair<Eigen::Index>,4>dimensions)

#define Reduce1 sum<std::array<Eigen::Index, 1>>
#define Reduce2 sum<std::array<Eigen::Index, 2>>
#define Reduce3 sum<std::array<Eigen::Index, 3>>
#define Reduce4 sum<std::array<Eigen::Index, 4>>

namespace SCESS {
    using Idx = Eigen::Index;
    using IdxP = Eigen::IndexPair<Idx>;

    template<Idx ... T>
    struct Indices {
        static constexpr Idx N = sizeof...(T);
        static constexpr std::array<Idx, N> values = {T...};
    };

    //
    //
    inline Eigen::Tensor<double, 2> MatrixToTensor(const Eigen::MatrixXd &matrix) {
        Eigen::Tensor<double, 2> tensor(matrix.rows(), matrix.cols());
        memcpy(tensor.data(), matrix.data(), matrix.size() * sizeof(double));
        return tensor;
    }

    //
    inline Eigen::MatrixXd TensorToMatrix(const Eigen::Tensor<double, 2> &tensor) {
        Eigen::MatrixXd matrix(tensor.dimension(0), tensor.dimension(1));
        memcpy(matrix.data(), tensor.data(), tensor.size() * sizeof(double));
        return matrix;
    }
}
