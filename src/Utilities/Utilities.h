#pragma once
#include <Eigen/Core>
#include<unsupported/Eigen/CXX11/Tensor>

namespace SCESS {
    class Utilities {
    public:
        static void PrintMatrix(const Eigen::MatrixXd& matrix);
        static void PrintMatrix(const Eigen::MatrixXd& matrix, std::stringstream& ss);
        static void PrintTensorMatrix(const Eigen::Tensor<double,2>& matrix);
    };
}

