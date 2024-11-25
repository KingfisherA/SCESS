#include "pch.h"

#include "Utilities.h"

#include <iomanip>
#include <iostream>
#include <Eigen/src/Core/util/Meta.h>

namespace SCESS {

    void Utilities::PrintMatrix(const Eigen::MatrixXd &matrix) {
        std::cout << std::fixed << std::setprecision(2);
        for (Eigen::Index row = 0; row < matrix.rows(); ++row) {
            std::cout << "| ";
            for (Eigen::Index col = 0; col < matrix.cols(); ++col) {
                if (matrix(row, col) >= 0.0) {
                    std::cout << ' ';
                }

                std::cout << matrix(row, col) << ' ';
            }
            std::cout << " |\n";
        }
    }

    void Utilities::PrintMatrix(const Eigen::MatrixXd &matrix, std::stringstream &ss) {
        for (Eigen::Index row = 0; row < matrix.rows(); ++row) {
            ss << "| ";
            for (Eigen::Index col = 0; col < matrix.cols(); ++col) {
                ss << matrix(row, col) << ' ';
            }
            ss << " |\n";
        }
    }

    void Utilities::PrintTensorMatrix(const Eigen::Tensor<double, 2> &matrix) {
        std::cout << std::fixed << std::setprecision(2);
        for (Eigen::Index row = 0; row < matrix.dimension(0); ++row) {
            std::cout << "| ";
            for (Eigen::Index col = 0; col < matrix.dimension(1); ++col) {
                std::cout << matrix(row, col) << ' ';
            }
            std::cout << " |\n";
        }
    }
}
