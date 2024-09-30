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
                std::cout << matrix(row, col) << ' ';
            }
            std::cout << " |\n";
        }
    }
}
