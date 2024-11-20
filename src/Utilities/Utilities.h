#pragma once
#include <Eigen/Core>
#include<unsupported/Eigen/CXX11/Tensor>

namespace SCESS {
    class Utilities {
    public:
        static void PrintMatrix(const Eigen::MatrixXd &matrix);
        static void PrintMatrix(const Eigen::MatrixXd &matrix, std::stringstream &ss);
        static void PrintTensorMatrix(const Eigen::Tensor<double, 2> &matrix);

    private:
        template<typename T, int rank>
        static void PrintTensorRecursive(const Eigen::Tensor<T, rank> &tensor, uint32_t dimensionIdx = 0,
                                         std::array<uint32_t, rank> *indices = nullptr) {
            for (uint32_t i = 0; i < tensor.dimension(dimensionIdx); i++) {
                if (dimensionIdx + 1 < rank) {
                    (*indices)[dimensionIdx] = i;
                    PrintTensorRecursive<T, rank>(tensor, dimensionIdx + 1, indices);
                    continue;
                }
                std::cout << "  [";
                (*indices)[dimensionIdx] = i;
                for (uint32_t j = 0; j < rank; j++) {
                    std::cout << (*indices)[j];
                    if (j + 1 < rank)
                        std::cout << ", ";
                }
                std::cout << "] = ";

                std::cout << tensor(*indices) << '\n';
            }
        }

    public:
        template<typename T, int rank>
        static void PrintTensor(const Eigen::Tensor<T, rank> &tensor) {
            //std::cout << "Rank " << rank << " Tensor:\n";
            std::array<uint32_t, rank> indicesArray;

            PrintTensorRecursive<T, rank>(tensor, 0, &indicesArray);
        }
    };
}
