#pragma once

namespace SCESS::BasisUtils {
    class F12BasisSet;
}

//Evaluate integrals with f(12) correlation kernels
namespace SCESS {
    using TEITensor = Eigen::Tensor<double, 4>;

    namespace Integrals::F12 {
        using LibIntCGF = std::vector<std::pair<double, double> >;

        LibIntCGF CalculateSlaterContraction(double gamma = 1.0);
        LibIntCGF CalculateSlaterSquaredContraction(const LibIntCGF& slaterContraction);

        /**
         * Evaluate the integral (ij | exp(-gamma * r12) | kl) for a Gaussian-contracted Slater type f(12)
         * @param basisSet Basis set to calculate matrix elements for
         * @param gamma The exponent of the f(12) correlation factor
         * @return A tensor of the required integrals
         */
        TEITensor TwoElectronCGFF12Tensor(const BasisUtils::F12BasisSet& basisSet, const LibIntCGF &cgf);

        /**
         * Evaluate the integral (ij | exp(-gamma * r12) * 1/r12 | kl) for a Gaussian-contracted Slater type f(12)
         * @param basisSet Basis set to calculate matrix elements for
         * @param gamma The exponent of the f(12) correlation factor
         * @return A tensor of the required integrals
         */
        TEITensor TwoElectronCGFF12RepulsionTensor(const BasisUtils::F12BasisSet& basisSet, const LibIntCGF &cgf);


        /**
         * Evaluate the integral (ij | |Del exp(-gamma * r12)|^2 | kl) for a Gaussian-contracted Slater type f(12)
         * @param basisSet
         * @param cgf
         * @return
         */
        TEITensor TwoElectronCGFDelF12SquaredTensor(const BasisUtils::F12BasisSet& basisSet, const LibIntCGF &cgf);

    }
}
