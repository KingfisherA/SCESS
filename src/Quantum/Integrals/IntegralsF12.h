#pragma once

//Evaluate integrals with f(12) correlation kernels
namespace SCESS {
    using TEITensor = Eigen::Tensor<double, 4>;

    namespace Integrals::F12 {
        /**
         * Evaluate the integral (ij | r12 * exp(-gamma * r12) | kl) for a Slater type f(12)
         * @param mol Molecule to calculate matrix elements for
         * @param gamma The exponent of the f(12) correlation factor
         * @return A tensor of the required integrals
         */
        TEITensor TwoElectronSlaterF12Tensor(const Molecule &mol, double gamma);

        /**
         * Evaluate the integral (ij | r12 * exp(-gamma * r12) * 1/r12 | kl) for a Slater type f(12)
         * @param mol Molecule to calculate matrix elements for
         * @param gamma The exponent of the f(12) correlation factor
         * @return A tensor of the required integrals
         */
        TEITensor TwoElectronSlaterF12RepulsionTensor(const Molecule &mol, double gamma);
    }
}
