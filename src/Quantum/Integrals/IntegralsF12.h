#pragma once

#include<libint2/atom.h>


//Evaluate integrals with f(12) correlation kernels
namespace SCESS {
    using TEITensor = Eigen::Tensor<double, 4>;
    class F12BasisSet;

    namespace Integrals {
        struct ConfigInts;
    }

    namespace Integrals::F12 {
        using LibIntCGF = std::vector<std::pair<double, double> >;
        using AtomListLI = std::vector<libint2::Atom>;

        struct F12Integrals {
        public:
            //eri = ( ij | 1/r12 | kl )
            TEITensor eri;
            //f12 = ( ij | f(12) | kl)
            TEITensor f12;
            //f12Sq = ( ij | f(12)^2 | kl)
            TEITensor f12Sq;
            //f12DelSq = ( ij | |Del f(12)|^2 | kl)
            TEITensor f12DelSq;
            //f12Repulsion = ( ij | f(12)/r12 | kl)
            TEITensor f12Repulsion;
        };

        ConfigInts CalculateConfiguration(const F12BasisSet &basisSet, const AtomListLI &atoms);

        Eigen::MatrixXd CalculateOverlapMatrix(const F12BasisSet &basisSet);
        Eigen::MatrixXd CalculateNuclearMatrix(const F12BasisSet &basisSet, const AtomListLI &atoms);
        Eigen::MatrixXd CalculateKineticMatrix(const F12BasisSet &basisSet);
        TEITensor TwoElectronRepulsionTensor(const F12BasisSet &basisSet);

        //Chemist's notation
        F12Integrals EvaluateF12Integrals(const F12BasisSet &basisSet, double gamma = 1.0);
        /**
         * Import an integral file from MRCC
         * @param fileName The file path of integral dump to read from
         * @param nbf Number of basis set functions in the OBS and CABS
         * @return The imported F12 integrals
         */
        F12Integrals ImportF12IntegralsMRCC(const std::string &fileName, uint32_t nbf);

        LibIntCGF CalculateSlaterContraction(double gamma = 1.0);
        LibIntCGF CalculateSlaterSquaredContraction(const LibIntCGF &slaterContraction);

        /**
         * Evaluate the integral (ij | exp(-gamma * r12) | kl) for a Gaussian-contracted Slater type f(12)
         * @param basisSet Basis set to calculate matrix elements for
         * @param gamma The exponent of the f(12) correlation factor
         * @return A tensor of the required integrals
         */
        TEITensor TwoElectronCGFF12Tensor(const F12BasisSet &basisSet, const LibIntCGF &cgf);

        /**
         * Evaluate the integral (ij | exp(-gamma * r12) * 1/r12 | kl) for a Gaussian-contracted Slater type f(12)
         * @param basisSet Basis set to calculate matrix elements for
         * @param gamma The exponent of the f(12) correlation factor
         * @return A tensor of the required integrals
         */
        TEITensor TwoElectronCGFF12RepulsionTensor(const F12BasisSet &basisSet, const LibIntCGF &cgf);

        /**
         * Evaluate the integral (ij | |Del exp(-gamma * r12)|^2 | kl) for a Gaussian-contracted Slater type f(12)
         * @param basisSet
         * @param cgf
         * @return
         */
        TEITensor TwoElectronCGFDelF12SquaredTensor(const F12BasisSet &basisSet, const LibIntCGF &cgf);
    }
}
