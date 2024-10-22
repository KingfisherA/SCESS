#pragma once

#include "Chemistry/Quantum/Molecule.h"
#include "Chemistry/Quantum/Integrals/Integrals.h"

namespace SCESS {
    class MolecularGrid;
    class Molecule;

    class DensityFunctional {
    public:
        static double DFT(const Molecule &mol);

        static Eigen::MatrixXd CreateColoumbMatrixTEI(const TEITensor &tei,
                                                   const Eigen::MatrixXd &Density);

        static Eigen::MatrixXd CreateDensityMatrix(const Eigen::MatrixXd& C, uint32_t Nelecs);

        static void CalculateXCMatrix(const Molecule& mol, const MolecularGrid& grid,
            Eigen::MatrixXd& ExchangeCorr, double& E_XC);
    };
}
