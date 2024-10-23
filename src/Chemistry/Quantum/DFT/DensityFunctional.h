#pragma once

#include "GraphicsEngine/Utilities/Memory.h"

#include "Chemistry/Quantum/Molecule.h"
#include "Chemistry/Quantum/Integrals/Integrals.h"
#include "Chemistry/Quantum/Integrals/Grid/MolecularGrid.h"

namespace SCESS {
    class Molecule;

    //A structure representing a Kohn-Sham wavefunction
    struct KSWavefunction {
    public:
        //The Kohn-Sham Fock matrix
        Eigen::MatrixXd F;
        //The converged KS orbital coefficients
        Eigen::MatrixXd C;
        //The converged density matrix
        Eigen::MatrixXd Density;

        //Total electronic energy
        double ElectronicE = 0.0;

        Eigen::VectorXd OrbitalEnergies;
        Eigen::VectorXi Occupation;
        std::optional<ConfigIntegrals> Integrals;
    };

    class DensityFunctional {
    public:
        static uRef<KSWavefunction> KohnShamDFT(const Molecule &mol, MolecularGrid::Resolution gridRes,
            std::stringstream* console = nullptr, std::string_view basisSetName = "");

        static Eigen::MatrixXd CreateColoumbMatrixTEI(const TEITensor &tei,
                                                      const Eigen::MatrixXd &Density);

        static Eigen::MatrixXd CreateDensityMatrix(const Eigen::MatrixXd &C, uint32_t Nelecs);

        static void CalculateXCMatrix(const Molecule &mol, const MolecularGrid &grid,
                                      Eigen::MatrixXd &ExchangeCorr, double &E_XC);
    };
}
