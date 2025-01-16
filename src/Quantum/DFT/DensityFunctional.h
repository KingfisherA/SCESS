#pragma once

#include "GraphicsEngine/Utilities/Memory.h"

#include  "Quantum/Molecule.h"
#include  "Quantum/Integrals/Integrals.h"
#include  "Quantum/Integrals/Grid/MolecularGrid.h"

namespace GEng {
    class Thread;
}

namespace SCESS::Integrals {
    struct ConfigInts;
}

namespace SCESS {
    class Molecule;
    class BasisSet;

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

        Eigen::VectorXd Energies;
        Eigen::VectorXi Occupation;
        std::optional<Integrals::ConfigInts> Integrals;
    };

    class DensityFunctional {
    public:
        static uRef<KSWavefunction> KohnShamDFT(const Molecule &mol,
                                                const BasisSet &basisSet,
                                                MolecularGrid::Resolution gridRes,
                                                const GEng::Thread *thread = nullptr,
                                                std::stringstream *console = nullptr);

        static Eigen::MatrixXd CreateColoumbMatrixTEI(const TEITensor &tei,
                                                      const Eigen::MatrixXd &Density);
        static Eigen::MatrixXd CreateDensityMatrix(const Eigen::MatrixXd &C, uint32_t Nelecs);
        static void CalculateXCMatrix(const MolecularGrid &grid,
                                      Eigen::MatrixXd &ExchangeCorr,
                                      double &E_XC);
    };
}
