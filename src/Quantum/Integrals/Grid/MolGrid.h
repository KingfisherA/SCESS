#pragma once

#include<vector>

#include "AtomicGrid.h"
#include "Quadrature.h"

namespace SCESS {
    class BasisSet;
    class Molecule;
}

namespace SCESS::Quadrature {
    //
    //
    class MolGrid {
    public:
        struct Resolution {
            uint32_t radialPoints;
            LebedevOrder lebedevOrder;
        };

        //Water/STO-3G LDA errors against the converged grid: LOW ~1e-3, MEDIUM ~1e-5, HIGH and ULTRAFINE < 1e-5
        //hartrees
        constexpr static Resolution LOW_RES{30, LebedevOrder::LEBEDEV_74};
        constexpr static Resolution MEDIUM_RES{50, LebedevOrder::LEBEDEV_110};
        constexpr static Resolution HIGH_RES{75, LebedevOrder::LEBEDEV_146};
        constexpr static Resolution ULTRAFINE_RES{100, LebedevOrder::LEBEDEV_194};

        explicit MolGrid(const Molecule *molecule, const BasisSet *basisSet);

        void Create(Resolution resolution);

        //Integration weights of every grid point (atomic grids in order)
        [[nodiscard]] const Eigen::VectorXd &GetWeights() const { return mWeights; }
        [[nodiscard]] Eigen::VectorXd GetDensities() const;
        //Values of the basis functions (BasisSet::GetNoFunctions(), the functions of the integrals) at every grid
        //point, one row per point
        [[nodiscard]] const Eigen::MatrixXd &GetAmplitudes() const { return mAmplitudes; }
        void CorrectDensities();

        //Sets rho = 2 phi^T D phi at every grid point (D is the density of one spin over the basis functions)
        void SetDensity(const Eigen::MatrixXd &D);
        [[nodiscard]] double CalculateDensity() const;

        [[nodiscard]] double GetBeckeWeight(uint32_t i, const Eigen::Vector3d &p0) const;

        std::vector<AtomicGrid> &GetAtomicGrids() { return mAtomicGrids; }
        [[nodiscard]] const std::vector<AtomicGrid> &GetAtomicGrids() const { return mAtomicGrids; }

        [[nodiscard]] size_t GetGridSize() const { return mGridSize; }

    private:
        std::vector<AtomicGrid> mAtomicGrids;

        uint32_t mAngularPoints = 0;
        size_t mGridSize = 0;

        //Basis function values (points x functions) and weights of every grid point
        Eigen::MatrixXd mAmplitudes;
        Eigen::VectorXd mWeights;

        const Molecule *mMolecule = nullptr;
        const BasisSet *mBasisSet = nullptr;
    };
    //
    //
} // SCESS
