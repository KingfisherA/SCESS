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

        constexpr static Resolution LOW_RES{10, LebedevOrder::LEBEDEV_50};
        constexpr static Resolution MEDIUM_RES{15, LebedevOrder::LEBEDEV_74};
        constexpr static Resolution HIGH_RES{20, LebedevOrder::LEBEDEV_110};
        constexpr static Resolution ULTRAFINE_RES{30, LebedevOrder::LEBEDEV_194};

        explicit MolGrid(const Molecule *molecule, const BasisSet *basisSet);

        void Create(Resolution resolution);

        [[nodiscard]] Eigen::VectorXd GetWeights() const;
        [[nodiscard]] Eigen::VectorXd GetDensities() const;
        [[nodiscard]] Eigen::MatrixXd GetAmplitudes() const;
        void CorrectDensities();

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

        const Molecule *mMolecule = nullptr;
        const BasisSet *mBasisSet = nullptr;
    };
    //
    //
} // SCESS
