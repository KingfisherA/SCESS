#pragma once

#include<vector>

#include "AtomicGrid.h"

namespace SCESS {
    class Molecule;

    class MolecularGrid {
    public:
        explicit MolecularGrid(const Molecule *mol);

        void Create(uint32_t radialPoints, uint32_t lebedevOrder);

        [[nodiscard]] Eigen::VectorXd GetWeights() const;
        [[nodiscard]] Eigen::VectorXd GetDensities() const;
        [[nodiscard]] Eigen::MatrixXd GetAmplitudes() const;
        void CorrectDensities();

        void SetDensity(const Eigen::MatrixXd &D);
        [[nodiscard]] double CalculateDensity() const;

        [[nodiscard]] double GetBeckeWeight(uint32_t i, const Eigen::Vector3d &p0) const;

    private:
        std::vector<AtomicGrid> mAtomicGrids;

        uint32_t mAngularPoints = 0;
        size_t mGridSize = 0;

        const Molecule *mMolecule = nullptr;
    };
} // SCESS
