#pragma once
#include<vector>
#include<Eigen/Core>
#include "GridPoint.h"

namespace SCESS {
    class Atom;
}

//
//
namespace SCESS::Quadrature {
    //
    //
    //Bragg radii for atoms
    static constexpr double BRAGG_RADII[] = {
        0.35,
        0.28,
        1.28,
        0.96,
        0.84,
        0.70,
        0.71,
        0.66,
        0.57,
        0.58,
        1.66,
        1.41,
        1.21,
        1.11,
        1.07,
        1.05,
        1.02,
        1.06,
        2.03,
        1.76
    };
    //
    //Omegab97-v
    class AtomicGrid {
    public:
        //Assign an atom to this grid
        AtomicGrid() = default;

        void Create(const Atom &atom,
                    const BasisSet &basisSet,
                    uint32_t angularPoints,
                    uint32_t radialPoints,
                    uint32_t lebedevOffset);

        void SetDensity(const Eigen::MatrixXd &D);

        void ScaleDensity(double s);
        void UpdateDensity() const;
        [[nodiscard]] double GetDensity() const;

        [[nodiscard]] Eigen::VectorXd GetDensities() const;
        [[nodiscard]] Eigen::MatrixXd GetAmplitudes() const;
        [[nodiscard]] Eigen::VectorXd GetWeights() const;
        [[nodiscard]] Eigen::MatrixXd GetPositions() const;

        [[nodiscard]] size_t GetNoOfGridPoints() const { return mGridPoints.size(); }

        void SetBeckeWeights(const Eigen::VectorXd &weights);

    private:
        std::vector<GridPoint> mGridPoints;

        uint32_t mNumAngularPoints = 0,
                mNumRadialPoints = 0,
                mLebedevOffset = 0,
                mNumBasisSetFuncs = 0;

        mutable bool mDensityInvalidated = true;
        mutable double mDensity = 0.0;

        Eigen::VectorXd mRadialPoints;
    };
    //
    //
} // SCESS
