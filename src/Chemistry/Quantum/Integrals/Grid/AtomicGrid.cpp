#include<pch.h>
#include<numbers>

#include "Chemistry/Quantum/Molecule.h"
#include "Chemistry/Quantum/Integrals/Grid/AtomicGrid.h"

#include "Chemistry/Quantum/Integrals/Grid/Quadrature.h"

namespace SCESS {
    using namespace Eigen;
    using std::numbers::pi;


    void AtomicGrid::Create(const Atom &atom, const Molecule &mol, uint32_t angularPoints,
                            uint32_t radialPoints, uint32_t lebedevOffset) {
        mNumAngularPoints = angularPoints;
        mNumRadialPoints = radialPoints;
        mLebedevOffset = lebedevOffset;
        mNumBasisSetFuncs = mol.GetNoOrbitals();

        const double f = pi / static_cast<double>(radialPoints + 1);
        mRadialPoints = VectorXd::Zero(radialPoints);

        // loop over all the radial points
        for (uint32_t p = 1; p <= radialPoints; p++) {
            // calculate weight function of the Gauss-Chebyshev integration
            double w = std::sin(f * static_cast<double>(p));
            w = f * w * w;

            // calculate x point in the interval [-1,1] for Gauss-Chebyshev integration
            double x = std::cos(f * static_cast<double>(p));

            // convert x to r so that the interval is converted from [-1,1] to [0, infinity]
            const double r = (1.0 + x) / (1.0 - x);
            mRadialPoints(p - 1) = r;

            // calculate weight function of the Gauss-Chebyshev integration
            w = w / std::sqrt(1.0 - x * x) *
                2.0 / std::pow(1.0 - x, 2.0);

            // loop over all the angular points
            for (uint32_t a = lebedevOffset; a < angularPoints + lebedevOffset; a++) {
                // create a new grid point given the coordinates on the unit sphere; these
                // coordinates are multiplied by the radius as calculated above
                Vector3d lebedevPoint = Vector3d(Quadrature::lebedevCoefficients[a][0],
                                                 Quadrature::lebedevCoefficients[a][1],
                                                 Quadrature::lebedevCoefficients[a][2]) * r;
                GridPoint &gridPoint =
                        mGridPoints.emplace_back(lebedevPoint + atom.GetPos(), lebedevPoint);

                // set the total integration weight by multiplying the Gauss-Chebyshev weight by the Lebedev weight
                // note that the 'weight' from the Jacobian due to spherical coordinates and the weight
                // due to Becke grid are calculated later
                gridPoint.SetWeight(w * Quadrature::lebedevCoefficients[a][3]);

                // set the amplitudes of the basis functions at this grid point
                gridPoint.SetBasisSetAmplitudes(mol);

                // correct the weights for the Jacobian r**2 and the 4PI from the Lebedev quadrature
                gridPoint.ScaleWeight(gridPoint.GetRelativePosition().squaredNorm() * 4.0 * pi);
            }
        }
    }

    void AtomicGrid::SetDensity(const Eigen::MatrixXd &D) {
        for (GridPoint &pnt: mGridPoints)
            pnt.SetDensity(D);

        mDensityInvalidated = true;
    }

    VectorXd AtomicGrid::GetDensities() const {
        VectorXd densities = VectorXd::Zero(static_cast<Eigen::Index>(mGridPoints.size()));

        for (uint32_t i = 0; i < mGridPoints.size(); i++) {
            densities(i) = mGridPoints[i].GetDensity();
        }

        return densities;
    }

    void AtomicGrid::ScaleDensity(double s) {
        for (GridPoint &pnt: mGridPoints)
            pnt.ScaleDensity(s);

        mDensityInvalidated = true;
    }

    void AtomicGrid::UpdateDensity() const {
        mDensity = 0.0;

        for (const GridPoint &pnt: mGridPoints)
            mDensity += pnt.GetWeight() * pnt.GetDensity();

        mDensityInvalidated = false;
    }

    double AtomicGrid::GetDensity() const {
        if (mDensityInvalidated) UpdateDensity();
        return mDensity;
    }

    MatrixXd AtomicGrid::GetAmplitudes() const {
        MatrixXd amplitudes = MatrixXd::Zero(mNumBasisSetFuncs,
                                             static_cast<Eigen::Index>(mGridPoints.size()));

        for (uint32_t i = 0; i < mGridPoints.size(); i++)
            amplitudes.col(i) = mGridPoints[i].GetBasisAmplitudes();

        return amplitudes;
    }

    VectorXd AtomicGrid::GetWeights() const {
        VectorXd weights = VectorXd(mGridPoints.size());

        for (uint32_t i = 0; i < mGridPoints.size(); i++)
            weights(i) = mGridPoints[i].GetWeight();

        return weights;
    }

    MatrixXd AtomicGrid::GetPositions() const {
        MatrixXd positions = MatrixXd::Zero(
            static_cast<Eigen::Index>(mGridPoints.size()),
            3);

        for (uint32_t i = 0; i < mGridPoints.size(); i++)
            positions.row(i) = mGridPoints[i].GetPosition();

        return positions;
    }

    void AtomicGrid::SetBeckeWeights(const Eigen::VectorXd &weights) {
        for (uint32_t i = 0; i < weights.size(); i++) {
            mGridPoints[i].ScaleWeight(weights(i));
            mGridPoints[i].SetBeckeWeight(weights(i));
        }
    }
} // SCESS
