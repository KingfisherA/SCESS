#include<pch.h>

#include "MolecularGrid.h"

#include "Chemistry/Quantum/Molecule.h"
#include "Chemistry/Quantum/Integrals/Grid/Quadrature.h"

namespace SCESS {
    using namespace Eigen;

    MolecularGrid::MolecularGrid(const Molecule *mol): mMolecule(mol) {
    }

    void MolecularGrid::Create(Resolution res) {
        mAngularPoints = Quadrature::nLebedevPoints[(int)res.lebedevOrder];

        // The Lebedev coefficients and radial points are stored in a large matrix.
        // Given a particular order for the integration, we have to start reading
        // this matrix from a specific position. Here, we calculate that position.
        uint32_t lebedevOffset = 0;
        for (uint32_t i = 0; i < (int)res.lebedevOrder; i++)
            lebedevOffset += Quadrature::nLebedevPoints[i];

        // generate for each atom an atomic grid
        for (uint32_t i = 0; i < mMolecule->GetNoAtoms(); i++) {
            AtomicGrid &grid =
                    mAtomicGrids.emplace_back();
            grid.Create(mMolecule->GetAtom(i),
                        *mMolecule,
                        mAngularPoints,
                        res.radialPoints,
                        lebedevOffset);

            mGridSize += grid.GetNoOfGridPoints();
        }

        // calculate the Becke weights for the atomic grids
        for (uint32_t i = 0; i < mMolecule->GetNoAtoms(); i++) {
            MatrixXd positions = mAtomicGrids[i].GetPositions();
            VectorXd beckeWeightCoeffs = VectorXd::Zero(positions.rows());

            // loop over positions
            for (uint32_t j = 0; j < positions.rows(); j++) {
                double denom = 0.0;
                double nom = 1.0;

                // loop over all atoms to get Pn(r)
                for (uint32_t k = 0; k < mMolecule->GetNoAtoms(); k++) {
                    double term = GetBeckeWeight(k, positions.row(j)); // obtain single Pn(r) term
                    denom += term;
                    if (i == k) {
                        nom = term;
                    }
                }

                // set weight from cell function: wn(r) = Pn(r) / SUM_m Pm(r) (eq. 22)
                beckeWeightCoeffs(j) = nom / denom;
            }

            mAtomicGrids[i].SetBeckeWeights(beckeWeightCoeffs);
        }
    }

    VectorXd MolecularGrid::GetWeights() const {
        VectorXd weights = VectorXd::Zero(static_cast<Eigen::Index>(mGridSize));
        uint32_t idx = 0;

        for (const AtomicGrid &grid: mAtomicGrids) {
            VectorXd atomicWeights = grid.GetWeights();
            for (const double &weight: atomicWeights)
                weights(idx++) = weight;
        }

        return weights;
    }

    VectorXd MolecularGrid::GetDensities() const {
        VectorXd densities = VectorXd::Zero(static_cast<Eigen::Index>(mGridSize));
        uint32_t idx = 0;

        for (const AtomicGrid &grid: mAtomicGrids) {
            VectorXd atomicDensities = grid.GetDensities();
            for (const double &density: atomicDensities)
                densities(idx++) = density;
        }

        return densities;
    }

    MatrixXd MolecularGrid::GetAmplitudes() const {
        MatrixXd amplitudes = MatrixXd::Zero(
            static_cast<Eigen::Index>(mMolecule->GetNoOrbitals()),
            static_cast<Eigen::Index>(mGridSize));
        uint32_t idx = 0;

        for (const AtomicGrid &grid: mAtomicGrids) {
            MatrixXd atomicAmplitudes = grid.GetAmplitudes();

            for (uint32_t k = 0; k < mMolecule->GetNoOrbitals(); k++)
                for (uint32_t j = 0; j < grid.GetNoOfGridPoints(); j++)
                    amplitudes(k, idx + j) = atomicAmplitudes(k, j);

            idx += grid.GetNoOfGridPoints();
        }

        return amplitudes;
    }

    void MolecularGrid::CorrectDensities() {
        double sum = 0.0;
        double charge = 0.0;

        for (uint32_t i = 0; i < mAtomicGrids.size(); i++) {
            sum += mAtomicGrids[i].GetDensity();
            charge += mMolecule->GetAtom(i).GetCharge();
        }

        double correctionFactor = charge / sum;

        for (AtomicGrid &grid: mAtomicGrids)
            grid.ScaleDensity(correctionFactor);
    }

    void MolecularGrid::SetDensity(const Eigen::MatrixXd &D) {
        for (AtomicGrid &grid: mAtomicGrids)
            grid.SetDensity(D);
    }

    double MolecularGrid::CalculateDensity() const {
        double sum = 0.0;

        for (const AtomicGrid &grid: mAtomicGrids)
            sum += grid.GetDensity();

        return sum;
    }

    double MolecularGrid::GetBeckeWeight(uint32_t i, const Eigen::Vector3d &p0) const {
        auto fK = [](uint32_t k, double mu) {
            for (uint32_t i = 0; i < k; i++)
                mu = 3.0 / 2.0 * mu - 0.5 * std::pow(mu, 3.0);

            return mu;
        };

        auto beckeCutoff = [&fK](double mu) {
            return 0.5 * (1.0 - fK(3, mu));
        };

        double wprod = 1.0;

        // get position of central atom in fuzzy cell
        const Vector3d &p1 = mMolecule->GetAtom(i).GetPos();

        // loop over all other atoms
        for (uint32_t j = 0; j < mMolecule->GetNoAtoms(); j++) {
            if (i == j) {
                continue;
            }

            // get position of other atom
            const Vector3d &p2 = mMolecule->GetAtom(j).GetPos();

            // calculate mu (eq. 11 in A.D. Becke J.Chem.Phys. 88, 2547)
            double mu = ((p0 - p1).norm()
                         - (p0 - p2).norm()) /
                        (p2 - p1).norm();

            // (eq. 13 in A.D. Becke J.Chem.Phys. 88, 2547)
            wprod *= beckeCutoff(mu);
        }

        return wprod;
    }
} // SCESS
