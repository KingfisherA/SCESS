#include<pch.h>

#include  "Quantum/Integrals/Grid/GridPoint.h"
#include  "Quantum/Molecule.h"

namespace SCESS {
    using namespace Eigen;

    GridPoint::GridPoint(const Eigen::Vector3d &absolutePos,
                         const Eigen::Vector3d &relativePos)
        : mPos(absolutePos), mRelativePos(relativePos) {
    }

    void GridPoint::SetDensity(const Eigen::MatrixXd &D) {
        mDensity = 2.0 * mBasisAmplitudes.dot(D * mBasisAmplitudes);
    }

    void GridPoint::SetGradient(const Eigen::MatrixXd &D) {
        auto x = mBasisGradients.col(0);
        auto y = mBasisGradients.col(1);
        auto z = mBasisGradients.col(2);

        // Apply product rule for the gradient
        double gx = 2.0 * mBasisAmplitudes.transpose().dot(D * x) +
                    2.0 * x.transpose().dot(D * mBasisAmplitudes);
        double gy = 2.0 * mBasisAmplitudes.transpose().dot(D * y) +
                    2.0 * y.transpose().dot(D * mBasisAmplitudes);
        double gz = 2.0 * mBasisAmplitudes.transpose().dot(D * z) +
                    2.0 * z.transpose().dot(D * mBasisAmplitudes);

        mGradient = {gx, gy, gz};
    }

    void GridPoint::SetBasisSetAmplitudes(const Molecule &mol) {
        mBasisAmplitudes = VectorXd::Zero(mol.GetNoOrbitals());

        for (uint32_t i = 0; i < mol.GetNoOrbitals(); i++)
            mBasisAmplitudes(i) = mol.GetOrbital(i).Evaluate(mPos);
    }

    void GridPoint::SetBasisSetGradients(const Molecule &mol) {
        mBasisGradients = MatrixXd::Zero(mol.GetNoOrbitals(), 3);

        for (uint32_t i = 0; i < mol.GetNoOrbitals(); i++)
            mBasisGradients.row(i) = mol.GetOrbital(i).EvaluateGradient(mPos);
    }
} // SCESS
