#include<pch.h>

#include  "Quantum/Integrals/Grid/GridPoint.h"
#include  "Quantum/Molecule.h"
#include "Quantum/Basis/BasisSet.h"

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

    void GridPoint::SetBasisSetAmplitudes(const BasisSet &basisSet) {
        mBasisAmplitudes = VectorXd::Zero(basisSet.GetNoCGFs());

        for (uint32_t i = 0; i < basisSet.GetNoCGFs(); i++)
            mBasisAmplitudes(i) = basisSet.Get(i).Evaluate(mPos);
    }

    void GridPoint::SetBasisSetGradients(const BasisSet &basisSet) {
        mBasisGradients = MatrixXd::Zero(basisSet.GetNoCGFs(), 3);

        for (uint32_t i = 0; i < basisSet.GetNoCGFs(); i++)
            mBasisGradients.row(i) = basisSet.Get(i).EvaluateGradient(mPos);
    }
} // SCESS
