#pragma once
#include <Eigen/Core>

namespace SCESS {
    class Atom;
    class Molecule;
    class BasisSet;
}

//
//
namespace SCESS::Quadrature {
    //
    //
    class GridPoint {
    public:
        GridPoint(const Eigen::Vector3d &absolutePos,
                  const Eigen::Vector3d &relativePos);

        //Set density from the density matrix
        void SetDensity(const Eigen::MatrixXd &D);
        void SetGradient(const Eigen::MatrixXd &D);
        void SetBasisSetAmplitudes(const BasisSet &basisSet);
        void SetBasisSetGradients(const BasisSet &basisSet);

        void SetWeight(double w) {
            mWeight = w;
        }

        void ScaleWeight(double s) {
            mWeight *= s;
        }

        [[nodiscard]] const double &GetWeight() const {
            return mWeight;
        }

        void SetBeckeWeight(double w) {
            mBeckeWeight = w;
        }

        void ScaleBeckeWeight(double s) {
            mBeckeWeight *= s;
        }

        [[nodiscard]] const double &GetBeckeWeight() const {
            return mBeckeWeight;
        }

        void ScaleDensity(double s) {
            mDensity *= s;
        }

        void SetDensityValue(double density) {
            mDensity = density;
        }

        [[nodiscard]] const double &GetDensity() const {
            return mDensity;
        }

        [[nodiscard]] const Eigen::Vector3d &GetPosition() const {
            return mPos;
        }

        [[nodiscard]] const Eigen::Vector3d &GetRelativePosition() const {
            return mRelativePos;
        }

        [[nodiscard]] Eigen::VectorXd GetAtomPos() const {
            return mPos - mRelativePos;
        }

        [[nodiscard]] const Eigen::Vector3d &GetGradient() const {
            return mGradient;
        }

        [[nodiscard]] const Eigen::VectorXd &GetBasisAmplitudes() const {
            return mBasisAmplitudes;
        }

        [[nodiscard]] const Eigen::MatrixXd &GetBasisGradients() const {
            return mBasisGradients;
        }

    private:
        double mDensity = 0.0;
        double mWeight = 0.0;
        double mBeckeWeight = 0.0;
        Eigen::Vector3d mPos;
        //Position relative to the atom this point is attached to
        Eigen::Vector3d mRelativePos;
        Eigen::Vector3d mGradient;
        Eigen::VectorXd mBasisAmplitudes;
        Eigen::MatrixXd mBasisGradients;
    };
    //
    //
} // SCESS
