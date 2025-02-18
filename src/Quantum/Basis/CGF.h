#pragma once
#include<Eigen/Core>
#include<boost/math/special_functions/factorials.hpp>
#include<vector>

namespace SCESS {
    class BasisSet;
    class Molecule;
}

namespace SCESS {
    //
    //
    //Represents a primitive gaussian function with three Cartesian harmonics
    class PGF {
    public:
        PGF(double exponent, double coeff, const Eigen::Vector3d &position = {0., 0., 0.},
            const Eigen::Vector3i &angular = {0, 0, 0}) : mPosition(position), mExponent(exponent),
                                                          mAngular(angular), mCoeff(coeff), mNorm(1.0) {
        }

        void CalculateNormalisationConstant();

        const Eigen::Vector3i &GetAngularMomentum() const {
            return mAngular;
        }

    public:
        Eigen::Vector3d mPosition;
        double mExponent;
        Eigen::Vector3i mAngular;
        double mCoeff;
        double mNorm;
    };

    //
    //Represents a collection of Gaussian functions with the same position
    //which form the basis for a single atom
    class CGF {
    public:
        CGF() : mPosition({0.0, 0.0, 0.0}) {
        }

        CGF(const CGF &col) = default;

        explicit CGF(const Eigen::Vector3d &position) : mPosition(position) {
        }

        [[nodiscard]] double Evaluate(const Eigen::Vector3d &pos) const;
        [[nodiscard]] Eigen::Vector3d EvaluateGradient(const Eigen::Vector3d &pos) const;

        PGF &Get(uint32_t i) {
            return mGaussians[i];
        }

        const PGF &Get(uint32_t i) const {
            return mGaussians[i];
        }

        void Add(const Eigen::Vector3i &angMom, const PGF &gaussianIn);

        void UpdatePosition(const Eigen::Vector3d &pos) {
            mPosition = pos;
            for (auto &g: mGaussians) {
                g.mPosition = pos;
            }
        }

        const Eigen::Vector3d &GetPos() const {
            return mPosition;
        }

        size_t Size() const {
            return mGaussians.size();
        }

        std::vector<PGF> &GetGaussians() {
            return mGaussians;
        }

        const std::vector<PGF> &GetGaussians() const {
            return mGaussians;
        }

    private:
        Eigen::Vector3d mPosition;
        std::vector<PGF> mGaussians;
    };
    //
    template<typename T>
    std::function<T(T, T, T)> GenerateMolecularWavefunction(const Eigen::VectorXd &coeffsIn,
                                                            const BasisSet &basisSet);
    //
}
