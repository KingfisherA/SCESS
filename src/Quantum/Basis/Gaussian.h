#pragma once
#include<Eigen/Core>
#include<numbers>
#include<unordered_map>
#include<boost/math/special_functions/factorials.hpp>
#include<vector>

namespace SCESS {
    class BasisSet;
}

namespace SCESS {
    class Molecule;

    //Represents a single gaussian function with three polynomial components
    class Gaussian {
    public:
        Gaussian(double exponent, double coeff, const Eigen::Vector3d &position = {0., 0., 0.},
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


    //Represents a collection of Gaussian functions with the same position
    //which form the basis for a single atom
    class ContractedGaussian {
    public:
        ContractedGaussian() : mPosition({0.0, 0.0, 0.0}) {
        }

        ContractedGaussian(const ContractedGaussian &col) = default;

        explicit ContractedGaussian(const Eigen::Vector3d &position) : mPosition(position) {
        }

        [[nodiscard]] double Evaluate(const Eigen::Vector3d& pos) const;
        [[nodiscard]] Eigen::Vector3d EvaluateGradient(const Eigen::Vector3d& pos) const;

        Gaussian &Get(uint32_t i) {
            return mGaussians[i];
        }

        const Gaussian &Get(uint32_t i) const {
            return mGaussians[i];
        }

        void Add(const Eigen::Vector3i &angMom, const Gaussian &gaussianIn);

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

        std::vector<Gaussian> &GetGaussians() {
            return mGaussians;
        }

        const std::vector<Gaussian> &GetGaussians() const {
            return mGaussians;
        }

    private:
        Eigen::Vector3d mPosition;
        std::vector<Gaussian> mGaussians;
    };

    template<typename T>
    std::function<T(T, T, T)> GenerateMolecularWavefunction(const Eigen::VectorXd &coeffsIn,
                                                            const BasisSet &basisSet);


}
