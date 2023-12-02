#pragma once
#include<Eigen/Core>
#include<numbers>
#include<unordered_map>
#include<boost/math/special_functions/factorials.hpp>
#include<vector>

namespace SCESS
{
    /*
    Represents a single gaussian function with three polynomial components

    ..math::
    (x - A_x)^ {
        a_x
    } *(y - A_y)^ { a_y } *(z - A_z)^ { a_z } *\exp{ (-b * ((x + y + z) - (A_x + A_y + A_z)) ^ 2) }

    ...

    Attributes
    ----------
    position : [float, 3]
    The x, y, z position of the gaussian function.Represented by Ax, Ay and Az
    exponent : float
    The exponent coefficient; represented by b
    angular : [float, 3]
    The angular exponents.Represented by ax, ay, az

    Methods
    ------ -
    */
    class Gaussian {
    public:
        Eigen::Vector3d mPosition;
        double mExponent;
        Eigen::Vector3i mAngular;
        double mCoeff;
        double mNorm;

        Gaussian(double exponent, double coeff, const Eigen::Vector3d& position = { 0., 0., 0. },
            const Eigen::Vector3i& angular = { 0, 0, 0 }) :mPosition(position), mExponent(exponent), mAngular(angular), mCoeff(coeff), mNorm(1.0)
        {
        }

        void CalculateNormalisationConstant()
        {
       /*     double numerator = std::pow(2.0, 2.0 * (mAngular[0] + mAngular[1] + mAngular[2]) + 3.0 / 2.0) *
                std::pow(mExponent, (mAngular[0] + mAngular[1] + mAngular[2]) + 3.0 / 2.0);*/

            double numerator = std::pow(2 * mExponent / std::numbers::pi, 0.75) * std::pow(4 * mExponent, (mAngular[0] + mAngular[1] + mAngular[2]) * 0.5);


            double denominator = std::sqrt((mAngular[0] < 1 ? 1 : boost::math::double_factorial<double>(2 * mAngular[0] - 1)) *
                (mAngular[1] < 1 ? 1 : boost::math::double_factorial<double>(2 * mAngular[1] - 1)) *
                (mAngular[2] < 1 ? 1 : boost::math::double_factorial<double>(2 * mAngular[2] - 1)));

            mNorm = numerator / denominator;
        }

        const Eigen::Vector3i& GetAngularMomentum() const {
            return mAngular;
        }
    };

    enum class OrbitalType {
        S,
        //
        P_X,
        P_Y,
        P_Z,
        //
        D_Z2,
        D_X2,
        D_XY,
        D_XZ,
        D_YZ,
        D_Y2,
        //
        F_Z3,
        F_Y3,
        F_X3,
        F_X2_Y,
        F_X2_Z,
        F_X_Y2,
        F_XYZ,
        F_X_Z2,
        F_Y2_Z,
        F_Y_Z2
    };


    //Represents a collection of Gaussian functions with the same position which form the basis for a single atom
    class GaussianCollection {
    private:

        Eigen::Vector3d mPosition;
        std::vector<Gaussian> mGaussians;
        static const std::unordered_map<OrbitalType,
            Eigen::Vector3i> sAngMomMap;
    public:
       

        GaussianCollection() : mPosition({ 0.0, 0.0, 0.0 }) {}
        GaussianCollection(const GaussianCollection& col) :mPosition(col.mPosition), mGaussians(col.mGaussians)
        {
        }
        GaussianCollection(const Eigen::Vector3d& position) :mPosition(position) {
        }

        static const Eigen::Vector3i& GetAngMomentum(SCESS::OrbitalType type)
        {
            return sAngMomMap.at(type);
        }

        Gaussian& Get(uint32_t i)
        {
            return mGaussians[i];
        }

        const Gaussian& Get(uint32_t i) const
        {
            return mGaussians[i];
        }

        void Add(const OrbitalType type, const Gaussian& gaussianIn) {
            Eigen::Vector3i angMom = sAngMomMap.at(type);

            Gaussian g(gaussianIn.mExponent, gaussianIn.mCoeff, this->mPosition, angMom);
            g.CalculateNormalisationConstant();
            mGaussians.push_back(g);
        }

        void UpdatePosition(const Eigen::Vector3d& pos)
        {
            mPosition = pos;
            for (auto& g : mGaussians)
            {
                g.mPosition = pos;
            }
        }

        const Eigen::Vector3d& GetPos() const {
            return mPosition;
        }

        size_t Size() const {
            return mGaussians.size();
        }

        std::vector<Gaussian>& GetGaussians() {
            return mGaussians;
        }

        const std::vector<Gaussian>& GetGaussians() const {
            return mGaussians;
        }
    };

}