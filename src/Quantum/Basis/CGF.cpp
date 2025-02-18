#include<Quantum/Q_PCH.h>
#include<type_traits>
//
#include "CGF.h"
#include "Quantum/Basis/BasisSet.h"

//
//
namespace SCESS {
    //
    //
    using namespace Eigen;
    using std::numbers::pi;
    using boost::math::double_factorial;
    //
    void PGF::CalculateNormalisationConstant() {
        double numerator = std::pow(2 * mExponent / pi, 0.75) *
                           std::pow(4 * mExponent, (mAngular[0] + mAngular[1] + mAngular[2]) * 0.5);

        double denominator = std::sqrt(
            (mAngular[0] < 1 ? 1 : double_factorial<double>(2 * mAngular[0] - 1)) *
            (mAngular[1] < 1 ? 1 : double_factorial<double>(2 * mAngular[1] - 1)) *
            (mAngular[2] < 1 ? 1 : double_factorial<double>(2 * mAngular[2] - 1)));

        mNorm = numerator / denominator;
    }
    //
    double CGF::Evaluate(const Vector3d &pos) const {
        double result = 0.0;
        for (const PGF &g: mGaussians) {
            Vector3d relativePos = pos - g.mPosition;

            result += g.mCoeff * g.mNorm
                    * std::pow(relativePos.x(), g.mAngular[0])
                    * std::pow(relativePos.y(), g.mAngular[1])
                    * std::pow(relativePos.z(), g.mAngular[2])
                    * std::exp(-g.mExponent * relativePos.squaredNorm());
        }
        return result;
    }
    //
    Vector3d CGF::EvaluateGradient(const Vector3d &pos) const {
        Vector3d gradientVector(0, 0, 0);

        for (const PGF &g: mGaussians) {
            const double ex = std::exp(-g.mExponent * std::pow(pos[0] - g.mPosition[0], 2));
            const double fx = std::pow(pos[0] - g.mPosition[0], g.mAngular[0]) * ex;

            const double ey = std::exp(-g.mExponent * std::pow(pos[1] - g.mPosition[1], 2));
            const double fy = std::pow(pos[1] - g.mPosition[1], g.mAngular[1]) * ey;

            const double ez = std::exp(-g.mExponent * std::pow(pos[2] - g.mPosition[2], 2));
            const double fz = std::pow(pos[2] - g.mPosition[2], g.mAngular[2]) * ez;

            double gx = -2.0 * g.mExponent * (pos[0] - g.mPosition[0]) * fx;
            double gy = -2.0 * g.mExponent * (pos[1] - g.mPosition[1]) * fy;
            double gz = -2.0 * g.mExponent * (pos[2] - g.mPosition[2]) * fz;

            if (g.mAngular[0] > 0) {
                gx += std::pow(pos[0] - g.mPosition[0], g.mAngular[0] - 1) * ex;
            }
            if (g.mAngular[1] > 0) {
                gy += std::pow(pos[1] - g.mPosition[1], g.mAngular[1] - 1) * ey;
            }
            if (g.mAngular[2] > 0) {
                gz += std::pow(pos[2] - g.mPosition[2], g.mAngular[2] - 1) * ez;
            }

            gradientVector += g.mCoeff * Vector3d(gx * fy * fz,
                                                  fx * gy * fz,
                                                  fx * fy * gz);
        }

        return gradientVector;
    }
    //
    void CGF::Add(const Eigen::Vector3i &angMom, const PGF &gaussianIn) {
        PGF g(gaussianIn.mExponent, gaussianIn.mCoeff, this->mPosition, angMom);
        g.CalculateNormalisationConstant();
        mGaussians.push_back(g);
    }
    //
    template<typename T>
    std::function<T(T, T, T)> GenerateMolecularWavefunction(const Eigen::VectorXd &coeffsIn,
                                                            const BasisSet &basisSet) {
        static_assert(std::is_same_v<T, float> || std::is_same_v<T, double>,
                      "Can only instantiate template with a floating point type");
        return [
                    coeffs = coeffsIn,
                    &basisSet,
                    gaussianData = std::vector<T>()
                ](T x, T y, T z) mutable -> T {
            if (gaussianData.size() == 0) {
                size_t dataSize = 0;
                for (const CGF &gc: basisSet.GetCGFs()) {
                    //3 for angular, 3 for position, 1 for exponent, 1 for normalisation
                    //1 for coefficient, 1 for MO coefficient
                    dataSize += gc.Size() * 10;
                }
                gaussianData.resize(dataSize);

                size_t idx = 0;
                for (size_t i = 0; i < basisSet.GetNoCGFs(); i++) {
                    const CGF &gc = basisSet.Get(i);
                    for (size_t j = 0; j < gc.Size(); j++) {
                        const PGF &g = gc.Get(j);
                        gaussianData[idx + 0] = (T) coeffs[i];
                        gaussianData[idx + 1] = (T) g.mCoeff;
                        gaussianData[idx + 2] = (T) g.mNorm;
                        gaussianData[idx + 3] = (T) g.mExponent;
                        gaussianData[idx + 4] = (T) g.mPosition[0];
                        gaussianData[idx + 5] = (T) g.mPosition[1];
                        gaussianData[idx + 6] = (T) g.mPosition[2];
                        gaussianData[idx + 7] = (T) g.mAngular[0];
                        gaussianData[idx + 8] = (T) g.mAngular[1];
                        gaussianData[idx + 9] = (T) g.mAngular[2];
                        idx += 10;
                    }
                }
            }

            T result = 0.0;
            for (size_t i = 0; i < gaussianData.size(); i += 10) {
                //Index 4 to 6 = Position (x,y,z) of gaussian
                T x0 = x - gaussianData[i + 4];
                T y0 = y - gaussianData[i + 5];
                T z0 = z - gaussianData[i + 6];

                //Index 0 = MO coefficient
                //Index 2 = Normalisation Constant
                //Index 1 = Gaussian Coefficient
                //Index 7 to 9 = Angular coefficients
                //Index 3 = Gaussian Exponent
                result += gaussianData[i + 0] * gaussianData[i + 2] * gaussianData[i + 1] *
                        std::pow(x0, gaussianData[i + 7])
                        * std::pow(y0, gaussianData[i + 8])
                        * std::pow(z0, gaussianData[i + 9])
                        * std::exp(-gaussianData[i + 3] * (x0 * x0 + y0 * y0 + z0 * z0));
            }
            return result;
        };
    }
    //
    template std::function<float(float, float, float)>
    GenerateMolecularWavefunction(const Eigen::VectorXd &coeffs, const BasisSet &basisSet);
    template std::function<double(double, double, double)>
    GenerateMolecularWavefunction(const Eigen::VectorXd &coeffs, const BasisSet &basisSet);
}
