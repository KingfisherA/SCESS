#include<pch.h>
#include "Gaussian.h"

#include "Chemistry/Quantum/Molecule.h"

#include "Platform/Intrinsics.h"
#include<immintrin.h>
#include<type_traits>


namespace SCESS {
    using namespace Eigen;
    using std::numbers::pi;
    using boost::math::double_factorial;

    void Gaussian::CalculateNormalisationConstant() {
        double numerator = std::pow(2 * mExponent / pi, 0.75) *
                           std::pow(4 * mExponent, (mAngular[0] + mAngular[1] + mAngular[2]) * 0.5);

        double denominator = std::sqrt(
            (mAngular[0] < 1 ? 1 : double_factorial<double>(2 * mAngular[0] - 1)) *
            (mAngular[1] < 1 ? 1 : double_factorial<double>(2 * mAngular[1] - 1)) *
            (mAngular[2] < 1 ? 1 : double_factorial<double>(2 * mAngular[2] - 1)));

        mNorm = numerator / denominator;
    }

    double ContractedGaussian::Evaluate(const Vector3d &pos) const {
        double result = 0.0;
        for (const Gaussian &g: mGaussians) {
            Vector3d relativePos = pos - g.mPosition;

            result += g.mCoeff * g.mNorm
                    * std::pow(relativePos.x(), g.mAngular[0])
                    * std::pow(relativePos.y(), g.mAngular[1])
                    * std::pow(relativePos.z(), g.mAngular[2])
                    * std::exp(-g.mExponent * relativePos.squaredNorm());
        }
        return result;
    }

    Vector3d ContractedGaussian::EvaluateGradient(const Vector3d &pos) const {
        Vector3d gradientVector(0, 0, 0);

        for (const Gaussian &g: mGaussians) {
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

    void ContractedGaussian::Add(const Eigen::Vector3i &angMom, const Gaussian &gaussianIn) {
        Gaussian g(gaussianIn.mExponent, gaussianIn.mCoeff, this->mPosition, angMom);
        g.CalculateNormalisationConstant();
        mGaussians.push_back(g);
    }

    template<typename T>
    std::function<T(T, T, T)> GenerateMolecularWavefunction(const Eigen::VectorXd &coeffsIn,
                                                            const Molecule &mol) {
        static_assert(std::is_same_v<T, float> || std::is_same_v<T, double>,
                      "Can only instantiate template with a floating point type");
        return [
                    coeffs = coeffsIn,
                    &molecule = static_cast<const Molecule &>(mol),
                    gaussianData = std::vector<T>()
                ](T x, T y, T z) mutable -> T {
            // #ifdef __AVX2__
            // 				if constexpr (std::is_same_v<T, float>)
            // 				{
            // 					if (gaussianData.size() == 0)
            // 					{
            // 						size_t dataSize = 0;
            // 						for (const GaussianCollection& gc : molecule.GetOrbitals())
            // 						{
            // 							//3 for angular, 3 for position, 1 for exponent, 1 for normalisation
            // 							//1 for coefficient, 1 for MO coefficient
            // 							dataSize += gc.Size() * 10;
            // 						}
            //
            // 						//Ensure data set is divisible by 8
            // 						dataSize = dataSize + (80 - dataSize % 80);
            // 						gaussianData.resize(dataSize, 0.0f);
            //
            // 						size_t idx = 0;
            // 						for (size_t i = 0; i < molecule.GetNoOrbitals(); i++)
            // 						{
            // 							const GaussianCollection& gc = molecule.GetOrbital(i);
            // 							for (size_t j = 0; j < gc.Size(); j++)
            // 							{
            // 								const Gaussian& g = gc.Get(j);
            // 								gaussianData[idx + 0] = (float)coeffs[i];
            // 								gaussianData[idx + 1] = (float)g.mCoeff;
            // 								gaussianData[idx + 2] = (float)g.mNorm;
            // 								gaussianData[idx + 3] = (float)g.mExponent;
            // 								gaussianData[idx + 4] = (float)g.mPosition[0];
            // 								gaussianData[idx + 5] = (float)g.mPosition[1];
            // 								gaussianData[idx + 6] = (float)g.mPosition[2];
            // 								gaussianData[idx + 7] = (float)g.mAngular[0];
            // 								gaussianData[idx + 8] = (float)g.mAngular[1];
            // 								gaussianData[idx + 9] = (float)g.mAngular[2];
            // 								idx += 10;
            // 							}
            // 						}
            // 					}
            // 					__m256 _x0, _y0, _z0, _moCoeff, _norm, _coeff, _ang1, _ang2, _ang3, _exponent1,
            // 						_exponent2, _result;
            // 					_result = _mm256_set1_ps(0.0f);
            // 					//Intrinsics can perform 8 operations at any one time.
            // 					for (size_t i = 0; i < gaussianData.size(); i += 80)
            // 					{
            // 						/*float x0 = x - gaussianData[idx + 4];
            // 						float y0 = y - gaussianData[idx + 5];
            // 						float z0 = z - gaussianData[idx + 6];*/
            //
            // 						//Index 4 to 6 = Position (x,y,z) of gaussian
            // 						_x0 = _mm256_set_ps(
            // 							gaussianData[i + 0 + 4],
            // 							gaussianData[i + 10 + 4],
            // 							gaussianData[i + 20 + 4],
            // 							gaussianData[i + 30 + 4],
            // 							gaussianData[i + 40 + 4],
            // 							gaussianData[i + 50 + 4],
            // 							gaussianData[i + 60 + 4],
            // 							gaussianData[i + 70 + 4]);
            // 						_y0 = _mm256_set_ps(
            // 							gaussianData[i + 0 + 5],
            // 							gaussianData[i + 10 + 5],
            // 							gaussianData[i + 20 + 5],
            // 							gaussianData[i + 30 + 5],
            // 							gaussianData[i + 40 + 5],
            // 							gaussianData[i + 50 + 5],
            // 							gaussianData[i + 60 + 5],
            // 							gaussianData[i + 70 + 5]);
            // 						_z0 = _mm256_set_ps(
            // 							gaussianData[i + 0 + 6],
            // 							gaussianData[i + 10 + 6],
            // 							gaussianData[i + 20 + 6],
            // 							gaussianData[i + 30 + 6],
            // 							gaussianData[i + 40 + 6],
            // 							gaussianData[i + 50 + 6],
            // 							gaussianData[i + 60 + 6],
            // 							gaussianData[i + 70 + 6]);
            //
            // 						_x0 = _mm256_sub_ps(_mm256_set1_ps(x), _x0);
            // 						_y0 = _mm256_sub_ps(_mm256_set1_ps(y), _y0);
            // 						_z0 = _mm256_sub_ps(_mm256_set1_ps(z), _z0);
            //
            // 						//Index 0 = MO coefficient
            // 						_moCoeff = _mm256_set_ps(
            // 							gaussianData[i + 0 + 0],
            // 							gaussianData[i + 10 + 0],
            // 							gaussianData[i + 20 + 0],
            // 							gaussianData[i + 30 + 0],
            // 							gaussianData[i + 40 + 0],
            // 							gaussianData[i + 50 + 0],
            // 							gaussianData[i + 60 + 0],
            // 							gaussianData[i + 70 + 0]);
            // 						//Index 2 = Normalisation Constant
            // 						_norm = _mm256_set_ps(
            // 							gaussianData[i + 0 + 2],
            // 							gaussianData[i + 10 + 2],
            // 							gaussianData[i + 20 + 2],
            // 							gaussianData[i + 30 + 2],
            // 							gaussianData[i + 40 + 2],
            // 							gaussianData[i + 50 + 2],
            // 							gaussianData[i + 60 + 2],
            // 							gaussianData[i + 70 + 2]);
            // 						//Index 1 = Gaussian Coefficient
            // 						_coeff = _mm256_set_ps(
            // 							gaussianData[i + 0 + 1],
            // 							gaussianData[i + 10 + 1],
            // 							gaussianData[i + 20 + 1],
            // 							gaussianData[i + 30 + 1],
            // 							gaussianData[i + 40 + 1],
            // 							gaussianData[i + 50 + 1],
            // 							gaussianData[i + 60 + 1],
            // 							gaussianData[i + 70 + 1]);
            // 						//Index 7 to 9 = Angular coefficients
            // 						_ang1 = _mm256_set_ps(
            // 							gaussianData[i + 0 + 7],
            // 							gaussianData[i + 10 + 7],
            // 							gaussianData[i + 20 + 7],
            // 							gaussianData[i + 30 + 7],
            // 							gaussianData[i + 40 + 7],
            // 							gaussianData[i + 50 + 7],
            // 							gaussianData[i + 60 + 7],
            // 							gaussianData[i + 70 + 7]);
            // 						_ang2 = _mm256_set_ps(
            // 							gaussianData[i + 0 + 8],
            // 							gaussianData[i + 10 + 8],
            // 							gaussianData[i + 20 + 8],
            // 							gaussianData[i + 30 + 8],
            // 							gaussianData[i + 40 + 8],
            // 							gaussianData[i + 50 + 8],
            // 							gaussianData[i + 60 + 8],
            // 							gaussianData[i + 70 + 8]);
            // 						_ang3 = _mm256_set_ps(
            // 							gaussianData[i + 0 + 9],
            // 							gaussianData[i + 10 + 9],
            // 							gaussianData[i + 20 + 9],
            // 							gaussianData[i + 30 + 9],
            // 							gaussianData[i + 40 + 9],
            // 							gaussianData[i + 50 + 9],
            // 							gaussianData[i + 60 + 9],
            // 							gaussianData[i + 70 + 9]);
            // 						//Index 3 = Gaussian Exponent
            // 						_exponent1 = _mm256_set_ps(
            // 							gaussianData[i + 0 + 3],
            // 							gaussianData[i + 10 + 3],
            // 							gaussianData[i + 20 + 3],
            // 							gaussianData[i + 30 + 3],
            // 							gaussianData[i + 40 + 3],
            // 							gaussianData[i + 50 + 3],
            // 							gaussianData[i + 60 + 3],
            // 							gaussianData[i + 70 + 3]);
            // 						_exponent2 = _mm256_set1_ps(0.0f);
            //
            // 						_moCoeff = _mm256_mul_ps(_moCoeff, _norm);
            // 						_moCoeff = _mm256_mul_ps(_coeff, _moCoeff);
            // 						_moCoeff = _mm256_mul_ps(_mm256_pow_ps(_x0, _ang1), _moCoeff);
            // 						_moCoeff = _mm256_mul_ps(_mm256_pow_ps(_y0, _ang2), _moCoeff);
            // 						_moCoeff = _mm256_mul_ps(_mm256_pow_ps(_z0, _ang3), _moCoeff);
            // 						_exponent2 = _mm256_fmadd_ps(_x0, _x0, _exponent2);
            // 						_exponent2 = _mm256_fmadd_ps(_y0, _y0, _exponent2);
            // 						_exponent2 = _mm256_fmadd_ps(_z0, _z0, _exponent2);
            // 						_exponent2 = _mm256_mul_ps(_mm256_sub_ps(_mm256_set1_ps(0.0f), _exponent1), _exponent2);
            // 						_exponent2 = _mm256_exp_ps(_exponent2);
            // 						_result = _mm256_fmadd_ps(_moCoeff, _exponent2, _result);
            //
            // 						/*	result += gaussianData[idx + 0] * gaussianData[idx + 2] * gaussianData[idx + 1] *
            // 								std::pow(x0, gaussianData[idx + 7])
            // 								* std::pow(y0, gaussianData[idx + 8])
            // 								* std::pow(z0, gaussianData[idx + 9])
            // 								* std::exp(-gaussianData[idx + 3] * (x0 * x0 + y0 * y0 + z0 * z0));*/
            // 					}
            // 					return __mm_sum8_ps(_result);
            // 				}
            // 				else if constexpr (std::is_same_v<T, double>)
            // 				{
            // 					if (gaussianData.size() == 0)
            // 					{
            // 						size_t dataSize = 0;
            // 						for (const GaussianCollection& gc : molecule.GetOrbitals())
            // 						{
            // 							//3 for angular, 3 for position, 1 for exponent, 1 for normalisation
            // 							//1 for coefficient, 1 for MO coefficient
            // 							dataSize += gc.Size() * 10;
            // 						}
            //
            // 						//Ensure data set is divisible by 4
            // 						dataSize = dataSize + (40 - dataSize % 40);
            // 						gaussianData.resize(dataSize, 0.0);
            //
            // 						for (size_t i = 0; i < molecule.GetOrbitals().size(); i++)
            // 						{
            // 							const GaussianCollection& gc = molecule.GetOrbital(i);
            // 							for (size_t j = 0; j < gc.Size(); j++)
            // 							{
            // 								const Gaussian& g = gc.Get(j);
            // 								size_t idx = i * gc.Size() * 10 + j * 10;
            // 								gaussianData[idx + 0] = (double)coeffs[i];
            // 								gaussianData[idx + 1] = (double)g.mCoeff;
            // 								gaussianData[idx + 2] = (double)g.mNorm;
            // 								gaussianData[idx + 3] = (double)g.mExponent;
            // 								gaussianData[idx + 4] = (double)g.mPosition[0];
            // 								gaussianData[idx + 5] = (double)g.mPosition[1];
            // 								gaussianData[idx + 6] = (double)g.mPosition[2];
            // 								gaussianData[idx + 7] = (double)g.mAngular[0];
            // 								gaussianData[idx + 8] = (double)g.mAngular[1];
            // 								gaussianData[idx + 9] = (double)g.mAngular[2];
            // 							}
            // 						}
            // 					}
            //
            // 					__m256d _x0, _y0, _z0, _moCoeff, _norm, _coeff, _ang1, _ang2, _ang3, _exponent1,
            // 						_exponent2, _result;
            // 					_result = _mm256_set1_pd(0.0);
            // 					//Intrinsics can perform 8 operations at any one time.
            // 					for (size_t i = 0; i < gaussianData.size(); i += 40)
            // 					{
            // 						/*double x0 = x - gaussianData[idx + 4];
            // 						double y0 = y - gaussianData[idx + 5];
            // 						double z0 = z - gaussianData[idx + 6];*/
            //
            // 						//Index 4 to 6 = Position (x,y,z) of gaussian
            // 						_x0 = _mm256_set_pd(
            // 							gaussianData[i + 0 + 4],
            // 							gaussianData[i + 10 + 4],
            // 							gaussianData[i + 20 + 4],
            // 							gaussianData[i + 30 + 4]);
            // 						_y0 = _mm256_set_pd(
            // 							gaussianData[i + 0 + 5],
            // 							gaussianData[i + 10 + 5],
            // 							gaussianData[i + 20 + 5],
            // 							gaussianData[i + 30 + 5]);
            // 						_z0 = _mm256_set_pd(
            // 							gaussianData[i + 0 + 6],
            // 							gaussianData[i + 10 + 6],
            // 							gaussianData[i + 20 + 6],
            // 							gaussianData[i + 30 + 6]);
            //
            // 						_x0 = _mm256_sub_pd(_mm256_set1_pd(x), _x0);
            // 						_y0 = _mm256_sub_pd(_mm256_set1_pd(y), _y0);
            // 						_z0 = _mm256_sub_pd(_mm256_set1_pd(z), _z0);
            //
            // 						//Index 0 = MO coefficient
            // 						_moCoeff = _mm256_set_pd(
            // 							gaussianData[i + 0 + 0],
            // 							gaussianData[i + 10 + 0],
            // 							gaussianData[i + 20 + 0],
            // 							gaussianData[i + 30 + 0]);
            // 						//Index 2 = Normalisation Constant
            // 						_norm = _mm256_set_pd(
            // 							gaussianData[i + 0 + 2],
            // 							gaussianData[i + 10 + 2],
            // 							gaussianData[i + 20 + 2],
            // 							gaussianData[i + 30 + 2]);
            // 						//Index 1 = Gaussian Coefficient
            // 						_coeff = _mm256_set_pd(
            // 							gaussianData[i + 0 + 1],
            // 							gaussianData[i + 10 + 1],
            // 							gaussianData[i + 20 + 1],
            // 							gaussianData[i + 30 + 1]);
            // 						//Index 7 to 9 = Angular coefficients
            // 						_ang1 = _mm256_set_pd(
            // 							gaussianData[i + 0 + 7],
            // 							gaussianData[i + 10 + 7],
            // 							gaussianData[i + 20 + 7],
            // 							gaussianData[i + 30 + 7]);
            // 						_ang2 = _mm256_set_pd(
            // 							gaussianData[i + 0 + 8],
            // 							gaussianData[i + 10 + 8],
            // 							gaussianData[i + 20 + 8],
            // 							gaussianData[i + 30 + 8]);
            // 						_ang3 = _mm256_set_pd(
            // 							gaussianData[i + 0 + 9],
            // 							gaussianData[i + 10 + 9],
            // 							gaussianData[i + 20 + 9],
            // 							gaussianData[i + 30 + 9]);
            // 						//Index 3 = Gaussian Exponent
            // 						_exponent1 = _mm256_set_pd(
            // 							gaussianData[i + 0 + 3],
            // 							gaussianData[i + 10 + 3],
            // 							gaussianData[i + 20 + 3],
            // 							gaussianData[i + 30 + 3]);
            // 						_exponent2 = _mm256_set1_pd(0.0);
            //
            // 						_moCoeff = _mm256_mul_pd(_moCoeff, _norm);
            // 						_moCoeff = _mm256_mul_pd(_coeff, _moCoeff);
            // 						_moCoeff = _mm256_mul_pd(_mm256_pow_pd(_x0, _ang1), _moCoeff);
            // 						_moCoeff = _mm256_mul_pd(_mm256_pow_pd(_y0, _ang2), _moCoeff);
            // 						_moCoeff = _mm256_mul_pd(_mm256_pow_pd(_z0, _ang3), _moCoeff);
            // 						_exponent2 = _mm256_fmadd_pd(_x0, _x0, _exponent2);
            // 						_exponent2 = _mm256_fmadd_pd(_y0, _y0, _exponent2);
            // 						_exponent2 = _mm256_fmadd_pd(_z0, _z0, _exponent2);
            // 						_exponent2 = _mm256_mul_pd(_mm256_sub_pd(_mm256_set1_pd(0.0), _exponent1), _exponent2);
            // 						_exponent2 = _mm256_exp_pd(_exponent2);
            // 						_result = _mm256_fmadd_pd(_moCoeff, _exponent2, _result);
            //
            // 						/*	result += gaussianData[idx + 0] * gaussianData[idx + 2] * gaussianData[idx + 1] *
            // 								std::pow(x0, gaussianData[idx + 7])
            // 								* std::pow(y0, gaussianData[idx + 8])
            // 								* std::pow(z0, gaussianData[idx + 9])
            // 								* std::exp(-gaussianData[idx + 3] * (x0 * x0 + y0 * y0 + z0 * z0));*/
            // 					}
            // 					return __mm_sum4_pd(_result);
            // 				}
            // #else
            if (gaussianData.size() == 0) {
                size_t dataSize = 0;
                for (const ContractedGaussian &gc: molecule.GetOrbitals()) {
                    //3 for angular, 3 for position, 1 for exponent, 1 for normalisation
                    //1 for coefficient, 1 for MO coefficient
                    dataSize += gc.Size() * 10;
                }
                gaussianData.resize(dataSize);

                size_t idx = 0;
                for (size_t i = 0; i < molecule.GetOrbitals().size(); i++) {
                    const ContractedGaussian &gc = molecule.GetOrbital(i);
                    for (size_t j = 0; j < gc.Size(); j++) {
                        const Gaussian &g = gc.Get(j);
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
            //#endif
        };
        /*	return[coeffs = coeffsM->col(mSelectedMO),
                &molecule = static_cast<const Molecule&>(*mMolecule)](float x, float y, float z) -> float
                {
                    float result = 0.0f;
                    for (size_t i = 0; i < molecule.GetOrbitals().size(); i++)
                    {
                        result += (float)coeffs[i] * EvaluateGaussianFunction(molecule.GetOrbital(i), x, y, z);
                    }
                    return result;
                };*/
    }

    template std::function<float(float, float, float)>
    GenerateMolecularWavefunction(const Eigen::VectorXd &coeffs, const Molecule &mol);
    template std::function<double(double, double, double)>
    GenerateMolecularWavefunction(const Eigen::VectorXd &coeffs, const Molecule &mol);
}
