#pragma once
#include<Eigen/Core>
#include<unsupported/Eigen/CXX11/Tensor>
#include<concepts>

#include "Quantum/Molecule.h"

namespace GEng {
    class Thread;
}

namespace SCESS {
    class Gaussian;
    class ContractedGaussian;
    class Atom;
    class Molecule;

    //Note that the tei is in chemists’ notation
    using TEITensor = Eigen::Tensor<double, 4>;

    struct ConfigIntegrals {
        ConfigIntegrals() = default;

        ConfigIntegrals(const Eigen::MatrixXd &S_,
                        const Eigen::MatrixXd &T_,
                        const Eigen::MatrixXd &V_,
                        const TEITensor &tei_
        ): S(S_), T(T_), V(V_), eri(tei_) {
        }

        ConfigIntegrals(const ConfigIntegrals &intIn): S(intIn.S),
                                                       T(intIn.T),
                                                       V(intIn.V),
                                                       eri(intIn.eri) {
        }

        ConfigIntegrals(ConfigIntegrals &&intIn) noexcept: S(std::move(intIn.S)),
                                                           T(std::move(intIn.T)),
                                                           V(std::move(intIn.V)),
                                                           eri(std::move(intIn.eri)) {
        }

        // S = Overlap, T = Kinetic, V = Nuclear
        Eigen::MatrixXd S, T, V;
        //Two-electron tensor
        TEITensor eri;
    };

    namespace Integrals {
        namespace Utils {
            //Helper functions
            /**
             * Raises -1 to the given integer power
             * @tparam intT
             * @param power
             * @return
             */
            template<std::integral intT>
            inline intT powMinus1(intT power) {
                if (power % 2 == 0)
                    return 1;
                return -1;
            }

            inline Eigen::Vector3d GaussianCenter(double exponent1, double exponent2,
                                                  const Eigen::Vector3d &pos1, const Eigen::Vector3d &pos2) {
                return (exponent1 * pos1 + exponent2 * pos2) / (exponent1 + exponent2);
            }

            template<std::integral intT>
            inline intT TEIYoshimineIdx(intT i, intT j, intT k, intT l) {
                if (i < j) {
                    std::swap(i, j);
                }
                if (k < l) {
                    std::swap(k, l);
                }

                intT ij = i * (i + 1) / 2 + j;
                intT kl = k * (k + 1) / 2 + l;

                if (ij < kl) {
                    std::swap(ij, kl);
                }

                return ij * (ij + 1) / 2 + kl;
            }

            //Evaluate the Boys function
            inline double BoysFunc(double m, double x) {
                //0.00000001
                x = std::max(std::abs(x), 1E-8);
                double val = boost::math::tgamma_lower(m + 0.5, x);
                return 0.5 * std::pow(x, -m - 0.5) * val;
            }

            template<std::integral intT>
            double Factorial(intT n) {
                double result = 1.0;

                for (size_t i = 0; i < n; i++)
                    result *= static_cast<double>(i + 1);

                return result;
            }

            template<std::integral intT>
            double BinomialCoefficient(intT n, intT k) {
                return Factorial(n) / (Factorial(n - k) * Factorial(k));
            }

            inline double PolynomialCoeff(int s, int ia, int ib, double xpa, double xpb) {
                double result = 0.0;

                for (int t = 0; t < s + 1; t++)
                    if ((s - ia <= t) && (t <= ib)) {
                        result += BinomialCoefficient(ia, s - t) *
                                BinomialCoefficient(ib, t) *
                                pow(xpa, ia - s + t) *
                                pow(xpb, ib - t);
                    }

                return result;
            }
        }

        /**
         * Evaluate all the required integrals for a qunatum calculation on a given molecule.
         * These include the overlap, kinetic, potential and two-electron integrals.
         * @param molecule The molecule to evaluate the integrals on
         * @return Configuration integrals
         */
        ConfigIntegrals CalculateConfiguration(const BasisSet &basisSet,
                                               const std::vector<Atom> &atoms,
                                               const GEng::Thread *thread = nullptr);

        //Overlap Integrals
        double Overlap1D(double gamma,
                         int angX1, double posX1,
                         int angX2, double posX2);
        double OverlapGaussian(const Gaussian &gto1, const Gaussian &gto2);
        double EvaluateOverlap(const ContractedGaussian &c1, const ContractedGaussian &c2);
        Eigen::MatrixXd CalculateOverlapMatrix(const BasisSet &basisSet);
        Eigen::MatrixXd CalculateMixedOverlapMatrix(const BasisSet &set1, const BasisSet &set2);

        using AArray = std::vector<double>;
        //Nuclear attraction integrals
        //https://csclub.uwaterloo.ca/~pbarfuss/jpsj.21.2313.pdf
        void A_Array(AArray &aArray, int ang1, int ang2,
                     double posA_P, double posB_P, double posNuc_P, double gamma);
        double NuclearGaussian(const Gaussian &gto1, const Gaussian &gto2, const Eigen::Vector3d &nuclearPos,
                               AArray &ax, AArray &ay, AArray &az);
        double EvaluateNuclear(const ContractedGaussian &c1, const ContractedGaussian &c2,
                               const Eigen::Vector3d &nuclearPos,
                               int charge, AArray &ax, AArray &ay, AArray &az);
        Eigen::MatrixXd CalculateNuclearMatrix(const BasisSet &basisSet, const std::vector<Atom> &atoms);

        //Kinetic energy integrals
        double KineticGaussian(const Gaussian &gto1, const Gaussian &gto2);
        double EvaluateKinetic(const ContractedGaussian &c1, const ContractedGaussian &c2);
        Eigen::MatrixXd CalculateKineticMatrix(const BasisSet &basisSet);

        using BArray = std::vector<double>;
        /**
         * Evaluate the two-electron repulsion integral (tei) tensor for a molecule in chemist's notation
         * @param molecule The molecule to evaluate the tei for
         * @return The tei in chemist's notation
         */
        TEITensor TwoElectronRepulsionTensor(const BasisSet &basisSet, const GEng::Thread *thread = nullptr);
        TEITensor TwoElectronRepulsionTensorMultithread(const BasisSet &basisSet,
                                                        const GEng::Thread *parentThread = nullptr);
        //Calculates <ij | 1/r12 | kl> (I1) for a contgacted
        double TwoElectronRepulsion(const ContractedGaussian &c1, const ContractedGaussian &c2,
                                    const ContractedGaussian &c3, const ContractedGaussian &c4,
                                    BArray &bx, BArray &by, BArray &bz);
        //Calculates <ij | 1/r12 | kl> (I1)
        double TwoElectronRepulsion(const Gaussian &gto1, const Gaussian &gto2,
                                    const Gaussian &gto3, const Gaussian &gto4,
                                    BArray &bx, BArray &by, BArray &bz);
        //https://csclub.uwaterloo.ca/~pbarfuss/jpsj.21.2313.pdf
        void B_Array(BArray &bArray, int l1, int l2, int l3, int l4,
                     double p, double a, double b, double q,
                     double c, double d,
                     double g1, double g2, double delta);
    };
}
