#pragma once
#include<Eigen/Core>
#include<unsupported/Eigen/CXX11/Tensor>

#include<unordered_map>
#include<concepts>

#include "Quantum/Molecule.h"

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
        ): S(S_), T(T_), V(V_), tei(tei_) {
        }

        ConfigIntegrals(const ConfigIntegrals &intIn): S(intIn.S),
                                                       T(intIn.T),
                                                       V(intIn.V),
                                                       tei(intIn.tei) {
        }

        ConfigIntegrals(ConfigIntegrals &&intIn) noexcept: S(std::move(intIn.S)),
                                                           T(std::move(intIn.T)),
                                                           V(std::move(intIn.V)),
                                                           tei(std::move(intIn.tei)) {
        }

        // S = Overlap, T = Kinetic, V = Nuclear
        Eigen::MatrixXd S, T, V;
        //Two-electron tensor
        TEITensor tei;
    };


    namespace Integrals {
        namespace Transforms {
            /**
             * Transform the two-electron integrals (in chemist's notation) from the atomic orbital to the molecular orbital basis set
             * @param tei Two-electron integral tensor in atomic orbital basis (in chemist's notation)
             * @param C Wavefunction coefficients
             * @return Transformed TEI tensor (in chemist's notation)
             */
            TEITensor ToMolecularTEI(const TEITensor &tei, const Eigen::MatrixXd &C);
            /**
             * Produce a tensor with < i j || k l>=<i j | k l> - <i j | l k>
             * Automatically converts chemist's to physicists notation
             * @param tei Two-electron integral tensor (in chemist's notation)
             * @return Two-electron integral tensor with antisymmetric integrals (in physicists' notation)
             */
            TEITensor AntisymmetriseTEI(const TEITensor &tei);
            /**
             * Transforms a molecular two electron integral tensor (in chemist's notation) into a tensor of antisymmetric
             * two-electron integrals with the spin functions integrated out.
             * A closed-shell wavefunction is assumed such that alpha and beta functions alternate.
             * @param moTEI The TEI tensor transformed into the MO basis (in chemist's notation)
             * @return A tensor with elements < i j || k l>= <i j | k l><alpha/beta ...> - <i j | l k><alpha/beta ...>
             */
            TEITensor AntisymmetriseIntegrateOutSpinClosedShell(const TEITensor &moTEI);
            //
            /**
             *Converts tei from chemists' to physicists' notation i.e. <i j | k l> = ( i k | j l)
             * @param tei The TEI tensor in chemist's notation
             * @return The TEI tensor in physicists' notation
             */
            TEITensor ConvertChemistsToPhysicists(const TEITensor &tei);
        };

        void InitLibInt();
        void CloseLibInt();

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

        /**
         * Evaluate all the required integrals for a qunatum calculation on a given molecule.
         * These include the overlap, kinetic, potential and two-electron integrals.
         * @param molecule The molecule to evaluate the integrals on
         * @return Configuration integrals
         */
        ConfigIntegrals CalculateConfiguration(const Molecule &molecule);

        //Overlap Integrals
        double Overlap1D(double gamma,
                         int angX1, double posX1,
                         int angX2, double posX2);
        double OverlapGaussianParam(double exponent1, double exponent2, const Eigen::Vector3i &angular1,
                                    const Eigen::Vector3i &angular2, const Eigen::Vector3d &pos1,
                                    const Eigen::Vector3d &pos2);
        double OverlapGaussian(const Gaussian &gto1, const Gaussian &gto2);
        double EvaluateOverlap(const ContractedGaussian &c1, const ContractedGaussian &c2);
        Eigen::MatrixXd CalculateOverlapMatrix(const Molecule &molecule);
        Eigen::MatrixXd CalculateMixedOverlapMatrix(const Molecule &mol1, const Molecule &mol2);

        using AArray = std::vector<double>;
        //Nuclear attraction integrals
        //A-matrix style integration
        //https://csclub.uwaterloo.ca/~pbarfuss/jpsj.21.2313.pdf
        void A_Array(AArray &aArray, int ang1, int ang2,
                     double posA_P, double posB_P, double posNuc_P, double gamma);
        // double NuclearGaussianParam(double exp1, double exp2, const Eigen::Vector3d &pos1,
        //                             const Eigen::Vector3d &pos2,
        //                             const Eigen::Vector3i &ang1, const Eigen::Vector3i &ang2,
        //                             const Eigen::Vector3d &nuclearPos);
        double NuclearGaussian(const Gaussian &gto1, const Gaussian &gto2, const Eigen::Vector3d &nuclearPos,
                               AArray &ax, AArray &ay, AArray &az);
        double EvaluateNuclear(const ContractedGaussian &c1, const ContractedGaussian &c2,
                               const Eigen::Vector3d &nuclearPos,
                               int charge, AArray &ax, AArray &ay, AArray &az);
        Eigen::MatrixXd CalculateNuclearMatrix(const Molecule &molecule);

        //Kinetic energy integrals
        double KineticGaussian(const Gaussian &gto1, const Gaussian &gto2);
        double EvaluateKinetic(const ContractedGaussian &c1, const ContractedGaussian &c2);
        Eigen::MatrixXd CalculateKineticMatrix(const Molecule &molecule);

        using BArray = std::vector<double>;
        /**
         * Evaluate the two-electron repulsion integral (tei) tensor for a molecule in chemist's notation
         * @param molecule The molecule to evaluate the tei for
         * @return The tei in chemist's notation
         */
        TEITensor TwoElectronRepulsionTensor(const Molecule &molecule);
        TEITensor TwoElectronRepulsionTensorLibInt(const Molecule &mol, std::string_view basisSet);
        TEITensor TwoElectronRepulsionTensorMultithread(const Molecule &molecule);
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
