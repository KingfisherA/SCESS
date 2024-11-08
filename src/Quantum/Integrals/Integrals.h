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
        double PolynomialCoeff(int s, int ia, int ib, double xpa, double xpb);
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
        //Nuclear attraction integrals
        //A-matrix style integration
        //https://csclub.uwaterloo.ca/~pbarfuss/jpsj.21.2313.pdf
        double AFactor(const int i, const int r, const int u,
                       const int ang1, const int ang2, const double posA_P, double posB_P, double posNuc_P,
                       double gamma);
        std::vector<double> GenerateACoeffs(int ang1, int ang2,
                                            double posA_P, double posB_P, double posNuc_P, double gamma);
        double NuclearGaussianParam(double exp1, double exp2, const Eigen::Vector3d &pos1,
                                    const Eigen::Vector3d &pos2,
                                    const Eigen::Vector3i &ang1, const Eigen::Vector3i &ang2,
                                    const Eigen::Vector3d &nuclearPos);
        double NuclearGaussian(const Gaussian &gto1, const Gaussian &gto2, const Eigen::Vector3d &nuclearPos);
        double EvaluateNuclear(const ContractedGaussian &c1, const ContractedGaussian &c2,
                               const Eigen::Vector3d &nuclearPos,
                               const int charge);
        Eigen::MatrixXd CalculateNuclearMatrix(const Molecule &molecule);

        //Kinetic energy integrals
        double KineticGaussian(const Gaussian &gto1, const Gaussian &gto2);
        double EvaluateKinetic(const ContractedGaussian &c1, const ContractedGaussian &c2);
        Eigen::MatrixXd CalculateKineticMatrix(const Molecule &molecule);

        /**
         * Evaluate the two-electron repulsion integral (tei) tensor for a molecule in chemist's notation
         * @param molecule The molecule to evaluate the tei for
         * @return The tei in chemist's notation
         */
        TEITensor TwoElectronRepulsionTensor(const Molecule &molecule);
        TEITensor TwoElectronRepulsionTensorLibInt(const Molecule &mol);
        TEITensor TwoElectronRepulsionTensorMultithread(const Molecule &molecule);
        //Calculates <ij | 1/r12 | kl> (I1) for a contgacted
        double TwoElectronRepulsion(const ContractedGaussian &c1, const ContractedGaussian &c2,
                                    const ContractedGaussian &c3, const ContractedGaussian &c4);
        //Calculates <ij | 1/r12 | kl> (I1)
        double TwoElectronRepulsion(const Gaussian &gto1, const Gaussian &gto2,
                                    const Gaussian &gto3, const Gaussian &gto4);
        //https://csclub.uwaterloo.ca/~pbarfuss/jpsj.21.2313.pdf
        std::vector<double> B_Array(const int l1, const int l2, const int l3, const int l4,
                                    const double p, const double a, const double b, const double q,
                                    const double c, const double d,
                                    const double g1, const double g2, const double delta);
        //   double B_term(const int i1, const int i2, const int r1, const int r2, const int u, const int l1,
        //                      const int l2, const int l3, const int l4,
        //                      const double px, const double ax, const double bx, const double qx, const double cx,
        //                      const double dx, const double gamma1,
        //                      const double gamma2, const double delta);
        //   double fB(const int i, const int l1, const int l2, const double p,
        //                  const double a, const double b, const int r, const double g);
        //   double B0(int i, int r, double g);
        //   double fact_ratio2(const int a, const int b);
    };
}
