#pragma once
#include<Eigen/Core>
#include <unsupported/Eigen/CXX11/Tensor>
#include<unordered_map>

#include "Chemistry/Quantum/Molecule.h"

namespace SCESS {
    class Gaussian;
    class ContractedGaussian;
    class Atom;
    class Molecule;

    //Note that the tei is in chemists’ notation
    using TEITensor = Eigen::Tensor<double, 4>;

    // struct TETensor
    // {
    // public:
    //
    //     TETensor(){}
    //     TETensor(const TETensor& tensorIn):mTEIntegrals(tensorIn.mTEIntegrals) {}
    //
    //     double& operator()(size_t i, size_t j, size_t k, size_t l)
    //     {
    //         size_t index = YoshimineIndex(i, j, k, l);
    //
    //         return mTEIntegrals.operator[](index);
    //     }
    //
    //     double operator()(size_t i, size_t j, size_t k, size_t l) const
    //     {
    //         size_t index = YoshimineIndex(i, j, k, l);
    //
    //         return mTEIntegrals.at(index);
    //     }
    //
    //     double& Get(size_t i, size_t j, size_t k, size_t l)
    //     {
    //         size_t index = YoshimineIndex(i, j, k, l);
    //
    //         return mTEIntegrals.operator[](index);
    //     }
    //
    //     double Get(size_t i, size_t j, size_t k, size_t l) const
    //     {
    //         size_t index = YoshimineIndex(i, j, k, l);
    //
    //         return mTEIntegrals.at(index);
    //     }
    //
    //     void Insert(size_t i, size_t j, size_t k, size_t l, double val)
    //     {
    //         mTEIntegrals[YoshimineIndex(i, j, k, l)] = val;
    //     }
    //
    //     bool Exists(size_t i, size_t j, size_t k, size_t l) const
    //     {
    //         size_t index = YoshimineIndex(i, j, k, l);
    //         return mTEIntegrals.find(index) != mTEIntegrals.end();
    //     }
    //
    //     static size_t YoshimineIndex(size_t i, size_t j, size_t k, size_t l)
    //     {
    //         if (i < j) {
    //             std::swap(i, j);
    //         }
    //         if (k < l) {
    //             std::swap(k, l);
    //         }
    //
    //         size_t ij = i * (i + 1) / 2 + j;
    //         size_t kl = k * (k + 1) / 2 + l;
    //
    //         if (ij < kl) {
    //             std::swap(ij, kl);
    //         }
    //
    //         return ij * (ij + 1) / 2 + kl;
    //     }
    //
    // private:
    //     std::unordered_map<size_t, double> mTEIntegrals;
    // };

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


    class Integrals {
    public:
        class Transforms {
        public:
            /**
             * Transform the two-electron integrals (in chemist's notation) from the atomic orbital to the molecular orbital basis set
             * @param tei Two-electron integral tensor in atomic orbital basis (in chemist's notation)
             * @param C Wavefunction coefficients
             * @return Transformed TEI tensor (in chemist's notation)
             */
            static TEITensor ToMolecularTEI(const TEITensor &tei, const Eigen::MatrixXd &C);
            /**
             * Produce a tensor with < i j || k l>=<i j | k l> - <i j | l k>
             * Automatically converts chemist's to physicists notation
             * @param tei Two-electron integral tensor (in chemist's notation)
             * @return Two-electron integral tensor with antisymmetric integrals (in physicists' notation)
             */
            static TEITensor AntisymmetriseTEI(const TEITensor &tei);
            /**
             * Transforms a molecular two electron integral tensor (in chemist's notation) into a tensor of antisymmetric
             * two-electron integrals with the spin functions integrated out.
             * A closed-shell wavefunction is assumed such that alpha and beta functions alternate.
             * @param moTEI The TEI tensor transformed into the MO basis (in chemist's notation)
             * @return A tensor with elements < i j || k l>= <i j | k l><alpha/beta ...> - <i j | l k><alpha/beta ...>
             */
            static TEITensor AntisymmetriseIntegrateOutSpinClosedShell(const TEITensor &moTEI);
            //
            /**
             *Converts tei from chemists' to physicists' notation i.e. <i j | k l> = ( i k | j l)
             * @param tei The TEI tensor in chemist's notation
             * @return The TEI tensor in physicists' notation
             */
            static TEITensor ConvertChemistsToPhysicists(const TEITensor &tei);
        };

        //Helper functions
        static double Binomial(int64_t a, int64_t b);
        static double PolynomialCoeff(int s, int ia, int ib,
                                      double xpa, double xpb);
        static Eigen::Vector3d GaussianCenter(double exponent1, double exponent2,
                                              const Eigen::Vector3d &pos1, const Eigen::Vector3d &pos2);

        static size_t TEIYoshimineIdx(size_t i, size_t j, size_t k, size_t l);

        /**
         * Evaluate all the required integrals for a qunatum calculation on a given molecule.
         * These include the overlap, kinetic, potential and two-electron integrals.
         * @param molecule The molecule to evaluate the integrals on
         * @return Configuration integrals
         */
        static ConfigIntegrals CalculateConfiguration(const Molecule &molecule);

        //Overlap Integrals
        static double Overlap1D(double gamma,
                                int angX1, double posX1,
                                int angX2, double posX2);
        static double OverlapGaussianParam(double exponent1, double exponent2, const Eigen::Vector3i &angular1,
                                           const Eigen::Vector3i &angular2, const Eigen::Vector3d &pos1,
                                           const Eigen::Vector3d &pos2);
        static double OverlapGaussian(const Gaussian &gto1, const Gaussian &gto2);
        static double EvaluateOverlap(const ContractedGaussian &c1, const ContractedGaussian &c2);
        static Eigen::MatrixXd CalculateOverlapMatrix(const Molecule &molecule);
        static Eigen::MatrixXd CalculateMixedOverlapMatrix(const Molecule &mol1, const Molecule &mol2);
        //Nuclear attraction integrals
        //A-matrix style integration
        //https://csclub.uwaterloo.ca/~pbarfuss/jpsj.21.2313.pdf
        static double AFactor(const int i, const int r, const int u,
                              const int ang1, const int ang2, const double posA_P, double posB_P, double posNuc_P,
                              double gamma);
        static std::vector<double> GenerateACoeffs(int ang1, int ang2,
                                                   double posA_P, double posB_P, double posNuc_P, double gamma);
        static double NuclearGaussianParam(double exp1, double exp2, const Eigen::Vector3d &pos1,
                                           const Eigen::Vector3d &pos2,
                                           const Eigen::Vector3i &ang1, const Eigen::Vector3i &ang2,
                                           const Eigen::Vector3d &nuclearPos);
        static double NuclearGaussian(const Gaussian &gto1, const Gaussian &gto2, const Eigen::Vector3d &nuclearPos);
        static double EvaluateNuclear(const ContractedGaussian &c1, const ContractedGaussian &c2,
                                      const Eigen::Vector3d &nuclearPos,
                                      const int charge);
        static Eigen::MatrixXd CalculateNuclearMatrix(const Molecule &molecule);

        //Kinetic energy integrals
        static double KineticGaussian(const Gaussian &gto1, const Gaussian &gto2);
        static double EvaluateKinetic(const ContractedGaussian &c1, const ContractedGaussian &c2);
        static Eigen::MatrixXd CalculateKineticMatrix(const Molecule &molecule);

        /**
         * Evaluate the two-electron integral (tei) tensor for a molecule in chemist's notation
         * @param molecule The molecule to evaluate the tei for
         * @return The tei in chemist's notation
         */
        static TEITensor CalculateTwoElectronTensor(const Molecule &molecule);
        static TEITensor CalculateTwoElectronTensorMultithread(const Molecule &molecule);
        static double EvaluateTwoElectron(const ContractedGaussian &c1, const ContractedGaussian &c2,
                                          const ContractedGaussian &c3, const ContractedGaussian &c4);
        static double TwoElectronGaussian(const Gaussian &gto1, const Gaussian &gto2,
                                          const Gaussian &gto3, const Gaussian &gto4);
        static double TwoElectronGaussian(const Eigen::Vector3d &a, const int la, const int ma, const int na,
                                          const double alphaa,
                                          const Eigen::Vector3d &b, const int lb, const int mb, const int nb,
                                          const double alphab,
                                          const Eigen::Vector3d &c, const int lc, const int mc, const int nc,
                                          const double alphac,
                                          const Eigen::Vector3d &d, const int ld, const int md, const int nd,
                                          const double alphad);
        //https://csclub.uwaterloo.ca/~pbarfuss/jpsj.21.2313.pdf
        static std::vector<double> B_Array(const int l1, const int l2, const int l3, const int l4,
                                           const double p, const double a, const double b, const double q,
                                           const double c, const double d,
                                           const double g1, const double g2, const double delta);
        // static double B_term(const int i1, const int i2, const int r1, const int r2, const int u, const int l1,
        //                      const int l2, const int l3, const int l4,
        //                      const double px, const double ax, const double bx, const double qx, const double cx,
        //                      const double dx, const double gamma1,
        //                      const double gamma2, const double delta);
        // static double fB(const int i, const int l1, const int l2, const double p,
        //                  const double a, const double b, const int r, const double g);
        // static double B0(int i, int r, double g);
        // static double fact_ratio2(const int a, const int b);
    };
}
