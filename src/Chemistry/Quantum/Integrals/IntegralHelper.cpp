//Class for performing integral calculations
#include<pch.h>
#include"IntegralHelper.h"

#include<boost/math/special_functions/binomial.hpp>
#include<boost/math/special_functions/factorials.hpp>

#include"Chemistry/Quantum/Integrals/Gamma.h"
#include"Chemistry/Quantum/Basis/Gaussian.h"
#include"Chemistry/Quantum/Molecule.h"

namespace SCESS {
#pragma region Helper Functions


    TEITensor IntegralHelper::CalculateMolecularTEI(const TEITensor &tei, const Eigen::MatrixXd &C) {
        //Note that the tei tensors are in chemists’ notation

        TEITensor transformed1(tei.dimensions());
        transformed1.setZero();
        TEITensor transformed2(tei.dimensions());
        transformed2.setZero();
        TEITensor transformed3(tei.dimensions());
        transformed3.setZero();
        TEITensor transformed4(tei.dimensions());
        transformed4.setZero();

        Eigen::Index dim = tei.dimension(0);

        for (Eigen::Index s = 0; s < dim; ++s) {
            for (Eigen::Index mu = 0; mu < dim; ++mu)
                for (Eigen::Index nu = 0; nu < dim; ++nu)
                    for (Eigen::Index kap = 0; kap < dim; ++kap)
                        for (Eigen::Index lam = 0; lam < dim; ++lam)
                            //From <mu nu | kap lam> to <mu nu | kap s>
                            //Einsum("mnkl,ls->mnks")
                            transformed1(mu, nu, kap, s) += C(lam, s) * tei(mu, nu, kap, lam);

            for (Eigen::Index r = 0; r < dim; ++r) {
                for (Eigen::Index mu = 0; mu < dim; ++mu)
                    for (Eigen::Index nu = 0; nu < dim; ++nu)
                        for (Eigen::Index kap = 0; kap < dim; ++kap)
                            //From <mu nu | kap s> to <mu nu | r s>
                            //Einsum("mnks,kr->mnrs")
                            transformed2(mu, nu, r, s) += C(kap, r) * transformed1(mu, nu, kap, s);

                for (Eigen::Index q = 0; q < dim; ++q) {
                    for (Eigen::Index mu = 0; mu < dim; ++mu)
                        for (Eigen::Index nu = 0; nu < dim; ++nu)
                            //From <mu nu | r s> to <mu q | r s>
                            //Einsum("mnrs,nq->mqrs")
                            transformed3(mu, q, r, s) += C(nu, q) * transformed2(mu, nu, r, s);

                    for (Eigen::Index p = 0; p < dim; ++p) {
                        for (Eigen::Index mu = 0; mu < dim; ++mu)
                            //From <mu q | r s> to <p q | r s>
                            //Einsum("mqrs,mp->pqrs")
                            transformed4(p, q, r, s) += C(mu, p) * transformed3(mu, q, r, s);
                    }
                }
            }
        }

        return transformed4;
    }

    TEITensor IntegralHelper::AntisymmetriseTEI(const TEITensor &tei) {
        TEITensor transformed(tei.dimensions());

        for (Eigen::Index p = 0; p < transformed.dimension(0); ++p)
            for (Eigen::Index q = 0; q < transformed.dimension(1); ++q)
                for (Eigen::Index r = 0; r < transformed.dimension(2); ++r)
                    for (Eigen::Index s = 0; s < transformed.dimension(3); ++s)
                        transformed(p, q, r, s) = tei(p, r, q, s) - tei(p, s, q, r);

        return transformed;
    }

    TEITensor IntegralHelper::AntisymmetriseIntegrateOutSpinClosedShell(const TEITensor &moTEI) {
        Eigen::DSizes<Eigen::Index, 4> dimensions = moTEI.dimensions();

        //Make the tensor twice as big to hold the spin-integrated integrals
        //Each integral in the moTEI actually represents two electrons
        for (Eigen::Index &dim: dimensions)
            dim *= 2;

        TEITensor transformed(dimensions);

        for (Eigen::Index p = 0; p < transformed.dimension(0); ++p)
            for (Eigen::Index q = 0; q < transformed.dimension(1); ++q)
                for (Eigen::Index r = 0; r < transformed.dimension(2); ++r)
                    for (Eigen::Index s = 0; s < transformed.dimension(3); ++s) {
                        //Only allow the integral to exist if the alpha is integrated with an alpha
                        //or beta with a beta (spin function orthogonality)
                        Eigen::Index coloumbMask = (p % 2 == r % 2) * (q % 2 == s % 2);
                        Eigen::Index exchangeMask = (p % 2 == s % 2) * (q % 2 == r % 2);

                        transformed(p, q, r, s) =
                                moTEI(p / 2, r / 2, q / 2, s / 2) * (double) coloumbMask
                                - moTEI(p / 2, s / 2, q / 2, r / 2) * (double) exchangeMask;
                    }

        return transformed;
    }

    TEITensor IntegralHelper::ConvertChemistsToPhysicists(const TEITensor &tei) {
        TEITensor transformed(tei.dimensions());
        Eigen::Index dim = tei.dimension(0);

        for (Eigen::Index mu = 0; mu < dim; ++mu)
            for (Eigen::Index nu = 0; nu < dim; ++nu)
                for (Eigen::Index kap = 0; kap < dim; ++kap)
                    for (Eigen::Index lam = 0; lam < dim; ++lam)
                        transformed(mu, nu, kap, lam) = tei(mu, kap, nu, lam);

        return transformed;
    }

    Eigen::Vector3d IntegralHelper::GaussianCenter(double exponent1, double exponent2,
                                                   const Eigen::Vector3d &pos1, const Eigen::Vector3d &pos2) {
        return (exponent1 * pos1 + exponent2 * pos2) / (exponent1 + exponent2);
    }

    ConfigIntegrals IntegralHelper::CalculateConfiguration(const Molecule &molecule) {
        Eigen::MatrixXd S =
                Eigen::MatrixXd::Zero(molecule.GetNoOrbitals(), molecule.GetNoOrbitals());
        Eigen::MatrixXd T =
                Eigen::MatrixXd::Zero(molecule.GetNoOrbitals(), molecule.GetNoOrbitals());
        Eigen::MatrixXd V =
                Eigen::MatrixXd::Zero(molecule.GetNoOrbitals(), molecule.GetNoOrbitals());
        TEITensor tei;

        for (size_t i = 0; i < molecule.GetNoOrbitals(); i++) {
            for (size_t j = i; j < molecule.GetNoOrbitals(); j++) {
                const GaussianCollection &c1 = molecule.GetOrbital(i);
                const GaussianCollection &c2 = molecule.GetOrbital(j);

                //Overlap Evaluation
                S(i, j) = S(j, i) = EvaluateOverlap(c1, c2);

                //Kinetic Evaluation
                T(i, j) = T(j, i) = EvaluateKinetic(c1, c2);

                double Vn = 0.0;
                for (const Atom &a: molecule.GetAtoms()) {
                    Vn += EvaluateNuclear(c1, c2, a.GetPos(), a.GetCharge());
                }
                V(i, j) = V(j, i) = Vn;
            }
        }

        tei = CalculateTwoElectronTensor(molecule);

        return ConfigIntegrals{S, T, V, tei};
    }

    double IntegralHelper::Binomial(int a, int b) {
        if ((a < 0) || (b < 0) || (a - b < 0)) {
            return 1.0;
        }
        return boost::math::factorial<double>(a) /
               (boost::math::factorial<double>(b) * boost::math::factorial<double>(a - b));
    }

    double IntegralHelper::PolynomialCoeff(int s, int ia, int ib, double xpa, double xpb) {
        double result = 0.0;

        for (int t = 0; t < s + 1; t++) {
            if ((s - ia <= t) && (t <= ib)) {
                result += Binomial(ia, s - t) *
                        Binomial(ib, t) *
                        std::pow(xpa, ia - s + t) *
                        std::pow(xpb, ib - t);
            }
        }

        return result;
    }

#pragma endregion

#pragma region Overlap

    double IntegralHelper::Overlap1D(double gamma,
                                     int angX1, double posX1,
                                     int angX2, double posX2) {
        // Only even values of i and j will persist
        // Note that(Px - Ax) terms reduce to zero if the posX coords are relative to the centre
        double result = 0.0;
        for (int i = 0; i < (1 + std::floor(0.5 * (angX1 + angX2))); i++) {
            result += PolynomialCoeff(2 * i, angX1, angX2, posX1, posX2) *
                    (i == 0 ? 1 : boost::math::double_factorial<double>(2 * i - 1)) /
                    std::pow(2 * gamma, i);
        }
        return result;
    }

    double IntegralHelper::OverlapGaussianParam(
        double exponent1, double exponent2,
        const Eigen::Vector3i &angular1, const Eigen::Vector3i &angular2,
        const Eigen::Vector3d &pos1, const Eigen::Vector3d &pos2) {
        double rab2 = (pos1 - pos2).squaredNorm();
        // Sum of the two gaussian exponents
        double gamma = exponent1 + exponent2;
        Eigen::Vector3d centre = GaussianCenter(exponent1, exponent2, pos1, pos2);

        double Eab = std::exp(-exponent1 * exponent2 * rab2 / gamma) * std::pow(std::numbers::pi / gamma, 1.5);
        double wx = Overlap1D(gamma, angular1[0], centre[0] - pos1[0],
                              angular2[0], centre[0] - pos2[0]);
        double wy = Overlap1D(gamma, angular1[1], centre[1] - pos1[1],
                              angular2[1], centre[1] - pos2[1]);
        double wz = Overlap1D(gamma, angular1[2], centre[2] - pos1[2],
                              angular2[2], centre[2] - pos2[2]);

        return Eab * wx * wy * wz;
    }

    double IntegralHelper::OverlapGaussian(const Gaussian &gto1, const Gaussian &gto2) {
        return OverlapGaussianParam(
            gto1.mExponent, gto2.mExponent,
            gto1.mAngular, gto2.mAngular,
            gto1.mPosition, gto2.mPosition);
    }

    double IntegralHelper::EvaluateOverlap(const GaussianCollection &c1, const GaussianCollection &c2) {
        double result = 0.0;

        for (size_t i = 0; i < c1.Size(); i++) {
            for (size_t j = 0; j < c2.Size(); j++) {
                const Gaussian &g1 = c1.Get(i);
                const Gaussian &g2 = c2.Get(j);
                result += g1.mNorm * g2.mNorm * g1.mCoeff *
                        g2.mCoeff * IntegralHelper::OverlapGaussian(g1, g2);
            }
        }

        return result;
    }

    Eigen::MatrixXd IntegralHelper::CalculateOverlapMatrix(const Molecule &molecule) {
        Eigen::MatrixXd S = Eigen::MatrixXd::Zero(molecule.GetOrbitals().size(), molecule.GetOrbitals().size());

        for (size_t i = 0; i < S.rows(); i++) {
            for (size_t j = i; j < S.cols(); j++) {
                S(i, j) = S(j, i) = EvaluateOverlap(molecule.GetOrbital(i), molecule.GetOrbital(j));
            }
        }

        return S;
    }

#pragma endregion

#pragma region Nuclear

    double IntegralHelper::AFactor(const int i, const int r, const int u, const int ang1, const int ang2,
                                   const double posA_P, double posB_P, double posNuc_P, double gamma) {
        double numerator = std::pow(-1, i) * PolynomialCoeff(i, ang1, ang2, posA_P, posB_P) * std::pow(-1, u) *
                           boost::math::factorial<double>(i) * std::pow(posNuc_P, i - 2 * r - 2 * u) * std::pow(
                               0.25 / gamma, r + u);
        double denominator = boost::math::factorial<double>(r) *
                             boost::math::factorial<double>(u) *
                             boost::math::factorial<double>(i - 2 * r - 2 * u);
        return numerator / denominator;
    }

    //https://csclub.uwaterloo.ca/~pbarfuss/jpsj.21.2313.pdf
    std::vector<double> IntegralHelper::GenerateACoeffs(int ang1,
                                                        int ang2, double posA_P, double posB_P, double posNuc_P,
                                                        double gamma) {
        int iRange = ang1 + ang2 + 1;
        std::vector<double> coeff(iRange, 0);

        for (size_t i = 0; i < iRange; i++) {
            for (size_t r = 0; r <= i / 2; r++) {
                for (size_t u = 0; u <= (i - 2 * r) / 2; u++) {
                    int index = i - 2 * r - u;
                    coeff[index] += AFactor(i, r, u, ang1, ang2, posA_P, posB_P, posNuc_P, gamma);
                }
            }
        }

        return coeff;
    }

    double IntegralHelper::NuclearGaussianParam(double exp1, double exp2, const Eigen::Vector3d &pos1,
                                                const Eigen::Vector3d &pos2,
                                                const Eigen::Vector3i &ang1, const Eigen::Vector3i &ang2,
                                                const Eigen::Vector3d &nuclearPos) {
        double gamma = exp1 + exp2;
        Eigen::Vector3d p = GaussianCenter(exp1, exp2, pos1, pos2);

        double r2 = (pos1 - pos2).squaredNorm();
        double rnuc2 = (nuclearPos - p).squaredNorm();

        std::vector<double> Ax = GenerateACoeffs(ang1[0], ang2[0], p[0] - pos1[0], p[0] - pos2[0], p[0] - nuclearPos[0],
                                                 gamma);
        std::vector<double> Ay = GenerateACoeffs(ang1[1], ang2[1], p[1] - pos1[1], p[1] - pos2[1], p[1] - nuclearPos[1],
                                                 gamma);
        std::vector<double> Az = GenerateACoeffs(ang1[2], ang2[2], p[2] - pos1[2], p[2] - pos2[2], p[2] - nuclearPos[2],
                                                 gamma);

        double result = 0.0;

        for (size_t i = 0; i <= ang1[0] + ang2[0]; i++) {
            for (size_t j = 0; j <= ang1[1] + ang2[1]; j++) {
                for (size_t k = 0; k <= ang1[2] + ang2[2]; k++) {
                    result += Ax[i] * Ay[j] * Az[k] * GammaFuncs::Fgamma(i + j + k, rnuc2 * gamma);
                }
            }
        }

        return -2.0 * std::numbers::pi / gamma * std::exp(-exp1 * exp2 * r2 / gamma) * result;
    }


    double IntegralHelper::NuclearGaussian(const Gaussian &gto1, const Gaussian &gto2,
                                           const Eigen::Vector3d &nuclearPos) {
        return NuclearGaussianParam(gto1.mExponent, gto2.mExponent,
                                    gto1.mPosition, gto2.mPosition, gto1.mAngular, gto2.mAngular, nuclearPos);
    }

    double IntegralHelper::EvaluateNuclear(const GaussianCollection &c1, const GaussianCollection &c2,
                                           const Eigen::Vector3d &nuclearPos, const int charge) {
        double result = 0.0;

        for (size_t i = 0; i < c1.Size(); i++) {
            for (size_t j = 0; j < c2.Size(); j++) {
                const Gaussian &g1 = c1.Get(i);
                const Gaussian &g2 = c2.Get(j);
                result += g1.mNorm * g2.mNorm * g1.mCoeff *
                        g2.mCoeff * IntegralHelper::NuclearGaussian(g1, g2, nuclearPos);
            }
        }

        return result * (double) charge;
    }

    Eigen::MatrixXd IntegralHelper::CalculateNuclearMatrix(const Molecule &molecule) {
        Eigen::MatrixXd V = Eigen::MatrixXd::Zero(molecule.GetNoOrbitals(), molecule.GetNoOrbitals());

        for (size_t i = 0; i < V.rows(); i++) {
            for (size_t j = i; j < V.rows(); j++) {
                double result = 0.0;
                for (const Atom &a: molecule.GetAtoms()) {
                    result += EvaluateNuclear(molecule.GetOrbital(i), molecule.GetOrbital(j),
                                              a.GetPos(), a.GetCharge());
                }
                V(i, j) = V(j, i) = result;
            }
        }

        return V;
    }
#pragma endregion

#pragma region Kinetic

    double IntegralHelper::KineticGaussian(const Gaussian &gto1, const Gaussian &gto2) {
        //double result = 0.0;
        //Eigen::Vector3i newAngular1;
        //Eigen::Vector3i newAngular2;

        ////Loop for x, y and z
        //for (unsigned char i = 0; i < 3; i++)
        //{
        //    newAngular1 = gto1.mAngular;
        //    newAngular2 = gto2.mAngular;

        //    newAngular1[i] += -1;
        //    newAngular2[i] += -1;

        //    result += 0.5 * gto1.mAngular[i] * gto2.mAngular[i] * OverlapGaussianParam(
        //        gto1.mExponent, gto2.mExponent,
        //        newAngular1, newAngular2, gto1.mPosition, gto2.mPosition);

        //    newAngular1 = gto1.mAngular;
        //    newAngular2 = gto2.mAngular;
        //    newAngular1[i] += 1;
        //    newAngular2[i] += 1;

        //    result += 2 * gto1.mExponent * gto2.mExponent * OverlapGaussianParam(
        //        gto1.mExponent, gto2.mExponent,
        //        newAngular1, newAngular2, gto1.mPosition, gto2.mPosition);

        //    newAngular1 = gto1.mAngular;
        //    newAngular2 = gto2.mAngular;
        //    newAngular1[i] += 1;
        //    newAngular2[i] += -1;

        //    result -= gto1.mExponent * gto2.mAngular[i] * OverlapGaussianParam(
        //        gto1.mExponent, gto2.mExponent,
        //        newAngular1, newAngular2, gto1.mPosition, gto2.mPosition);

        //    newAngular1 = gto1.mAngular;
        //    newAngular2 = gto2.mAngular;
        //    newAngular1[i] += -1;
        //    newAngular2[i] += 1;

        //    result -= gto2.mExponent * gto1.mAngular[i] * OverlapGaussianParam(
        //        gto1.mExponent, gto2.mExponent,
        //        newAngular1, newAngular2, gto1.mPosition, gto2.mPosition);
        //}

        //return result;

        double term0 = gto2.mExponent *
                       (2.0 * (gto2.mAngular[0] + gto2.mAngular[1] + gto2.mAngular[2]) + 3.0) *
                       OverlapGaussian(gto1, gto2);

        double term1 = -2.0 * std::pow(gto2.mExponent, 2.0) * (
                           OverlapGaussianParam(gto1.mExponent, gto2.mExponent, gto1.mAngular,
                                                gto2.mAngular + Eigen::Vector3i(2, 0, 0), gto1.mPosition,
                                                gto2.mPosition) +
                           OverlapGaussianParam(gto1.mExponent, gto2.mExponent, gto1.mAngular,
                                                gto2.mAngular + Eigen::Vector3i(0, 2, 0), gto1.mPosition,
                                                gto2.mPosition) +
                           OverlapGaussianParam(gto1.mExponent, gto2.mExponent, gto1.mAngular,
                                                gto2.mAngular + Eigen::Vector3i(0, 0, 2), gto1.mPosition,
                                                gto2.mPosition)
                       );
        double term2 = -0.5 * (gto2.mAngular[0] * (gto2.mAngular[0] - 1) *
                               OverlapGaussianParam(gto1.mExponent, gto2.mExponent, gto1.mAngular,
                                                    gto2.mAngular - Eigen::Vector3i(2, 0, 0), gto1.mPosition,
                                                    gto2.mPosition) +
                               gto2.mAngular[1] * (gto2.mAngular[1] - 1) *
                               OverlapGaussianParam(gto1.mExponent, gto2.mExponent, gto1.mAngular,
                                                    gto2.mAngular - Eigen::Vector3i(0, 2, 0), gto1.mPosition,
                                                    gto2.mPosition) +
                               gto2.mAngular[2] * (gto2.mAngular[2] - 1) *
                               OverlapGaussianParam(gto1.mExponent, gto2.mExponent, gto1.mAngular,
                                                    gto2.mAngular - Eigen::Vector3i(0, 0, 2), gto1.mPosition,
                                                    gto2.mPosition));

        return term0 + term1 + term2;
    }

    double IntegralHelper::EvaluateKinetic(const GaussianCollection &c1, const GaussianCollection &c2) {
        double result = 0.0;

        for (size_t i = 0; i < c1.Size(); i++) {
            for (size_t j = 0; j < c2.Size(); j++) {
                const Gaussian &g1 = c1.Get(i);
                const Gaussian &g2 = c2.Get(j);
                result += g1.mNorm * g2.mNorm * g1.mCoeff *
                        g2.mCoeff * IntegralHelper::KineticGaussian(g1, g2);
            }
        }

        return result;
    }

    Eigen::MatrixXd IntegralHelper::CalculateKineticMatrix(const Molecule &molecule) {
        Eigen::MatrixXd K = Eigen::MatrixXd::Zero(molecule.GetNoOrbitals(), molecule.GetNoOrbitals());

        //i=Rows, j=Cols
        for (size_t i = 0; i < K.rows(); i++) {
            for (size_t j = i; j < K.cols(); j++) {
                K(i, j) = K(j, i) = EvaluateKinetic(molecule.GetOrbital(i), molecule.GetOrbital(j));
            }
        }

        return K;
    }

#pragma endregion

#pragma region Two Electron Integrals

    TEITensor IntegralHelper::CalculateTwoElectronTensor(const Molecule &molecule) {
        TEITensor returnTensor((Eigen::Index) molecule.GetNoOrbitals(),
                               (Eigen::Index) molecule.GetNoOrbitals(),
                               (Eigen::Index) molecule.GetNoOrbitals(),
                               (Eigen::Index) molecule.GetNoOrbitals());

        for (Eigen::Index i = 0; i < molecule.GetNoOrbitals(); i++) {
            for (Eigen::Index j = 0; j < molecule.GetNoOrbitals(); j++) {
                for (Eigen::Index k = 0; k < molecule.GetNoOrbitals(); k++) {
                    for (Eigen::Index l = 0; l < molecule.GetNoOrbitals(); l++) {
                        returnTensor(i, j, k, l) = EvaluateTwoElectron(
                            molecule.GetOrbital(i),
                            molecule.GetOrbital(j),
                            molecule.GetOrbital(k),
                            molecule.GetOrbital(l));
                    }
                }
            }
        }

        return returnTensor;
    }

    double IntegralHelper::EvaluateTwoElectron(const GaussianCollection &c1, const GaussianCollection &c2,
                                               const GaussianCollection &c3, const GaussianCollection &c4) {
        double result = 0.0;

        for (unsigned int i = 0; i < c1.Size(); i++) {
            for (unsigned int j = 0; j < c2.Size(); j++) {
                for (unsigned int k = 0; k < c3.Size(); k++) {
                    for (unsigned int l = 0; l < c4.Size(); l++) {
                        const Gaussian &g1 = c1.Get(i);
                        const Gaussian &g2 = c2.Get(j);
                        const Gaussian &g3 = c3.Get(k);
                        const Gaussian &g4 = c4.Get(l);

                        result += g1.mNorm * g2.mNorm * g3.mNorm * g4.mNorm *
                                g1.mCoeff * g2.mCoeff * g3.mCoeff * g4.mCoeff * TwoElectronGaussian(g1, g2, g3, g4);
                    }
                }
            }
        }

        return result;
    }

    double IntegralHelper::TwoElectronGaussian(const Gaussian &gto1, const Gaussian &gto2,
                                               const Gaussian &gto3, const Gaussian &gto4) {
        return TwoElectronGaussian(
            gto1.mPosition, gto1.mAngular[0], gto1.mAngular[1], gto1.mAngular[2], gto1.mExponent,
            gto2.mPosition, gto2.mAngular[0], gto2.mAngular[1], gto2.mAngular[2], gto2.mExponent,
            gto3.mPosition, gto3.mAngular[0], gto3.mAngular[1], gto3.mAngular[2], gto3.mExponent,
            gto4.mPosition, gto4.mAngular[0], gto4.mAngular[1], gto4.mAngular[2], gto4.mExponent);
    }

    double IntegralHelper::TwoElectronGaussian(
        const Eigen::Vector3d &a, const int la, const int ma, const int na, const double alphaa,
        const Eigen::Vector3d &b, const int lb, const int mb, const int nb, const double alphab,
        const Eigen::Vector3d &c, const int lc, const int mc, const int nc, const double alphac,
        const Eigen::Vector3d &d, const int ld, const int md, const int nd, const double alphad) {
        static const double pi = boost::math::constants::pi<double>();
        double rab2 = (a - b).squaredNorm();
        double rcd2 = (c - d).squaredNorm();

        Eigen::Vector3d p = GaussianCenter(alphaa, alphab, a, b);
        Eigen::Vector3d q = GaussianCenter(alphac, alphad, c, d);

        double rpq2 = (p - q).squaredNorm();

        double gamma1 = alphaa + alphab;
        double gamma2 = alphac + alphad;
        double delta = 0.25 * (1.0 / gamma1 + 1.0 / gamma2);

        std::vector<double> bx = B_array(la, lb, lc, ld, p[0], a[0], b[0], q[0], c[0], d[0], gamma1, gamma2, delta);
        std::vector<double> by = B_array(ma, mb, mc, md, p[1], a[1], b[1], q[1], c[1], d[1], gamma1, gamma2, delta);
        std::vector<double> bz = B_array(na, nb, nc, nd, p[2], a[2], b[2], q[2], c[2], d[2], gamma1, gamma2, delta);

        double sum = 0.0;
        for (int i = 0; i <= (la + lb + lc + ld); i++) {
            for (int j = 0; j <= (ma + mb + mc + md); j++) {
                for (int k = 0; k <= (na + nb + nc + nd); k++) {
                    sum += bx[i] * by[j] * bz[k] * GammaFuncs::Fgamma(i + j + k, 0.25 * rpq2 / delta);
                }
            }
        }

        return 2.0 * std::pow(pi, 2.5) / (gamma1 * gamma2 * std::sqrt(gamma1 + gamma2)) *
               std::exp(-alphaa * alphab * rab2 / gamma1) *
               std::exp(-alphac * alphad * rcd2 / gamma2) * sum;
    }

    std::vector<double> IntegralHelper::B_array(const int l1, const int l2, const int l3, const int l4,
                                                const double p, const double a, const double b, const double q,
                                                const double c, const double d,
                                                const double g1, const double g2, const double delta) {
        int imax = l1 + l2 + l3 + l4 + 1;
        std::vector<double> arrB(imax, 0);

        for (int i1 = 0; i1 < l1 + l2 + 1; i1++) {
            for (int i2 = 0; i2 < l3 + l4 + 1; i2++) {
                for (int r1 = 0; r1 < i1 / 2 + 1; r1++) {
                    for (int r2 = 0; r2 < i2 / 2 + 1; r2++) {
                        for (int u = 0; u < (i1 + i2) / 2 - r1 - r2 + 1; u++) {
                            int i = i1 + i2 - 2 * (r1 + r2) - u;
                            arrB[i] += B_term(i1, i2, r1, r2, u, l1, l2, l3, l4,
                                              p, a, b, q, c, d, g1, g2, delta);
                        }
                    }
                }
            }
        }
        return arrB;
    }

    double IntegralHelper::B_term(const int i1, const int i2, const int r1, const int r2, const int u, const int l1,
                                  const int l2, const int l3, const int l4,
                                  const double px, const double ax, const double bx, const double qx, const double cx,
                                  const double dx, const double gamma1,
                                  const double gamma2, const double delta) {
        return fB(i1, l1, l2, px, ax, bx, r1, gamma1) *
               pow(-1, i2) * fB(i2, l3, l4, qx, cx, dx, r2, gamma2) *
               pow(-1, u) * fact_ratio2(i1 + i2 - 2 * (r1 + r2), u) *
               pow(qx - px, i1 + i2 - 2 * (r1 + r2) - 2 * u) /
               pow(delta, i1 + i2 - 2 * (r1 + r2) - u);
    }

    double IntegralHelper::fB(const int i, const int l1, const int l2, const double p, const double a, const double b,
                              const int r, const double g) {
        return PolynomialCoeff(i, l1, l2, p - a, p - b) * B0(i, r, g);
    }

    double IntegralHelper::B0(int i, int r, double g) {
        return fact_ratio2(i, r) * pow(4 * g, r - i);
    }

    double IntegralHelper::fact_ratio2(const int a, const int b) {
        return boost::math::factorial<double>(a) / boost::math::factorial<double>(b) / boost::math::factorial<
                   double>(a - 2 * b);
    }

#pragma endregion
}
