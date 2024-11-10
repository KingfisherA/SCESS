//Class for performing integral calculations
//Code modified from https://github.com/ifilot/pyqint?tab=readme-ov-file

#include<pch.h>
#include<libint2.hpp>

#include"Integrals.h"
#include "Quantum/Basis/Gaussian.h"
#include "Quantum/Molecule.h"

namespace SCESS {
    void Integrals::InitLibInt() {
        libint2::initialize(true);
    }

    void Integrals::CloseLibInt() {
        libint2::finalize();
    }


#pragma region Transforms
    //
    //
    TEITensor Integrals::Transforms::ToMolecularTEI(const TEITensor &tei, const Eigen::MatrixXd &C) {
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

    TEITensor Integrals::Transforms::AntisymmetriseTEI(const TEITensor &tei) {
        TEITensor transformed(tei.dimensions());

        for (Eigen::Index p = 0; p < transformed.dimension(0); ++p)
            for (Eigen::Index q = 0; q < transformed.dimension(1); ++q)
                for (Eigen::Index r = 0; r < transformed.dimension(2); ++r)
                    for (Eigen::Index s = 0; s < transformed.dimension(3); ++s)
                        transformed(p, q, r, s) = tei(p, r, q, s) - tei(p, s, q, r);

        return transformed;
    }

    TEITensor Integrals::Transforms::AntisymmetriseIntegrateOutSpinClosedShell(const TEITensor &moTEI) {
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
                        //One if the spin functions commute. Zero otherwise
                        Eigen::Index exchangeMask = (p % 2 == s % 2) * (q % 2 == r % 2);

                        //Convert from chemist's to physicists' notation and apply spin integration
                        transformed(p, q, r, s) =
                                moTEI(p / 2, r / 2, q / 2, s / 2) * static_cast<double>(coloumbMask)
                                - moTEI(p / 2, s / 2, q / 2, r / 2) * static_cast<double>(exchangeMask);
                    }

        return transformed;
    }

    TEITensor Integrals::Transforms::ConvertChemistsToPhysicists(const TEITensor &tei) {
        TEITensor transformed(tei.dimensions());
        Eigen::Index dim = tei.dimension(0);

        for (Eigen::Index mu = 0; mu < dim; ++mu)
            for (Eigen::Index nu = 0; nu < dim; ++nu)
                for (Eigen::Index kap = 0; kap < dim; ++kap)
                    for (Eigen::Index lam = 0; lam < dim; ++lam)
                        transformed(mu, nu, kap, lam) = tei(mu, kap, nu, lam);

        return transformed;
    }

#pragma endregion

#pragma region Helper Functions

    ConfigIntegrals Integrals::CalculateConfiguration(const Molecule &molecule) {
        Eigen::MatrixXd S =
                Eigen::MatrixXd::Zero(molecule.GetNoOrbitals(), molecule.GetNoOrbitals());
        Eigen::MatrixXd T =
                Eigen::MatrixXd::Zero(molecule.GetNoOrbitals(), molecule.GetNoOrbitals());
        Eigen::MatrixXd V =
                Eigen::MatrixXd::Zero(molecule.GetNoOrbitals(), molecule.GetNoOrbitals());

        AArray ax, ay, az;
        for (size_t i = 0; i < molecule.GetNoOrbitals(); i++)
            for (size_t j = i; j < molecule.GetNoOrbitals(); j++) {
                const ContractedGaussian &c1 = molecule.GetOrbital(i);
                const ContractedGaussian &c2 = molecule.GetOrbital(j);

                //Overlap Evaluation
                S(i, j) = S(j, i) = EvaluateOverlap(c1, c2);

                //Kinetic Evaluation
                T(i, j) = T(j, i) = EvaluateKinetic(c1, c2);

                double Vn = 0.0;
                for (const Atom &a: molecule.GetAtoms()) {
                    Vn += EvaluateNuclear(c1, c2, a.GetPos(), a.GetCharge(), ax, ay, az);
                }
                V(i, j) = V(j, i) = Vn;
            }


        TEITensor tei = TwoElectronRepulsionTensorMultithread(molecule);

        return ConfigIntegrals{S, T, V, tei};
    }

#pragma endregion

#pragma region Overlap

    double Integrals::Overlap1D(double gamma,
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

    double Integrals::OverlapGaussianParam(
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

    double Integrals::OverlapGaussian(const Gaussian &gto1, const Gaussian &gto2) {
        return OverlapGaussianParam(
            gto1.mExponent, gto2.mExponent,
            gto1.mAngular, gto2.mAngular,
            gto1.mPosition, gto2.mPosition);
    }

    double Integrals::EvaluateOverlap(const ContractedGaussian &c1, const ContractedGaussian &c2) {
        double result = 0.0;

        for (size_t i = 0; i < c1.Size(); i++) {
            for (size_t j = 0; j < c2.Size(); j++) {
                const Gaussian &g1 = c1.Get(i);
                const Gaussian &g2 = c2.Get(j);
                result += g1.mNorm * g2.mNorm * g1.mCoeff *
                        g2.mCoeff * Integrals::OverlapGaussian(g1, g2);
            }
        }

        return result;
    }

    Eigen::MatrixXd Integrals::CalculateOverlapMatrix(const Molecule &molecule) {
        Eigen::MatrixXd S = Eigen::MatrixXd::Zero(molecule.GetOrbitals().size(), molecule.GetOrbitals().size());

        for (Eigen::Index i = 0; i < S.rows(); i++)
            for (Eigen::Index j = i; j < S.cols(); j++)
                S(i, j) = S(j, i) = EvaluateOverlap(molecule.GetOrbital(i), molecule.GetOrbital(j));

        return S;
    }


    Eigen::MatrixXd Integrals::CalculateMixedOverlapMatrix(const Molecule &mol1, const Molecule &mol2) {
        Eigen::MatrixXd S = Eigen::MatrixXd::Zero(mol1.GetOrbitals().size(), mol2.GetOrbitals().size());

        for (Eigen::Index i = 0; i < S.rows(); i++)
            for (Eigen::Index j = i; j < S.cols(); j++)
                S(i, j) = S(j, i) = EvaluateOverlap(mol1.GetOrbital(i), mol2.GetOrbital(j));

        return S;
    }

#pragma endregion

#pragma region Nuclear

    void Integrals::A_Array(AArray &aArray, int ang1, int ang2, double posA_P, double posB_P, double posNuc_P,
                            double gamma) {
        auto AFactor = [](int i, int r, int u, int ang1, int ang2,
                          double posA_P, double posB_P, double posNuc_P, double gamma)-> double {
            double numerator = powMinus1(i) * PolynomialCoeff(i, ang1, ang2, posA_P, posB_P) * powMinus1(u) *
                               Factorial(i) * pow(posNuc_P, i - 2 * r - 2 * u) * std::pow(
                                   0.25 / gamma, r + u);
            double denominator = Factorial(r) *
                                 Factorial(u) *
                                 Factorial(i - 2 * r - 2 * u);
            return numerator / denominator;
        };

        for (size_t i = 0; i < aArray.size(); i++)
            for (size_t r = 0; r <= i / 2; r++)
                for (size_t u = 0; u <= (i - 2 * r) / 2; u++) {
                    int index = i - 2 * r - u;
                    aArray[index] += AFactor(i, r, u, ang1, ang2, posA_P, posB_P, posNuc_P, gamma);
                }
    }

    double Integrals::NuclearGaussian(const Gaussian &gto1, const Gaussian &gto2,
                                      const Eigen::Vector3d &nuclearPos,
                                      AArray &ax, AArray &ay, AArray &az) {
        using std::numbers::pi;

        const auto &exp1 = gto1.mExponent;
        const auto &exp2 = gto2.mExponent;
        const auto &pos1 = gto1.mPosition;
        const auto &pos2 = gto2.mPosition;
        const auto &ang1 = gto1.mAngular;
        const auto &ang2 = gto2.mAngular;

        double gamma = exp1 + exp2;
        Eigen::Vector3d p = GaussianCenter(exp1, exp2, pos1, pos2);

        double r2 = (gto1.mPosition - gto2.mPosition).squaredNorm();
        double rnuc2 = (nuclearPos - p).squaredNorm();

        ax.clear();
        ay.clear();
        az.clear();
        int lMax = ang1[0] + ang2[0] + 1;
        int mMax = ang1[1] + ang2[1] + 1;
        int nMax = ang1[2] + ang2[2] + 1;
        ax.reserve(std::max(100, lMax));
        ay.reserve(std::max(100, mMax));
        az.reserve(std::max(100, nMax));
        ax.resize(lMax, 0.0);
        ay.resize(mMax, 0.0);
        az.resize(nMax, 0.0);

        A_Array(ax, ang1[0], ang2[0], p[0] - pos1[0], p[0] - pos2[0], p[0] - nuclearPos[0],
                gamma);
        A_Array(ay, ang1[1], ang2[1], p[1] - pos1[1], p[1] - pos2[1], p[1] - nuclearPos[1],
                gamma);
        A_Array(az, ang1[2], ang2[2], p[2] - pos1[2], p[2] - pos2[2], p[2] - nuclearPos[2],
                gamma);

        double result = 0.0;

        for (size_t i = 0; i <= ang1[0] + ang2[0]; i++)
            for (size_t j = 0; j <= ang1[1] + ang2[1]; j++)
                for (size_t k = 0; k <= ang1[2] + ang2[2]; k++) {
                    result += ax[i] * ay[j] * az[k] * BoysFunc(i + j + k, rnuc2 * gamma);
                }

        return -2.0 * pi / gamma * exp(-gto1.mExponent * gto2.mExponent * r2 / gamma) * result;
    }

    double Integrals::EvaluateNuclear(const ContractedGaussian &c1, const ContractedGaussian &c2,
                                      const Eigen::Vector3d &nuclearPos, const int charge,
                                      AArray &ax, AArray &ay, AArray &az) {
        double result = 0.0;

        for (size_t i = 0; i < c1.Size(); i++) {
            for (size_t j = 0; j < c2.Size(); j++) {
                const Gaussian &g1 = c1.Get(i);
                const Gaussian &g2 = c2.Get(j);
                result += g1.mNorm * g2.mNorm * g1.mCoeff *
                        g2.mCoeff * Integrals::NuclearGaussian(g1, g2, nuclearPos, ax, ay, az);
            }
        }

        return result * (double) charge;
    }

    Eigen::MatrixXd Integrals::CalculateNuclearMatrix(const Molecule &molecule) {
        Eigen::MatrixXd V = Eigen::MatrixXd::Zero(molecule.GetNoOrbitals(), molecule.GetNoOrbitals());

        AArray ax, ay, az;
        for (size_t i = 0; i < V.rows(); i++)
            for (size_t j = i; j < V.rows(); j++) {
                double result = 0.0;
                for (const Atom &a: molecule.GetAtoms()) {
                    result += EvaluateNuclear(molecule.GetOrbital(i), molecule.GetOrbital(j),
                                              a.GetPos(), a.GetCharge(), ax, ay, az);
                }
                V(i, j) = V(j, i) = result;
            }

        return V;
    }
#pragma endregion

#pragma region Kinetic

    double Integrals::KineticGaussian(const Gaussian &gto1, const Gaussian &gto2) {
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

    double Integrals::EvaluateKinetic(const ContractedGaussian &c1, const ContractedGaussian &c2) {
        double result = 0.0;

        for (size_t i = 0; i < c1.Size(); i++) {
            for (size_t j = 0; j < c2.Size(); j++) {
                const Gaussian &g1 = c1.Get(i);
                const Gaussian &g2 = c2.Get(j);
                result += g1.mNorm * g2.mNorm * g1.mCoeff *
                        g2.mCoeff * Integrals::KineticGaussian(g1, g2);
            }
        }

        return result;
    }

    Eigen::MatrixXd Integrals::CalculateKineticMatrix(const Molecule &molecule) {
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

    TEITensor Integrals::TwoElectronRepulsionTensor(const Molecule &molecule) {
        TEITensor returnTensor((Eigen::Index) molecule.GetNoOrbitals(),
                               (Eigen::Index) molecule.GetNoOrbitals(),
                               (Eigen::Index) molecule.GetNoOrbitals(),
                               (Eigen::Index) molecule.GetNoOrbitals());

        std::vector<double> bx;
        std::vector<double> by;
        std::vector<double> bz;

        for (Eigen::Index i = 0; i < molecule.GetNoOrbitals(); i++) {
            for (Eigen::Index j = 0; j < molecule.GetNoOrbitals(); j++) {
                for (Eigen::Index k = 0; k < molecule.GetNoOrbitals(); k++) {
                    for (Eigen::Index l = 0; l < molecule.GetNoOrbitals(); l++) {
                        returnTensor(i, j, k, l) = TwoElectronRepulsion(
                            molecule.GetOrbital(i),
                            molecule.GetOrbital(j),
                            molecule.GetOrbital(k),
                            molecule.GetOrbital(l),
                            bx, by, bz);
                    }
                }
            }
        }

        return returnTensor;
    }

    TEITensor Integrals::TwoElectronRepulsionTensorMultithread(const Molecule &molecule) {
        std::unordered_map<size_t, std::vector<Eigen::Vector4i> > yoshIdxMap;

        for (int i = 0; i < molecule.GetNoOrbitals(); i++)
            for (int j = 0; j < molecule.GetNoOrbitals(); j++)
                for (int k = 0; k < molecule.GetNoOrbitals(); k++)
                    for (int l = 0; l < molecule.GetNoOrbitals(); l++) {
                        size_t yoshIdx = TEIYoshimineIdx(i, j, k, l);

                        auto itPair = yoshIdxMap.find(yoshIdx);
                        std::vector<Eigen::Vector4i> *indxArray = nullptr;
                        if (itPair == yoshIdxMap.end()) {
                            indxArray = &yoshIdxMap.emplace(yoshIdx, std::vector<Eigen::Vector4i>{}).first->second;
                            indxArray->reserve(4);
                        } else {
                            indxArray = &itPair->second;
                        }

                        indxArray->emplace_back((int) i, (int) j, (int) k, (int) l);
                    }

        constexpr size_t NUM_THREADS = 3;
        std::thread threads[NUM_THREADS];
        std::unordered_map<size_t, double> integralResultMaps[NUM_THREADS + 1];

        auto threadWorker = [](const std::unordered_map<size_t, std::vector<Eigen::Vector4i> > &yoshIdxMap,
                               const Molecule &molecule, std::unordered_map<size_t, double> &integralResult,
                               size_t startingIdx,
                               size_t chunkSize) {
            auto it = yoshIdxMap.begin();
            for (size_t i = 0; i < startingIdx; i++)
                ++it;

            std::vector<double> bx;
            std::vector<double> by;
            std::vector<double> bz;

            for (size_t i = startingIdx; i < startingIdx + chunkSize; i++) {
                const Eigen::Vector4i &indices = it->second[0];
                integralResult[it->first] = TwoElectronRepulsion(molecule.GetOrbital(indices.x()),
                                                                 molecule.GetOrbital(indices.y()),
                                                                 molecule.GetOrbital(indices.z()),
                                                                 molecule.GetOrbital(indices.w()),
                                                                 bx, by, bz);
                ++it;
            }
        };

        size_t chunkSize = yoshIdxMap.size() / (NUM_THREADS + 1);

        if (chunkSize >= 20) {
            for (uint32_t i = 0; i < NUM_THREADS; i++)
                threads[i] = std::thread(threadWorker, std::ref(yoshIdxMap), std::ref(molecule),
                                         std::ref(integralResultMaps[i + 1]),
                                         i * chunkSize, chunkSize);

            //Dispatch work on current thread as well
            threadWorker(yoshIdxMap, molecule, integralResultMaps[0], NUM_THREADS * chunkSize,
                         yoshIdxMap.size() - NUM_THREADS * chunkSize);

            for (auto &thread: threads)
                thread.join();
        } else
            threadWorker(yoshIdxMap, molecule, integralResultMaps[0], 0, yoshIdxMap.size());


        TEITensor teiTensor((Eigen::Index) molecule.GetNoOrbitals(),
                            (Eigen::Index) molecule.GetNoOrbitals(),
                            (Eigen::Index) molecule.GetNoOrbitals(),
                            (Eigen::Index) molecule.GetNoOrbitals());

        for (auto &[yoshIdx, indices]: yoshIdxMap) {
            double integralResult = [&yoshIdx, &integralResultMaps, NUM_THREADS]() {
                for (auto &integralResultMap: integralResultMaps)
                    if (integralResultMap.contains(yoshIdx))
                        return integralResultMap.at(yoshIdx);
                return 0.0;
            }();

            for (const Eigen::Vector4i &idxVec: indices)
                teiTensor(idxVec[0], idxVec[1], idxVec[2], idxVec[3]) = integralResult;
        }


        return teiTensor;
    }

    double Integrals::TwoElectronRepulsion(const ContractedGaussian &c1, const ContractedGaussian &c2,
                                           const ContractedGaussian &c3, const ContractedGaussian &c4,
                                           BArray &bx, BArray &by, BArray &bz) {
        double result = 0.0;

        for (size_t i = 0; i < c1.Size(); i++)
            for (size_t j = 0; j < c2.Size(); j++)
                for (size_t k = 0; k < c3.Size(); k++)
                    for (size_t l = 0; l < c4.Size(); l++) {
                        const Gaussian &g1 = c1.Get(i),
                                g2 = c2.Get(j),
                                g3 = c3.Get(k),
                                g4 = c4.Get(l);

                        result += g1.mNorm * g2.mNorm * g3.mNorm * g4.mNorm *
                                g1.mCoeff * g2.mCoeff * g3.mCoeff * g4.mCoeff * TwoElectronRepulsion(
                                    g1, g2, g3, g4, bx, by, bz);
                    }

        return result;
    }

    double Integrals::TwoElectronRepulsion(const Gaussian &gto1, const Gaussian &gto2,
                                           const Gaussian &gto3, const Gaussian &gto4,
                                           BArray &bx, BArray &by, BArray &bz) {
        constexpr double pi = std::numbers::pi_v<double>;

        const Eigen::Vector3d &a = gto1.mPosition;
        const int &la = gto1.mAngular[0], ma = gto1.mAngular[1], na = gto1.mAngular[2];
        const double &alphaa = gto1.mExponent;

        const Eigen::Vector3d &b = gto2.mPosition;
        const int &lb = gto2.mAngular[0], mb = gto2.mAngular[1], nb = gto2.mAngular[2];
        const double &alphab = gto2.mExponent;

        const Eigen::Vector3d &c = gto3.mPosition;
        const int &lc = gto3.mAngular[0], mc = gto3.mAngular[1], nc = gto3.mAngular[2];
        const double &alphac = gto3.mExponent;

        const Eigen::Vector3d &d = gto4.mPosition;
        const int &ld = gto4.mAngular[0], md = gto4.mAngular[1], nd = gto4.mAngular[2];
        const double &alphad = gto4.mExponent;

        double rab2 = (a - b).squaredNorm();
        double rcd2 = (c - d).squaredNorm();

        Eigen::Vector3d p = GaussianCenter(alphaa, alphab, a, b);
        Eigen::Vector3d q = GaussianCenter(alphac, alphad, c, d);

        double rpq2 = (p - q).squaredNorm();

        double gamma1 = alphaa + alphab;
        double gamma2 = alphac + alphad;
        double delta = 0.25 * (1.0 / gamma1 + 1.0 / gamma2);

        bx.clear();
        by.clear();
        bz.clear();
        int lMax = la + lb + lc + ld + 1;
        int mMax = ma + mb + mc + md + 1;
        int nMax = na + nb + nc + nd + 1;
        bx.reserve(std::max(100, lMax));
        by.reserve(std::max(100, mMax));
        bz.reserve(std::max(100, nMax));
        bx.resize(lMax, 0.0);
        by.resize(mMax, 0.0);
        bz.resize(nMax, 0.0);

        B_Array(bx, la, lb, lc, ld, p[0], a[0], b[0], q[0], c[0], d[0], gamma1, gamma2, delta);
        B_Array(by, ma, mb, mc, md, p[1], a[1], b[1], q[1], c[1], d[1], gamma1, gamma2, delta);
        B_Array(bz, na, nb, nc, nd, p[2], a[2], b[2], q[2], c[2], d[2], gamma1, gamma2, delta);

        // std::vector<double> bx = B_Array(la, lb, lc, ld, p[0], a[0], b[0], q[0], c[0], d[0], gamma1, gamma2, delta);
        // std::vector<double> by = B_Array(ma, mb, mc, md, p[1], a[1], b[1], q[1], c[1], d[1], gamma1, gamma2, delta);
        // std::vector<double> bz = B_Array(na, nb, nc, nd, p[2], a[2], b[2], q[2], c[2], d[2], gamma1, gamma2, delta);

        double sum = 0.0;
        for (int i = 0; i <= (la + lb + lc + ld); i++)
            for (int j = 0; j <= (ma + mb + mc + md); j++)
                for (int k = 0; k <= (na + nb + nc + nd); k++) {
                    // Boys(i + j + k, 0.25 * rpq2 / delta)
                    sum += bx[i] * by[j] * bz[k] * BoysFunc(i + j + k, 0.25 * rpq2 / delta);
                }

        return 2.0 * pow(pi, 2.5) / (gamma1 * gamma2 * sqrt(gamma1 + gamma2)) *
               exp(-alphaa * alphab * rab2 / gamma1) *
               exp(-alphac * alphad * rcd2 / gamma2) * sum;
    }

    void Integrals::B_Array(BArray &bArray, int l1, int l2, int l3, int l4,
                            double p, double a, double b, double q,
                            double c, double d,
                            double g1, double g2, double delta) {
        auto factRatio2 = [](int a, int b)-> int {
            return Factorial(a) / Factorial(b) /
                   Factorial(a - 2 * b);
        };

        auto B0 = [&factRatio2](int i, int r, double g)-> double {
            return factRatio2(i, r) * pow(4 * g, r - i);
        };

        auto fB = [&B0](int i, int l1, int l2, double p, double a, double b,
                        int r, double g)-> double {
            return PolynomialCoeff(i, l1, l2, p - a, p - b) * B0(i, r, g);
        };

        auto BArrayTerm = [&fB, &factRatio2](int i1, int i2, int r1, int r2, int u,
                                             int l1,
                                             int l2, int l3, int l4,
                                             double px, double ax, double bx, double qx,
                                             double cx,
                                             double dx, double gamma1,
                                             double gamma2, double delta)-> double {
            return fB(i1, l1, l2, px, ax, bx, r1, gamma1) *
                   powMinus1(i2) * fB(i2, l3, l4, qx, cx, dx, r2, gamma2) *
                   powMinus1(u) * factRatio2(i1 + i2 - 2 * (r1 + r2), u) *
                   pow(qx - px, i1 + i2 - 2 * (r1 + r2) - 2 * u) /
                   pow(delta, i1 + i2 - 2 * (r1 + r2) - u);
        };


        for (int i1 = 0; i1 < l1 + l2 + 1; i1++)
            for (int i2 = 0; i2 < l3 + l4 + 1; i2++)
                for (int r1 = 0; r1 < i1 / 2 + 1; r1++)
                    for (int r2 = 0; r2 < i2 / 2 + 1; r2++)
                        for (int u = 0; u < (i1 + i2) / 2 - r1 - r2 + 1; u++) {
                            int i = i1 + i2 - 2 * (r1 + r2) - u;
                            bArray[i] += BArrayTerm(i1, i2, r1, r2, u, l1, l2, l3, l4,
                                                    p, a, b, q, c, d, g1, g2, delta);
                        }
    }


#pragma endregion
    //
}
