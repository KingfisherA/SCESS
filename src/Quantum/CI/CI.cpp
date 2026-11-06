#include <array>
#include <cppitertools/combinations.hpp>
#include <boost/range/irange.hpp>
#include <boost/math/special_functions/binomial.hpp>
//
#include "CI.h"
#include "Quantum/Integrals/IntegralTransforms.h"
#include "Tensor/TensorParallel.h"
//
#include <Eigen/Sparse>
#include "Dahlia/Tools/Threading.h"


//
//
namespace SCESS::CI {
    //
    using namespace Eigen;
    //
    //
#pragma region Slater Determinant Code
    //
    //
    SlaterDet SlaterDet::MakeReference(uint16_t nAlphaOcc, uint16_t nBetaOcc) {
        SlaterDet reference;
        uint16_t i = 0;

        for (; i < nAlphaOcc; i++)
            reference.GetAlphaOrbitals() |= 1ULL << i;

        for (i = 0; i < nBetaOcc; i++)
            reference.GetBetaOrbitals() |= 1ULL << i;

        return reference;
    }
    //
    static uint16_t CountSetBits(uint64_t bitStr, uint32_t nBits) {
        uint32_t nBitsSet = 0;
        for (uint32_t i = 0; i < nBits; i++)
            if ((1ULL << i) & bitStr)
                nBitsSet++;

        return nBitsSet;
    }
    //
    void SlaterDet::CalculateDifferences(SlaterDet I,
                                         SlaterDet J,
                                         uint16_t numOfAlphaOrbitals,
                                         uint16_t numOfBetaOrbitals,
                                         uint16_t *uniqueIArray,
                                         uint16_t *uniqueISize,
                                         uint16_t *uniqueJArray,
                                         uint16_t *uniqueJSize,
                                         uint16_t *commonIJArray,
                                         uint16_t *commonIJSize,
                                         int8_t *sign) {
        //Arrays for uniqueIArray, uniqueJArray and commonIJArray must be at least
        //numOfAlphaOrbitals + numOfBetaOrbitals in size as they hold spin orbitals of both spins

        uint16_t nElectrons = CountSetBits(I.GetAlphaOrbitals(), numOfAlphaOrbitals) +
                              CountSetBits(I.GetBetaOrbitals(), numOfBetaOrbitals);
        uint16_t excitationLvl = CountSetBits(I.GetAlphaOrbitals() ^ J.GetAlphaOrbitals(), numOfAlphaOrbitals) +
                                 CountSetBits(I.GetBetaOrbitals() ^ J.GetBetaOrbitals(), numOfBetaOrbitals);

        uint16_t shift = nElectrons - excitationLvl;
        uint16_t flips = 0;

        uint16_t commonIJPtr = 0,
                uniqueIPtr = 0,
                uniqueJPtr = 0;

        auto CalcDifferences = [&flips, &shift,
                    &commonIJArray, &commonIJPtr,
                    &uniqueIArray, &uniqueIPtr,
                    &uniqueJArray, &uniqueJPtr](OrbBitStr I,
                                                OrbBitStr J,
                                                size_t numOfOrbs,
                                                uint8_t (*func)(uint8_t)) {
            uint32_t iel1 = 0, iel2 = 0, iexcit1 = 0, iexcit2 = 0;

            for (uint32_t i = 0; i < numOfOrbs; i++) {
                //Occupied in J
                if ((1ULL << i) & J) iel2++;

                //Occupied in I
                if ((1ULL << i) & I) {
                    iel1++;
                    //Occupied in J
                    if ((1ULL << i) & J) {
                        if (commonIJArray)
                            commonIJArray[commonIJPtr++] = func(i);
                        //commonIJArray.push_back(func(i));
                    }
                    //Not occupied in J
                    else {
                        iexcit1++;
                        flips += (shift - iel1 + iexcit1);
                        // uniqueIArray.push_back(func(i));
                        if (uniqueIArray)
                            uniqueIArray[uniqueIPtr++] = func(i);
                    }
                }
                //Not occupied in I
                else {
                    //Occupied in J
                    if ((1ULL << i) & J) {
                        iexcit2++;
                        //uniqueJArray.push_back(func(i));
                        flips += (shift - iel2 + iexcit2);
                        if (uniqueJArray)
                            uniqueJArray[uniqueJPtr++] = func(i);
                    }
                }
            }
        };


        CalcDifferences(I.GetAlphaOrbitals(), J.GetAlphaOrbitals(), numOfAlphaOrbitals,
                        [](uint8_t i)-> uint8_t { return i * 2; });
        CalcDifferences(I.GetBetaOrbitals(), J.GetBetaOrbitals(), numOfBetaOrbitals,
                        [](uint8_t i)-> uint8_t { return i * 2 + 1; });

        if (commonIJSize)
            *commonIJSize = commonIJPtr;
        if (uniqueISize)
            *uniqueISize = uniqueIPtr;
        if (uniqueJSize)
            *uniqueJSize = uniqueJPtr;

        if (sign) {
            *sign = 1;
            if (flips % 2 == 1)
                *sign = -1;
        }
    }
    //
    int8_t SlaterDet::CalculateSign(SlaterDet I, SlaterDet J,
                                    uint16_t numOfAlphaOrbitals, uint16_t numOfBetaOrbitals) {
        uint16_t nElectrons = CountSetBits(I.GetAlphaOrbitals(), numOfAlphaOrbitals) +
                              CountSetBits(I.GetBetaOrbitals(), numOfBetaOrbitals);
        uint16_t excitationLvl = CountSetBits(I.GetAlphaOrbitals() ^ J.GetAlphaOrbitals(), numOfAlphaOrbitals) +
                                 CountSetBits(I.GetBetaOrbitals() ^ J.GetBetaOrbitals(), numOfBetaOrbitals);

        uint16_t shift = nElectrons - excitationLvl;
        uint16_t flips = 0;

        auto CalcSign = [&flips, &shift](OrbBitStr I,
                                         OrbBitStr J,
                                         size_t numOfOrbs) {
            uint16_t iel1 = 0, iel2 = 0, iexcit1 = 0, iexcit2 = 0;

            for (uint16_t i = 0; i < numOfOrbs; i++) {
                //Occupied in J
                if ((1ULL << i) & J) iel2++;

                //Occupied in I
                if ((1ULL << i) & I) {
                    iel1++;
                    //Not occupied in J
                    if (!((1ULL << i) & J)) {
                        iexcit1++;
                        flips += (shift - iel1 + iexcit1);
                    }
                }
                //Not occupied in I
                else {
                    //Occupied in J
                    if ((1ULL << i) & J) {
                        iexcit2++;
                        flips += (shift - iel2 + iexcit2);
                    }
                }
            }
        };

        CalcSign(I.GetAlphaOrbitals(), J.GetAlphaOrbitals(), numOfAlphaOrbitals);
        CalcSign(I.GetBetaOrbitals(), J.GetBetaOrbitals(), numOfBetaOrbitals);

        int8_t sign = 1;
        if (flips % 2 == 1)
            sign = -1;

        return sign;
    }

    //
    //
#pragma endregion
    //
    //
#pragma region CI Utility Functions
    //
    //
    MatrixXd GenerateHCoreSpinTransform(
        const MatrixXd &h,
        const MatrixXd &C) {
        // Transform the one-electron Hamiltonian matrix into the MO basis set
        // and integrate out spin functions
        MatrixXd HSpin = MatrixXd::Zero(C.rows() * 2, C.cols() * 2);
        MatrixXd H_Transformed = C.transpose() * h * C;

        for (uint32_t i = 0; i < HSpin.rows(); i++)
            for (uint32_t j = i; j < HSpin.cols(); j++) {
                //Integrate out alpha and beta spin
                if ((i + j) % 2 != 0) {
                    //Alpha with beta spin function so zero
                    HSpin(i, j) = HSpin(j, i) = 0.0;
                    continue;
                }

                uint32_t ii = i / 2, jj = j / 2;

                HSpin(i, j) = HSpin(j, i) = H_Transformed(ii, jj);
            }

        return HSpin;
    }
    //
    MatrixXd GenerateHCoreSpin(const Eigen::MatrixXd &h) {
        MatrixXd HSpin = MatrixXd::Zero(h.rows() * 2, h.cols() * 2);

        for (uint32_t i = 0; i < HSpin.rows(); i++)
            for (uint32_t j = i; j < HSpin.cols(); j++) {
                //Integrate out alpha and beta spin
                if ((i + j) % 2 != 0) {
                    //Alpha with beta spin function so zero
                    HSpin(i, j) = HSpin(j, i) = 0.0;
                    continue;
                }

                uint32_t ii = i / 2, jj = j / 2;

                HSpin(i, j) = HSpin(j, i) = h(ii, jj);
            }

        return HSpin;
    }
    //
    double CalculateHamiltonianElement(SlaterDet detI, SlaterDet detJ,
                                       size_t numOfAlphaOrbitals,
                                       size_t numOfBetaOrbitals,
                                       const MatrixXd &HSpin,
                                       const TEITensor &spinERI) {
        return CalculateHamiltonianElement(detI, detJ, numOfAlphaOrbitals, numOfBetaOrbitals,
                                           DenseSpinOrbitalIntegrals{HSpin, spinERI});
    }
    //
    SpinOrbitalIntegrals MakeSpinOrbitalIntegrals(const HFWfn &rhfWfn, const Integrals::MolInts &mints) {
        const MatrixXd &C = rhfWfn.C_Alpha;
        return SpinOrbitalIntegrals{
            C.transpose() * (mints.T + mints.V) * C,
            Integrals::Transforms::TransformTEI(mints.eri, C)
        };
    }
    //
    std::vector<SlaterDet> ConnectedDeterminants(SlaterDet det, uint32_t nSpatialOrbitals) {
        const OrbBitStr mask = nSpatialOrbitals >= 64 ? ~0ULL : (1ULL << nSpatialOrbitals) - 1;
        auto Occupied = [nSpatialOrbitals](OrbBitStr str) {
            std::vector<uint32_t> orbs;
            for (uint32_t i = 0; i < nSpatialOrbitals; i++)
                if (str >> i & 1ULL)
                    orbs.push_back(i);
            return orbs;
        };
        const OrbBitStr alpha = det.GetAlphaOrbitals(), beta = det.GetBetaOrbitals();
        const std::vector<uint32_t> occA = Occupied(alpha), occB = Occupied(beta);
        const std::vector<uint32_t> virA = Occupied(~alpha & mask), virB = Occupied(~beta & mask);

        //Single and double replacements of one string
        auto SameSpin = [](OrbBitStr str, const std::vector<uint32_t> &occ, const std::vector<uint32_t> &vir,
                           std::vector<OrbBitStr> &singles, std::vector<OrbBitStr> &doubles) {
            for (uint32_t i: occ)
                for (uint32_t a: vir)
                    singles.push_back(str ^ (1ULL << i) ^ (1ULL << a));
            for (size_t i = 0; i < occ.size(); i++)
                for (size_t j = i + 1; j < occ.size(); j++)
                    for (size_t a = 0; a < vir.size(); a++)
                        for (size_t b = a + 1; b < vir.size(); b++)
                            doubles.push_back(str ^ (1ULL << occ[i]) ^ (1ULL << occ[j]) ^
                                              (1ULL << vir[a]) ^ (1ULL << vir[b]));
        };
        std::vector<OrbBitStr> singlesA, doublesA, singlesB, doublesB;
        SameSpin(alpha, occA, virA, singlesA, doublesA);
        SameSpin(beta, occB, virB, singlesB, doublesB);

        std::vector<SlaterDet> connected;
        connected.reserve(singlesA.size() + doublesA.size() + singlesB.size() + doublesB.size() +
                          singlesA.size() * singlesB.size());
        for (OrbBitStr a: singlesA)
            connected.emplace_back(a, beta);
        for (OrbBitStr a: doublesA)
            connected.emplace_back(a, beta);
        for (OrbBitStr b: singlesB)
            connected.emplace_back(alpha, b);
        for (OrbBitStr b: doublesB)
            connected.emplace_back(alpha, b);
        for (OrbBitStr a: singlesA)
            for (OrbBitStr b: singlesB)
                connected.emplace_back(a, b);
        return connected;
    }
    //
    OrbBitStr MixAlphaAndBeta(OrbBitStr alphaList,
                              size_t numAlpha,
                              OrbBitStr betaList,
                              size_t numBeta) {
        OrbBitStr combinedIndices = 0;
        size_t i = 0;

        for (i = 0; i < numAlpha; i++)
            combinedIndices |= (alphaList >> i) & 0x1;

        for (i = 0; i < numBeta; i++)
            combinedIndices |= (betaList >> i) & 0x1;

        return combinedIndices;
    }
    //
    std::vector<SlaterDet> EnumerateDet(uint32_t nSpinOrbs,
                                        uint32_t nAlphaElectrons,
                                        uint32_t nBetaElectrons,
                                        uint32_t excitorLevel) {
        //Even spin orbs are alpha spin
        //Odd spin orbs are beta spin

        std::vector<uint32_t> occupiedSpinOrbs;
        std::vector<uint32_t> virtualSpinOrbs;
        for (uint32_t i = 0; i < nSpinOrbs; i++) {
            const bool isAlpha = i % 2 == 0;
            if ((isAlpha && i / 2 < nAlphaElectrons) || (!isAlpha && i / 2 < nBetaElectrons))
                occupiedSpinOrbs.push_back(i);
            else
                virtualSpinOrbs.push_back(i);
        }

        if (excitorLevel > occupiedSpinOrbs.size())
            excitorLevel = occupiedSpinOrbs.size();

        if (excitorLevel > virtualSpinOrbs.size())
            excitorLevel = virtualSpinOrbs.size();

        auto CountAlpha = [](const auto &spinOrbs) {
            return std::ranges::count_if(spinOrbs, [](uint32_t i) { return i % 2 == 0; });
        };

        std::vector<SlaterDet> determinants;
        determinants.emplace_back(SlaterDet::MakeReference(nAlphaElectrons, nBetaElectrons));

        for (uint32_t ex = 1; ex <= excitorLevel; ex++)
            for (auto &&virtualSpinOrbComb: iter::combinations(virtualSpinOrbs, ex))
                for (auto &&occupiedSpinOrbComb: iter::combinations(occupiedSpinOrbs, ex)) {
                    //Only keep spin-conserving excitations (the same number of alpha and beta electrons as the
                    //reference) as the Hamiltonian cannot couple determinants of different M_s
                    if (CountAlpha(virtualSpinOrbComb) != CountAlpha(occupiedSpinOrbComb))
                        continue;

                    SlaterDet &det = determinants.emplace_back();
                    det = SlaterDet::MakeReference(nAlphaElectrons, nBetaElectrons);

                    for (const auto oSpinOrbIdx: occupiedSpinOrbComb) {
                        if (oSpinOrbIdx % 2 == 0) {
                            // If even then alpha spin
                            det.GetAlphaOrbitals() ^= 1ULL << (oSpinOrbIdx / 2);
                        } else {
                            // If odd then beta spin
                            det.GetBetaOrbitals() ^= 1ULL << (oSpinOrbIdx / 2);
                        }
                    }
                    for (const auto vSpinOrbIdx: virtualSpinOrbComb) {
                        if (vSpinOrbIdx % 2 == 0) {
                            // If even then alpha spin
                            det.GetAlphaOrbitals() |= 1ULL << (vSpinOrbIdx / 2);
                        } else {
                            // If odd then beta spin
                            det.GetBetaOrbitals() |= 1ULL << (vSpinOrbIdx / 2);
                        }
                    }
                }

        return determinants;
    }
    //
    std::vector<SlaterDet> EnumerateAllDet(uint32_t nSpinOrbs, uint32_t nElectrons) {
        auto nDeterminants = static_cast<size_t>(
            boost::math::binomial_coefficient<double>(nSpinOrbs, nElectrons));
        std::vector<SlaterDet> determinants;
        determinants.reserve(nDeterminants);

        for (auto &spinOrbCombination: iter::combinations(
                 boost::irange<uint8_t>(0, nSpinOrbs), nElectrons)) {
            SlaterDet &det = determinants.emplace_back();

            for (const auto spinOrbIndx: spinOrbCombination) {
                // If even then alpha spin
                if (spinOrbIndx % 2 == 0) {
                    det.GetAlphaOrbitals() |= 1ULL << (spinOrbIndx / 2);
                    continue;
                }
                // If odd then beta spin
                det.GetBetaOrbitals() |= 1ULL << (spinOrbIndx / 2);
            }
        }

        return determinants;
    }
    //
    std::vector<SlaterDet> EnumerateAllDet(uint32_t nSpinOrbs, uint32_t nAlphaElectrons,
                                           uint32_t nBetaElectrons) {
        const uint32_t nSpatialOrbs = nSpinOrbs / 2;

        auto MakeBitStrings = [nSpatialOrbs](uint32_t nElectrons) {
            std::vector<OrbBitStr> bitStrings;
            for (auto &&orbCombination: iter::combinations(boost::irange<uint32_t>(0, nSpatialOrbs), nElectrons)) {
                OrbBitStr bitStr = 0;
                for (const auto orbIdx: orbCombination)
                    bitStr |= 1ULL << orbIdx;
                bitStrings.push_back(bitStr);
            }
            return bitStrings;
        };

        const std::vector<OrbBitStr> alphaStrings = MakeBitStrings(nAlphaElectrons);
        const std::vector<OrbBitStr> betaStrings = MakeBitStrings(nBetaElectrons);

        //The first determinant is the reference as the lowest orbitals form the first combination
        std::vector<SlaterDet> determinants;
        determinants.reserve(alphaStrings.size() * betaStrings.size());
        for (const OrbBitStr alpha: alphaStrings)
            for (const OrbBitStr beta: betaStrings)
                determinants.emplace_back(alpha, beta);

        return determinants;
    }
    //
    size_t GetNumberOfOccupiedOrbitals(OrbBitStr orbitalBitStr,
                                       size_t numOfOrbitals) {
        size_t numOccupied = 0;

        for (size_t i = 0; i < numOfOrbitals; i++)
            if ((1ULL << i) & orbitalBitStr)
                numOccupied++;

        return numOccupied;
    }
    //
    void PrintDeterminants(const std::vector<SlaterDet> &determinants,
                           size_t numOfOrbitals) {
        for (auto det: determinants) {
            std::cout << "[Alpha, ";

            for (uint32_t i = 0; i < numOfOrbitals; i++)
                std::cout << (uint32_t) ((det.GetAlphaOrbitals() >> i) & 0x1) << ", ";

            std::cout << "Beta, ";

            for (uint32_t i = 0; i < numOfOrbitals; i++)
                std::cout << (uint32_t) ((det.GetBetaOrbitals() >> i) & 0x1) << ", ";


            std::cout << "] \n";
        }
    }
    //
    //
#pragma endregion
    //
    //
#pragma region Configuration Interaction Impl
    //
    //
    DavidsonResult Davidson(const std::function<MatrixXd(const MatrixXd &)> &sigma,
                            const VectorXd &diagonal,
                            const MatrixXd &guess,
                            const double residualTolerance,
                            const size_t maxIterations,
                            const DH::Thread *thread) {
        const Index n = diagonal.size();
        const Index nRoots = guess.cols();
        //The subspace collapses onto the current Ritz vectors when it would exceed this size
        const Index maxSubspace = std::min<Index>(n, std::max<Index>(4 * nRoots, 12));

        //Orthonormal basis of the subspace and H times it
        MatrixXd V = HouseholderQR<MatrixXd>(guess).householderQ() * MatrixXd::Identity(n, nRoots);
        MatrixXd S = sigma(V);

        DavidsonResult result;
        for (size_t iteration = 0; iteration < maxIterations; iteration++) {
            if (thread && thread->ShouldStop())
                return {};

            const MatrixXd subspaceH = V.transpose() * S;
            const SelfAdjointEigenSolver<MatrixXd> eig(0.5 * (subspaceH + subspaceH.transpose()));
            const VectorXd theta = eig.eigenvalues().head(nRoots);
            const MatrixXd Y = eig.eigenvectors().leftCols(nRoots);
            MatrixXd X = V * Y;
            MatrixXd SX = S * Y;
            const MatrixXd R = SX - X * theta.asDiagonal();

            result.eigenvalues = theta;
            result.eigenvectors = X;

            //Diagonally preconditioned corrections of the unconverged roots
            std::vector<VectorXd> corrections;
            for (Index k = 0; k < nRoots; k++) {
                if (R.col(k).norm() < residualTolerance)
                    continue;
                VectorXd t(n);
                for (Index i = 0; i < n; i++) {
                    double denominator = theta(k) - diagonal(i);
                    if (std::abs(denominator) < 1e-8)
                        denominator = denominator < 0.0 ? -1e-8 : 1e-8;
                    t(i) = R(i, k) / denominator;
                }
                corrections.push_back(std::move(t));
            }
            if (corrections.empty()) {
                result.converged = true;
                return result;
            }

            if (V.cols() + (Index) corrections.size() > maxSubspace) {
                V = std::move(X);
                S = std::move(SX);
            }

            //Orthogonalise the corrections against the subspace (twice for numerical stability) and append them
            std::vector<VectorXd> added;
            for (VectorXd &t: corrections) {
                for (int pass = 0; pass < 2; pass++) {
                    t -= V * (V.transpose() * t);
                    for (const VectorXd &a: added)
                        t -= a.dot(t) * a;
                }
                const double norm = t.norm();
                if (norm > 1e-10)
                    added.push_back(t / norm);
            }
            if (added.empty()) {
                result.converged = true;
                return result;
            }

            MatrixXd newV(n, (Index) added.size());
            for (Index k = 0; k < newV.cols(); k++)
                newV.col(k) = added[k];
            const MatrixXd newS = sigma(newV);

            const Index m = V.cols();
            V.conservativeResize(n, m + newV.cols());
            V.rightCols(newV.cols()) = newV;
            S.conservativeResize(n, m + newS.cols());
            S.rightCols(newS.cols()) = newS;
        }
        return result;
    }
    //
    namespace {
        //E_pq |string> = sign |target> (a+_p a_q) over the orbitals of one spin
        struct StringReplacement {
            uint32_t target;
            //p + n q, the row of (pq|rs) in the n^2 x n^2 integral matrix
            uint32_t pq;
            double sign;
        };
        //
        //Every string of nElectrons in nOrbitals orbitals, in the order of EnumerateAllDet, with every replacement
        //E_pq of each string (including p == q)
        struct StringSpace {
            std::vector<OrbBitStr> strings;
            std::vector<std::vector<StringReplacement> > replacements;

            StringSpace(uint32_t nOrbitals, uint32_t nElectrons) {
                for (auto &&combination: iter::combinations(boost::irange<uint32_t>(0, nOrbitals), nElectrons)) {
                    OrbBitStr str = 0;
                    for (const auto orb: combination)
                        str |= 1ULL << orb;
                    strings.push_back(str);
                }
                std::unordered_map<OrbBitStr, uint32_t> index;
                index.reserve(strings.size());
                for (uint32_t i = 0; i < strings.size(); i++)
                    index.emplace(strings[i], i);

                replacements.resize(strings.size());
                for (uint32_t i = 0; i < strings.size(); i++) {
                    const OrbBitStr str = strings[i];
                    for (uint32_t q = 0; q < nOrbitals; q++) {
                        if (!(str >> q & 1ULL))
                            continue;
                        for (uint32_t p = 0; p < nOrbitals; p++) {
                            if (p != q && (str >> p & 1ULL))
                                continue;
                            const OrbBitStr target = str ^ (1ULL << q) ^ (1ULL << p);
                            //The electrons between p and q which the moved electron passes
                            const uint32_t lo = std::min(p, q), hi = std::max(p, q);
                            const OrbBitStr between = hi - lo > 1 ? ((1ULL << hi) - 1) & ~((1ULL << (lo + 1)) - 1) : 0;
                            const double sign = std::popcount(str & between) % 2 ? -1.0 : 1.0;
                            replacements[i].push_back({index.at(target), p + nOrbitals * q, sign});
                        }
                    }
                }
            }
        };
        //
        //Full CI Hamiltonian products over the alpha x beta string products of a closed-shell determinant space
        //(Knowles & Handy, Chem. Phys. Lett. 111, 315 (1984); Olsen et al., J. Chem. Phys. 89, 2185 (1988)). A CI
        //vector is the matrix C(I_beta, I_alpha), the determinant order of EnumerateAllDet.
        class FCISigma {
        public:
            FCISigma(const SpinOrbitalIntegrals &ints, uint32_t nOrbitals, uint32_t nElectronsPerSpin)
                : mSpace(nOrbitals, nElectronsPerSpin), mN(nOrbitals) {
                const Index n = nOrbitals;
                mG = Map<const MatrixXd>(ints.eri.data(), n * n, n * n);
                //k_pq = h_pq - 1/2 sum_r (pr|rq)
                mK = ints.h;
                for (Index p = 0; p < n; p++)
                    for (Index q = 0; q < n; q++)
                        for (Index r = 0; r < n; r++)
                            mK(p, q) -= 0.5 * ints.eri(p, r, r, q);
            }

            [[nodiscard]] Index NumStrings() const { return (Index) mSpace.strings.size(); }

            MatrixXd operator()(const MatrixXd &vectors) const {
                MatrixXd result(vectors.rows(), vectors.cols());
                for (Index v = 0; v < vectors.cols(); v++)
                    result.col(v) = Apply(vectors.col(v));
                return result;
            }

        private:
            VectorXd Apply(const VectorXd &c) const {
                const Index nStr = NumStrings();
                const Map<const MatrixXd> C(c.data(), nStr, nStr);
                MatrixXd sigma = MatrixXd::Zero(nStr, nStr);

                //Alpha-alpha (columns are the alpha strings) and beta-beta, which is the same on the transpose as the
                //alpha and beta strings are the same for a closed-shell space
                SameSpin(C, sigma);
                const MatrixXd CT = C.transpose();
                MatrixXd sigmaT = MatrixXd::Zero(nStr, nStr);
                SameSpin(CT, sigmaT);
                sigma += sigmaT.transpose();
                OppositeSpin(C, sigma);

                return Map<const VectorXd>(sigma.data(), sigma.size());
            }

            //out(:, I) += sum_J <J| sum_pq k_pq E_pq + 1/2 sum_pqrs (pq|rs) E_pq E_rs |I> in(:, J), the excitations
            //acting on the strings of the columns
            void SameSpin(const MatrixXd &in, MatrixXd &out) const {
                const Index nStr = NumStrings();
                std::vector<VectorXd> F(NumComputeThreads(), VectorXd::Zero(nStr));
                std::vector<std::vector<uint32_t> > touched(NumComputeThreads());

                ParallelForDynamic((size_t) nStr, [&](size_t I, int thread) {
                    VectorXd &f = F[thread];
                    std::vector<uint32_t> &t = touched[thread];
                    auto Add = [&](uint32_t J, double value) {
                        if (f(J) == 0.0)
                            t.push_back(J);
                        f(J) += value;
                    };
                    for (const StringReplacement &r1: mSpace.replacements[I]) {
                        const Index p = r1.pq % mN, q = r1.pq / mN;
                        Add(r1.target, r1.sign * mK(p, q));
                        for (const StringReplacement &r2: mSpace.replacements[r1.target])
                            Add(r2.target, 0.5 * r1.sign * r2.sign * mG(r1.pq, r2.pq));
                    }
                    for (uint32_t J: t) {
                        if (f(J) != 0.0)
                            out.col((Index) I) += f(J) * in.col(J);
                        f(J) = 0.0;
                    }
                    t.clear();
                });
            }

            //out(Ib, Ia) += sum (pq|rs) <Ia|E_pq|Ja> <Ib|E_rs|Jb> C(Jb, Ja)
            void OppositeSpin(const Map<const MatrixXd> &C, MatrixXd &out) const {
                const Index nStr = NumStrings();
                ParallelForDynamic((size_t) nStr, [&](size_t Ia, int) {
                    for (const StringReplacement &a: mSpace.replacements[Ia]) {
                        const auto column = C.col(a.target);
                        const auto g = mG.row(a.pq);
                        for (Index Ib = 0; Ib < nStr; Ib++) {
                            double value = 0.0;
                            for (const StringReplacement &b: mSpace.replacements[Ib])
                                value += b.sign * g(b.pq) * column(b.target);
                            out(Ib, (Index) Ia) += a.sign * value;
                        }
                    }
                });
            }

            StringSpace mSpace;
            Index mN;
            MatrixXd mG, mK;
        };
    }
    //
    CI_Wfn CI(const HFWfn &rhfWfn,
              const Integrals::MolInts &mints,
              double Vnn,
              int excitorLevel,
              const DH::Thread *thread,
              std::stringstream *console) {
        if (console)
            *console << "--- Configuration Interaction --- \n\n    Transforming the ERI integrals...\n";

        //Spatial MO integrals; the spin-orbital elements are formed on demand
        const SpinOrbitalIntegrals ints = MakeSpinOrbitalIntegrals(rhfWfn, mints);
        if (console)
            *console << "   Transformation of integrals completed!\n";

        //Must be the same number of alpha and beta spin electrons for a closed-shell RHF determinant
        const uint32_t nElectrons = std::distance(rhfWfn.Occ_Alpha.begin(),
                                                  std::ranges::find(rhfWfn.Occ_Alpha, 0)) * 2;
        const uint32_t nSpatialOrbitals = rhfWfn.Occ_Alpha.size();
        const uint32_t nSpinOrbs = nSpatialOrbitals * 2;
        const bool fullCI = excitorLevel < 0;

        if (console)
            *console << "   Generating the Slater determinants...\n";
        auto genDetTime1 = std::chrono::high_resolution_clock::now();
        //Only the M_s = 0 determinants can couple to the closed-shell ground state
        const std::vector<SlaterDet> determinants =
                fullCI
                    ? EnumerateAllDet(nSpinOrbs, nElectrons / 2, nElectrons / 2)
                    : EnumerateDet(nSpinOrbs, nElectrons / 2, nElectrons / 2, static_cast<uint32_t>(excitorLevel));
        auto genDetTime2 = std::chrono::high_resolution_clock::now();
        const Index nDet = (Index) determinants.size();
        if (console)
            *console << "   Generated " << nDet << " determinants in " << std::chrono::duration_cast<
                std::chrono::milliseconds>(genDetTime2 - genDetTime1).count() << " milliseconds\n";

        if (thread && thread->ShouldStop())
            return {};

        auto Element = [&](SlaterDet I, SlaterDet J) {
            return CalculateHamiltonianElement(I, J, nSpatialOrbitals, nSpatialOrbitals, ints);
        };

        VectorXd energies;
        MatrixXd wfns;
        if ((size_t) nDet <= DENSE_CI_LIMIT) {
            if (console)
                *console << "   Forming and diagonalising the full CI Hamiltonian matrix\n";
            MatrixXd H = MatrixXd::Zero(nDet, nDet);
            ParallelForDynamic((size_t) nDet, [&](size_t i, int) {
                for (Index j = (Index) i; j < nDet; j++)
                    H(j, (Index) i) = Element(determinants[i], determinants[j]);
            });
            if (thread && thread->ShouldStop())
                return {};

            const SelfAdjointEigenSolver<MatrixXd> eigensolver(H.selfadjointView<Lower>());
            energies = eigensolver.eigenvalues();
            wfns = eigensolver.eigenvectors();
        } else {
            //The diagonal preconditions Davidson's method; the reference is the initial guess
            VectorXd diagonal(nDet);
            ParallelForDynamic((size_t) nDet, [&](size_t i, int) {
                diagonal((Index) i) = Element(determinants[i], determinants[i]);
            });
            MatrixXd guess = MatrixXd::Zero(nDet, 1);
            guess(0, 0) = 1.0;

            DavidsonResult davidson;
            if (fullCI) {
                if (console)
                    *console << "   Davidson diagonalisation with string-driven Hamiltonian products\n";
                const FCISigma sigma(ints, nSpatialOrbitals, nElectrons / 2);
                davidson = Davidson(sigma, diagonal, guess, 1e-7, 200, thread);
            } else {
                if (console)
                    *console << "   Davidson diagonalisation with a sparse Hamiltonian matrix\n";
                //Every row from the single and double excitations of its determinant present in the space
                const DeterminantIndex index(determinants);
                std::vector<std::vector<Triplet<double> > > rows(nDet);
                ParallelForDynamic((size_t) nDet, [&](size_t i, int) {
                    rows[i].emplace_back((Index) i, (Index) i, diagonal((Index) i));
                    for (const SlaterDet J: ConnectedDeterminants(determinants[i], nSpatialOrbitals)) {
                        const uint32_t j = index.Find(J);
                        if (j != index.Size())
                            rows[i].emplace_back((Index) i, (Index) j, Element(determinants[i], J));
                    }
                });
                std::vector<Triplet<double> > triplets;
                for (auto &row: rows) {
                    triplets.insert(triplets.end(), row.begin(), row.end());
                    std::vector<Triplet<double> >().swap(row);
                }
                SparseMatrix<double> H(nDet, nDet);
                H.setFromTriplets(triplets.begin(), triplets.end());
                std::vector<Triplet<double> >().swap(triplets);

                davidson = Davidson([&H](const MatrixXd &x) { return MatrixXd(H * x); },
                                    diagonal, guess, 1e-7, 200, thread);
            }
            if (thread && thread->ShouldStop())
                return {};
            if (!davidson.converged && console)
                *console << "!!! Davidson diagonalisation failed to converge !!!\n";
            energies = davidson.eigenvalues;
            wfns = davidson.eigenvectors;
        }

        for (double &energy: energies)
            energy += Vnn;

        if (console) {
            *console << "\n\n\n--- " << (fullCI ? "Full" : "Truncated") << " CI procedure completed! ---\n\n";
            *console << "   E(CI Correlation) = " << energies(0) - rhfWfn.E_HF << " hartrees\n";
            *console << "   E(Total) = E(Hartree Fock) + E(Correlation) = " << energies(0) << " hartrees\n";
            *console << "   Found " << wfns.cols() << " CI wavefunctions\n";
        }

        return CI_Wfn{energies, wfns};
    }
    //
    //
#pragma endregion
    //
    //
}
