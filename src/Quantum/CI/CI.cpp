#include <cppitertools/combinations.hpp>
#include <boost/range/irange.hpp>
#include <boost/math/special_functions/binomial.hpp>
//
#include "CI.h"
#include "Quantum/Integrals/IntegralTransforms.h"
//
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
            reference.GetAlphaOrbitals() |= 0x1 << i;

        for (i = 0; i < nBetaOcc; i++)
            reference.GetBetaOrbitals() |= 0x1 << i;

        return reference;
    }
    //
    static uint16_t CountSetBits(uint64_t bitStr, uint32_t nBits) {
        uint32_t nBitsSet = 0;
        for (uint32_t i = 0; i < nBits; i++)
            if ((0x1 << i) & bitStr)
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
        //Arrays for uniqueIArray, uniqueJArray and commonIJArray must be at least as larger as
        //max(numOfAlphaOrbitals, numOfBetaOrbitals) in size

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
                if ((0x1 << i) & J) iel2++;

                //Occupied in I
                if ((0x1 << i) & I) {
                    iel1++;
                    //Occupied in J
                    if ((0x1 << i) & J) {
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
                    if ((0x1 << i) & J) {
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
                if ((0x1 << i) & J) iel2++;

                //Occupied in I
                if ((0x1 << i) & I) {
                    iel1++;
                    //Not occupied in J
                    if (!((0x1 << i) & J)) {
                        iexcit1++;
                        flips += (shift - iel1 + iexcit1);
                    }
                }
                //Not occupied in I
                else {
                    //Occupied in J
                    if ((0x1 << i) & J) {
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
        if (numOfAlphaOrbitals != numOfBetaOrbitals)
            throw std::runtime_error(
                "Only RHF determinants with equal numbers of alpha and beta number of electrons are supported currently");

        // std::vector<uint32_t> uniqueI, uniqueJ, commonIJ;
        // uniqueI.reserve(numOfAlphaOrbitals);
        // uniqueJ.reserve(numOfAlphaOrbitals);
        // commonIJ.reserve(numOfAlphaOrbitals);
        auto *uniqueI = reinterpret_cast<uint16_t *>(
            alloca(std::max(numOfAlphaOrbitals, numOfBetaOrbitals) * sizeof(uint16_t)));
        auto *uniqueJ = reinterpret_cast<uint16_t *>(
            alloca(std::max(numOfAlphaOrbitals, numOfBetaOrbitals) * sizeof(uint16_t)));
        auto *commonIJ = reinterpret_cast<uint16_t *>(
            alloca(std::max(numOfAlphaOrbitals, numOfBetaOrbitals) * sizeof(uint16_t)));
        uint16_t uniqueISize, uniqueJSize, commonIJSize;
        int8_t sign;
        SlaterDet::CalculateDifferences(detI, detJ,
                                        numOfAlphaOrbitals, numOfBetaOrbitals,
                                        uniqueI, &uniqueISize,
                                        uniqueJ, &uniqueJSize,
                                        commonIJ, &commonIJSize,
                                        &sign);

        switch (uniqueISize) {
            case 0: {
                double element = 0.0;

                for (uint16_t m = 0; m < commonIJSize; m++)
                    element += HSpin(commonIJ[m], commonIJ[m]);

                for (uint16_t m = 0; m < commonIJSize - 1; m++)
                    for (uint16_t n = m + 1; n < commonIJSize; n++)
                        element += spinERI(commonIJ[m], commonIJ[n], commonIJ[m], commonIJ[n]);

                return element;
            }

            case 1: {
                auto m = uniqueI[0];
                auto p = uniqueJ[0];

                double element = HSpin(m, p);

                for (uint16_t k = 0; k < commonIJSize; k++) {
                    auto n = commonIJ[k];
                    element += spinERI(m, n, p, n);
                }

                return sign * element;
            }
            case 2: {
                double element = spinERI(uniqueI[0],
                                         uniqueI[1],
                                         uniqueJ[0],
                                         uniqueJ[1]);

                return sign * element;
            }

            default: {
                return 0.0;
            }
        }
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

        std::unordered_set<uint32_t> occupiedSpinOrbs;
        std::unordered_set<uint32_t> virtualSpinOrbs;
        for (uint32_t i = 0; i < nSpinOrbs; i++)
            virtualSpinOrbs.insert(i);
        for (uint32_t i = 0; i < nAlphaElectrons; i++) {
            virtualSpinOrbs.erase(2 * i);
            occupiedSpinOrbs.insert(2 * i);
        }
        for (uint32_t i = 0; i < nBetaElectrons; i++) {
            virtualSpinOrbs.erase(2 * i + 1);
            occupiedSpinOrbs.insert(2 * i + 1);
        }

        if (excitorLevel > occupiedSpinOrbs.size())
            excitorLevel = occupiedSpinOrbs.size();

        if (excitorLevel > virtualSpinOrbs.size())
            excitorLevel = virtualSpinOrbs.size();

        uint32_t ex = 1;
        size_t nDeterminants = 1;
        for (ex = 1; ex <= excitorLevel; ex++) {
            nDeterminants += static_cast<size_t>(
                boost::math::binomial_coefficient<double>(virtualSpinOrbs.size(), ex) *
                boost::math::binomial_coefficient<double>(occupiedSpinOrbs.size(), ex));
        }

        std::vector<SlaterDet> determinants;
        determinants.reserve(nDeterminants);

        determinants.emplace_back(SlaterDet::MakeReference(nAlphaElectrons, nBetaElectrons));

        for (ex = 1; ex <= excitorLevel; ex++)
            for (auto &virtualSpinOrbComb: iter::combinations(virtualSpinOrbs, ex))
                for (auto &occupiedSpinOrbComb: iter::combinations(occupiedSpinOrbs, ex)) {
                    SlaterDet &det = determinants.emplace_back();
                    det = SlaterDet::MakeReference(nAlphaElectrons, nBetaElectrons);

                    for (const auto oSpinOrbIdx: occupiedSpinOrbComb) {
                        if (oSpinOrbIdx % 2 == 0) {
                            // If even then alpha spin
                            det.GetAlphaOrbitals() ^= 1 << (oSpinOrbIdx / 2);
                        } else {
                            // If odd then beta spin
                            det.GetBetaOrbitals() ^= 1 << (oSpinOrbIdx / 2);
                        }
                    }
                    for (const auto vSpinOrbIdx: virtualSpinOrbComb) {
                        if (vSpinOrbIdx % 2 == 0) {
                            // If even then alpha spin
                            det.GetAlphaOrbitals() |= 1 << (vSpinOrbIdx / 2);
                        } else {
                            // If odd then beta spin
                            det.GetBetaOrbitals() |= 1 << (vSpinOrbIdx / 2);
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
                    det.GetAlphaOrbitals() |= 1 << (spinOrbIndx / 2);
                    continue;
                }
                // If odd then beta spin
                det.GetBetaOrbitals() |= 1 << (spinOrbIndx / 2);
            }
        }

        return determinants;
    }
    //
    size_t GetNumberOfOccupiedOrbitals(OrbBitStr orbitalBitStr,
                                       size_t numOfOrbitals) {
        size_t numOccupied = 0;

        for (size_t i = 0; i < numOfOrbitals; i++)
            if ((0x1 << i) & orbitalBitStr)
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
    CI_Wfn CI(const HFWfn &rhfWfn,
              const Integrals::MolInts &mints,
              double Vnn,
              int excitorLevel,
              const DH::Thread *thread,
              std::stringstream *console) {
        const MatrixXd &C = rhfWfn.C_Alpha;

        if (console)
            *console <<
                    "--- Configuration Interaction --- \n\n    Transforming and antisymmetrising ERI integrals...\n";

        // Transform tei tensor into MO basis in physicists' notation and
        // integrate out spin functions
        TEITensor moSpinTEI =
                Integrals::Transforms::AntisymmetriseIntegrateOutSpinClosedShell(
                    Integrals::Transforms::TransformTEI(mints.eri, C));

        if (console)
            *console << "   Transformation of integrals completed!\n";


        // Transform the one-electron Hamiltonian matrix into the MO basis set
        // and integrate out spin functions
        auto h = mints.T + mints.V;
        MatrixXd HSpin = GenerateHCoreSpinTransform(h, rhfWfn.C_Alpha);

        //Must be the same number of alpha and beta spin electrons for a closed-shell RHF determinant
        uint32_t nElectrons = std::distance(rhfWfn.Occ_Alpha.begin(),
                                            std::ranges::find(rhfWfn.Occ_Alpha, 0)) * 2;
        uint32_t nSpatialOrbitals = rhfWfn.Occ_Alpha.size();
        uint32_t nSpinOrbs = nSpatialOrbitals * 2;


        if (console)
            *console << "   Generating all Slater determinants in Hilbert space...\n";
        auto genDetTime1 = std::chrono::high_resolution_clock::now();
        std::vector<SlaterDet> determinants = EnumerateAllDet(
            nSpinOrbs, excitorLevel < 0 ? nElectrons : static_cast<uint32_t>(excitorLevel));
        auto genDetTime2 = std::chrono::high_resolution_clock::now();
        if (console)
            *console << "   Generated " << determinants.size() << " determinants in " << std::chrono::duration_cast<
                std::chrono::milliseconds>(genDetTime2 - genDetTime1).count() << " milliseconds\n";


        if (console) {
            *console << "   Forming Full CI Hamiltonian Matrix\n";
            *console << "   !!! No symmetry specified; will have to evaluate all " << determinants.size() <<
                    " matrix elements !!!\n";
        }
        MatrixXd H = MatrixXd::Zero((Index) determinants.size(),
                                    (Index) determinants.size());

        if (thread)
            if (thread->ShouldStop())
                return {};

        for (uint32_t i = 0; i < H.rows(); i++)
            for (uint32_t j = i; j < H.cols(); j++) {
                H(i, j) = H(j, i) =
                          CalculateHamiltonianElement(determinants[i], determinants[j],
                                                      nSpatialOrbitals,
                                                      nSpatialOrbitals, HSpin, moSpinTEI);
            }

        if (console)
            *console << "   Diagonalising FCI matrix...\n";

        if (thread)
            if (thread->ShouldStop())
                return {};

        //Diagonalise the Full CI matrix
        SelfAdjointEigenSolver<MatrixXd> eigensolver(H);
        VectorXd E_FCI = eigensolver.eigenvalues().real();

        for (double &energy: E_FCI)
            energy += Vnn;

        MatrixXd wfns = eigensolver.eigenvectors().real();

        if (console) {
            *console << "\n\n\n--- Full CI procedure completed! ---\n\n";
            *console << "   E(FCI Correlation) = " << E_FCI(0) - rhfWfn.E_HF << " hartrees\n";
            *console << "   E(Total) = E(Hartree Fock) + E(Correlation) = " << E_FCI(0) << " hartrees\n";
            *console << "   Generated " << wfns.rows() << " FCI wavefunctions\n";
        }

        return CI_Wfn{E_FCI, wfns};
    }
    //
    //
#pragma endregion
    //
    //
}
