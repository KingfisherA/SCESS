#include <cppitertools/combinations.hpp>
#include <boost/range/irange.hpp>
#include <boost/math/special_functions/binomial.hpp>
//
#include "CI.h"
#include "Quantum/Integrals/IntegralTransforms.h"
//
#include "GraphicsEngine/Tools/Threading.h"


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
    SlaterDet SlaterDet::MakeReference(uint32_t nAlphaOcc, uint32_t nBetaOcc) {
        SlaterDet reference;
        uint32_t i = 0;

        for (; i < nAlphaOcc; i++)
            reference.GetAlphaOrbitals() |= 0x1 << i;

        for (i = 0; i < nBetaOcc; i++)
            reference.GetBetaOrbitals() |= 0x1 << i;

        return reference;
    }
    //
    static uint32_t CountSetBits(uint64_t bitStr, uint32_t nBits) {
        uint32_t nBitsSet = 0;
        for (uint32_t i = 0; i < nBits; i++)
            if ((0x1 << i) & bitStr)
                nBitsSet++;

        return nBitsSet;
    }
    //
    void SlaterDet::CalculateDifferences(SlaterDet I,
                                         SlaterDet J,
                                         size_t numOfAlphaOrbitals,
                                         size_t numOfBetaOrbitals,
                                         std::vector<uint32_t> &uniqueI,
                                         std::vector<uint32_t> &uniqueJ,
                                         std::vector<uint32_t> &commonIJ,
                                         int8_t &sign) {
        uint32_t nElectrons = CountSetBits(I.GetAlphaOrbitals(), numOfAlphaOrbitals) +
                              CountSetBits(I.GetBetaOrbitals(), numOfBetaOrbitals);
        uint32_t excitationLvl = CountSetBits(I.GetAlphaOrbitals() ^ J.GetAlphaOrbitals(), numOfAlphaOrbitals) +
                                 CountSetBits(I.GetBetaOrbitals() ^ J.GetBetaOrbitals(), numOfBetaOrbitals);

        uint32_t shift = nElectrons - excitationLvl;

        uint32_t flips = 0;

        auto CalcDifferences = [&flips, &shift, &commonIJ, &uniqueI, &uniqueJ](OrbBitStr I,
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
                        commonIJ.push_back(func(i));
                    }
                    //Not occupied in J
                    else {
                        iexcit1++;
                        flips += (shift - iel1 + iexcit1);
                        uniqueI.push_back(func(i));
                    }
                }
                //Not occupied in I
                else {
                    //Occupied in J
                    if ((0x1 << i) & J) {
                        iexcit2++;
                        uniqueJ.push_back(func(i));
                        flips += (shift - iel2 + iexcit2);
                    }
                }
            }
        };


        CalcDifferences(I.GetAlphaOrbitals(), J.GetAlphaOrbitals(), numOfAlphaOrbitals,
                        [](uint8_t i)-> uint8_t { return i * 2; });
        CalcDifferences(I.GetBetaOrbitals(), J.GetBetaOrbitals(), numOfBetaOrbitals,
                        [](uint8_t i)-> uint8_t { return i * 2 + 1; });

        sign = 1;
        if (flips % 2 == 1)
            sign = -1;
    }

    int8_t SlaterDet::CalculateSign(SlaterDet I, SlaterDet J, size_t numOfAlphaOrbitals,
                                    size_t numOfBetaOrbitals) {
        uint32_t nElectrons = CountSetBits(I.GetAlphaOrbitals(), numOfAlphaOrbitals) +
                              CountSetBits(I.GetBetaOrbitals(), numOfBetaOrbitals);
        uint32_t excitationLvl = CountSetBits(I.GetAlphaOrbitals() ^ J.GetAlphaOrbitals(), numOfAlphaOrbitals) +
                                 CountSetBits(I.GetBetaOrbitals() ^ J.GetBetaOrbitals(), numOfBetaOrbitals);

        uint32_t shift = nElectrons - excitationLvl;

        uint32_t flips = 0;

        auto CalcSign = [&flips, &shift](OrbBitStr I,
                                         OrbBitStr J,
                                         size_t numOfOrbs) {
            uint32_t iel1 = 0, iel2 = 0, iexcit1 = 0, iexcit2 = 0;

            for (uint32_t i = 0; i < numOfOrbs; i++) {
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

        std::vector<uint32_t> uniqueI, uniqueJ, commonIJ;
        uniqueI.reserve(numOfAlphaOrbitals);
        uniqueJ.reserve(numOfAlphaOrbitals);
        commonIJ.reserve(numOfAlphaOrbitals);

        int8_t sign = 1;
        SlaterDet::CalculateDifferences(detI, detJ, numOfAlphaOrbitals, numOfBetaOrbitals,
                                        uniqueI, uniqueJ, commonIJ, sign);

        switch (uniqueI.size()) {
            case 0: {
                double element = 0.0;

                for (const auto m: commonIJ)
                    element += HSpin(m, m);

                for (uint32_t m = 0; m < commonIJ.size() - 1; m++)
                    for (uint32_t n = m + 1; n < commonIJ.size(); n++)
                        element += spinERI(commonIJ[m], commonIJ[n], commonIJ[m], commonIJ[n]);

                return element;
            }

            case 1: {
                uint32_t m = uniqueI[0], p = uniqueJ[0];

                double element = HSpin(m, p);

                for (auto &n: commonIJ)
                    element += spinERI(m, n, p, n);

                return sign * element;
            }
            case 2: {
                double element = spinERI(uniqueI[0], uniqueI[1], uniqueJ[0],
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
    std::vector<SlaterDet> EnumerateDet(uint32_t nSpinOrbs, uint32_t excitorLevel) {
        auto nDeterminants = static_cast<size_t>(
            boost::math::binomial_coefficient<double>(nSpinOrbs, excitorLevel));
        std::vector<SlaterDet> determinants;
        determinants.reserve(nDeterminants);

        for (auto &&spinOrbCombination: iter::combinations(
                 boost::irange<uint8_t>(0, nSpinOrbs), excitorLevel)) {
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
              const GEng::Thread *thread,
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
        std::vector<SlaterDet> determinants = EnumerateDet(
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
