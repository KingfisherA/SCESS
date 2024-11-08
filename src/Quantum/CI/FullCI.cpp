#include <iostream>

#include "ConfigurationInteraction.h"
#include  "Quantum/HF/HartreeFock.h"
#include  "Quantum/Integrals/Integrals.h"
#include "Utilities/Utilities.h"

#include <cppitertools/combinations.hpp>
#include <boost/range/irange.hpp>
#include <boost/math/special_functions/binomial.hpp>

namespace SCESS {
    CIDeterminant CIDeterminant::MakeReference(uint32_t nAlphaOcc, uint32_t nBetaOcc) {
        CIDeterminant reference;
        uint32_t i = 0;

        for (; i < nAlphaOcc; i++)
            reference.GetAlphaOrbitals() |= 0x1 << i;

        for (i = 0; i < nBetaOcc; i++)
            reference.GetBetaOrbitals() |= 0x1 << i;

        return reference;
    }

    static uint32_t CountSetBits(uint64_t bitStr, uint32_t nBits) {
        uint32_t nBitsSet = 0;
        for (uint32_t i = 0; i < nBits; i++)
            if ((0x1 << i) & bitStr)
                nBitsSet++;

        return nBitsSet;
    }

    void CIDeterminant::CalculateDifferences(CIDeterminant I,
                                             CIDeterminant J,
                                             size_t numOfAlphaOrbitals,
                                             size_t numOfBetaOrbitals,
                                             std::vector<uint32_t> &uniqueI,
                                             std::vector<uint32_t> &uniqueJ,
                                             std::vector<uint32_t> &commonIJ,
                                             int &sign) {
        uint32_t nElectrons = CountSetBits(I.GetAlphaOrbitals(), numOfAlphaOrbitals) +
                              CountSetBits(I.GetBetaOrbitals(), numOfBetaOrbitals);
        uint32_t excitationLvl = CountSetBits(I.GetAlphaOrbitals() ^ J.GetAlphaOrbitals(), numOfAlphaOrbitals) +
                                 CountSetBits(I.GetBetaOrbitals() ^ J.GetBetaOrbitals(), numOfBetaOrbitals);

        uint32_t shift = nElectrons - excitationLvl;

        uint32_t flips = 0;

        auto CalcDifferences = [&flips, &shift, &commonIJ, &uniqueI, &uniqueJ](OrbitalBitStr I,
                                                                               OrbitalBitStr J,
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

    uRef<CIWavefunction> ConfigurationInteraction::FCI(const RHFWavefunction &rhfWavefunction, double Vnn,
                                                       std::stringstream *console) {
        const Eigen::MatrixXd &C = rhfWavefunction.Coeffs;
        const ConfigIntegrals &ints = rhfWavefunction.Integrals;

        if (console)
            *console <<
                    "--- Full Configuration Interaction --- \n\n    Transforming and antisymmetrising ERI integrals...\n";

        // Transform tei tensor into MO basis in physicists' notation and
        // integrate out spin functions
        TEITensor moSpinTEI =
                Integrals::Transforms::AntisymmetriseIntegrateOutSpinClosedShell(
                    Integrals::Transforms::ToMolecularTEI(ints.tei, C));

        if (console)
            *console << "   Transformation of integrals completed!\n";


        // Transform the one-electron Hamiltonian matrix into the MO basis set
        // and integrate out spin functions
        Eigen::MatrixXd HSpin = GenerateOneElectronHamiltonianSpin(rhfWavefunction);

        //Must be the same number of alpha and beta spin electrons for a closed-shell RHF determinant
        uint32_t nElectrons = std::distance(rhfWavefunction.Occupation.begin(),
                                            std::ranges::find(rhfWavefunction.Occupation, 0)) * 2,
                nSpatialOrbitals = rhfWavefunction.Occupation.size(),
                nSpinOrbs = nSpatialOrbitals * 2;


        if (console)
            *console << "   Generating all Slater determinants in Hilbert space...\n";
        auto genDetTime1 = std::chrono::high_resolution_clock::now();
        std::vector<CIDeterminant> determinants = GenerateAllDeterminants(nSpinOrbs, nElectrons);
        auto genDetTime2 = std::chrono::high_resolution_clock::now();
        if (console)
            *console << "   Generated " << determinants.size() << " determinants in " << std::chrono::duration_cast<
                std::chrono::milliseconds>(genDetTime2 - genDetTime1).count() << " milliseconds\n";


        if (console) {
            *console << "   Forming Full CI Hamiltonian Matrix\n";
            *console << "   !!! No symmetry specified; will have to evaluate all " << determinants.size() <<
                    " matrix elements !!!\n";
        }
        Eigen::MatrixXd H_FCI = Eigen::MatrixXd::Zero((Eigen::Index) determinants.size(),
                                                      (Eigen::Index) determinants.size());

        for (uint32_t i = 0; i < H_FCI.rows(); i++)
            for (uint32_t j = i; j < H_FCI.cols(); j++) {
                H_FCI(i, j) = H_FCI(j, i) =
                              CalculateHamiltonianElement(determinants[i], determinants[j],
                                                          nSpatialOrbitals,
                                                          nSpatialOrbitals, HSpin, moSpinTEI);
            }

        if (console)
            *console << "   Diagonalising FCI matrix...\n";
        //Diagonalise the Full CI matrix
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(H_FCI);
        Eigen::VectorXd Energies = eigensolver.eigenvalues().real();

        for (double &energy: Energies)
            energy += Vnn;

        Eigen::MatrixXd Wavefunctions = eigensolver.eigenvectors().real();

        if (console) {
            *console << "\n\n\n--- Full CI procedure completed! ---\n\n";
            *console << "   E(FCI Correlation) = " << Energies(0) - rhfWavefunction.Energy << " hartrees\n";
            *console << "   E(Total) = E(Hartree Fock) + E(Correlation) = " << Energies(0) << " hartrees\n";
            *console << "   Generated " << Wavefunctions.rows() << " FCI wavefunctions\n";
        }

        return std::make_unique<CIWavefunction>(Energies, Wavefunctions, std::move(determinants));
    }

    Eigen::MatrixXd ConfigurationInteraction::GenerateOneElectronHamiltonianSpin(const RHFWavefunction &rhfWfn) {
        const Eigen::MatrixXd &C = rhfWfn.Coeffs;
        // Transform the one-electron Hamiltonian matrix into the MO basis set
        // and integrate out spin functions
        Eigen::MatrixXd HSpin = Eigen::MatrixXd::Zero(C.rows() * 2, C.cols() * 2);
        Eigen::MatrixXd H_Transformed = C.transpose() * rhfWfn.Hcore * C;

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

    double ConfigurationInteraction::CalculateHamiltonianElement(CIDeterminant detI, CIDeterminant detJ,
                                                                 size_t numOfAlphaOrbitals,
                                                                 size_t numOfBetaOrbitals,
                                                                 const Eigen::MatrixXd &HSpin,
                                                                 const TEITensor &moSpinTEI) {
        if (numOfAlphaOrbitals != numOfBetaOrbitals)
            throw std::runtime_error(
                "Only RHF determinants with equal numbers of alpha and beta number of electrons are supported currently");


        std::vector<uint32_t> uniqueI, uniqueJ, commonIJ;
        uniqueI.reserve(numOfAlphaOrbitals);
        uniqueJ.reserve(numOfAlphaOrbitals);
        commonIJ.reserve(numOfAlphaOrbitals);

        int sign = 0;
        CIDeterminant::CalculateDifferences(detI, detJ, numOfAlphaOrbitals, numOfBetaOrbitals,
                                            uniqueI, uniqueJ, commonIJ, sign);

        switch (uniqueI.size()) {
            case 0: {
                double element = 0.0;

                for (const auto m: commonIJ)
                    element += HSpin(m, m);

                for (uint32_t m = 0; m < commonIJ.size() - 1; m++)
                    for (uint32_t n = m + 1; n < commonIJ.size(); n++)
                        element += moSpinTEI(commonIJ[m], commonIJ[n], commonIJ[m], commonIJ[n]);

                return element;
            }

            case 1: {
                uint32_t m = uniqueI[0], p = uniqueJ[0];

                double element = HSpin(m, p);

                for (auto &n: commonIJ)
                    element += moSpinTEI(m, n, p, n);

                return sign * element;
            }
            case 2: {
                double element = moSpinTEI(uniqueI[0], uniqueI[1], uniqueJ[0],
                                           uniqueJ[1]);

                return sign * element;
            }

            default: {
                return 0.0;
            }
        }
    }

    CIDeterminant::OrbitalBitStr ConfigurationInteraction::MixAlphaAndBeta(const CIDeterminant::OrbitalBitStr alphaList,
                                                                           size_t numAlpha,
                                                                           const CIDeterminant::OrbitalBitStr betaList,
                                                                           size_t numBeta) {
        CIDeterminant::OrbitalBitStr combinedIndices = 0;
        size_t i = 0;

        for (i = 0; i < numAlpha; i++)
            combinedIndices |= (alphaList >> i) & 0x1;

        for (i = 0; i < numBeta; i++)
            combinedIndices |= (betaList >> i) & 0x1;

        return combinedIndices;
    }

    std::vector<CIDeterminant> ConfigurationInteraction::GenerateAllDeterminants(
        uint32_t nSpinOrbs, uint32_t nElectrons) {
        auto nDeterminants = static_cast<size_t>(
            boost::math::binomial_coefficient<double>(nSpinOrbs, nElectrons));
        std::vector<CIDeterminant> determinants;
        determinants.reserve(nDeterminants);

        for (auto &&spinOrbCombination: iter::combinations(
                 boost::irange<uint8_t>(0, nSpinOrbs), nElectrons)) {
            CIDeterminant &det = determinants.emplace_back();

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

    size_t ConfigurationInteraction::GetNumberOfOccupiedOrbitals(CIDeterminant::OrbitalBitStr orbitalBitStr,
                                                                 size_t numOfOrbitals) {
        size_t numOccupied = 0;

        for (size_t i = 0; i < numOfOrbitals; i++)
            if ((0x1 << i) & orbitalBitStr)
                numOccupied++;

        return numOccupied;
    }


    void ConfigurationInteraction::PrintDeterminants(const std::vector<CIDeterminant> &determinants,
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
} // SCESS
