#include <iostream>

#include "ConfigurationInteraction.h"
#include "Chemistry/Quantum/HF/HartreeFock.h"
#include "Chemistry/Quantum/Integrals/Integrals.h"
#include "Utilities/Utilities.h"

#include <cppitertools/combinations.hpp>
#include <boost/range/irange.hpp>

namespace SCESS {
    void CIDeterminant::CalculateDifferences(const CIDeterminant &I, const CIDeterminant &J,
                                             CIDeterminant::OrbitalStr &uniqueI, CIDeterminant::OrbitalStr &uniqueJ,
                                             CIDeterminant::OrbitalStr &commonIJ,
                                             int &sign) {
        if (I.GetAlphaOrbitals().size() != J.GetAlphaOrbitals().size() ||
            I.GetBetaOrbitals().size() != J.GetBetaOrbitals().size())
            throw std::runtime_error("Alpha and beta electrons must be the same when calculating difference");

        uint32_t flips = 0;
        auto CalcDifferences = [&flips, &commonIJ, &uniqueI, &uniqueJ](const CIDeterminant::OrbitalStr &I,
                                                                       const CIDeterminant::OrbitalStr &J,
                                                                       uint8_t (*func)(uint8_t)) {
            int i = 0;
            int j = 0;
            int k = 0;
            int icnt = 0; /* number of counts in I array */
            int jcnt = 0; /* number of counts in J array */
            int flipI = 0; /* number of flips in I string to have max coincidence */
            int flipJ = 0; /* number of flips in J string to have max coincidence */
            int ndoI = 0; /* number of different orbitals in I string */
            int ndoJ = 0; /* number of different orbitals in J string */

            while ((i < I.size()) && (j < I.size())) {
                if (I[i] == J[j]) {
                    commonIJ.push_back(func(I[i]));
                    i++;
                    j++;
                    k++;
                } else if (I[i] < J[j]) {
                    uniqueI.push_back(func(I[i]));
                    flipI += i - ndoI;
                    ndoI++;
                    i++;
                } else if (I[i] > J[j]) {
                    uniqueJ.push_back(func(J[j]));
                    flipJ += j - ndoJ;
                    ndoJ++;
                    j++;
                }
            }

            if (i != j) {
                if (i < j) {
                    while (i < I.size()) {
                        uniqueI.push_back(func(I[i]));
                        flipI += i - ndoI;
                        ndoI++;
                        i++;
                    }
                } else {
                    while (j < I.size()) {
                        uniqueJ.push_back(func(J[j]));
                        flipJ += j - ndoJ;
                        ndoJ++;
                        j++;
                    }
                }
            }

            flips += flipI + flipJ;
        };

        CalcDifferences(I.GetAlphaOrbitals(), J.GetAlphaOrbitals(), [](uint8_t i)-> uint8_t { return i * 2; });
        CalcDifferences(I.GetBetaOrbitals(), J.GetBetaOrbitals(), [](uint8_t i)-> uint8_t { return i * 2 + 1; });

        sign = 1;
        if (flips % 2 == 0)
            sign = -1;
    }

    uRef<CIResult> ConfigurationInteraction::FCI(const RHFResult &rhfWavefunction, double Vnn) {
        const Eigen::MatrixXd &C = rhfWavefunction.Coeffs;
        const ConfigIntegrals &ints = rhfWavefunction.Integrals;

        // Transform tei tensor into MO basis in physicists' notation and
        // integrate out spin functions
        TEITensor moSpinTEI =
                Integrals::Transforms::AntisymmetriseIntegrateOutSpinClosedShell(
                    Integrals::Transforms::ToMolecularTEI(ints.tei, C));


        // Transform the one-electron Hamiltonian matrix into the MO basis set
        // and integrate out spin functions
        Eigen::MatrixXd HSpin = Eigen::MatrixXd::Zero(C.rows() * 2, C.cols() * 2); {
            Eigen::MatrixXd H_Transformed = C.transpose() * rhfWavefunction.Hcore * C;

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
        }

        //Must be the same number of alpha and beta spin electrons for a closed-shell RHF determinant
        uint32_t nElectrons = std::distance(rhfWavefunction.Occupation.begin(),
                                            std::ranges::find(rhfWavefunction.Occupation, 0)) * 2,
                nSpatialOrbitals = rhfWavefunction.Occupation.size(),
                nSpinOrbs = nSpatialOrbitals * 2;


        std::vector<CIDeterminant> determinants = GenerateAllDeterminants(nSpinOrbs, nElectrons);

        Eigen::MatrixXd H_FCI = Eigen::MatrixXd::Zero((Eigen::Index) determinants.size(),
                                                      (Eigen::Index) determinants.size());

        for (uint32_t i = 0; i < H_FCI.rows(); i++)
            for (uint32_t j = i; j < H_FCI.cols(); j++) {
                H_FCI(i, j) = H_FCI(j, i) =
                              CalculateHamiltonianElement(determinants[i], determinants[j], HSpin, moSpinTEI);
            }

        //Diagonalise the Full CI matrix
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(H_FCI);
        Eigen::VectorXd Energies = eigensolver.eigenvalues().real();

        for (double &energy: Energies)
            energy += Vnn;

        Eigen::MatrixXd Wavefunctions = eigensolver.eigenvectors().real();

        return std::make_unique<CIResult>(Energies, Wavefunctions, std::move(determinants));
    }

    double ConfigurationInteraction::CalculateHamiltonianElement(const CIDeterminant &detI, const CIDeterminant &detJ,
                                                                 const Eigen::MatrixXd &HSpin, const TEITensor &moTEI) {
        CIDeterminant::OrbitalStr uniqueI, uniqueJ, commonIJ;

        if (detI.GetAlphaOrbitals().size() != detJ.GetAlphaOrbitals().size() ||
            detI.GetBetaOrbitals().size() != detJ.GetBetaOrbitals().size())
            return 0.0;

        int sign = 0;
        CIDeterminant::CalculateDifferences(detI, detJ, uniqueI, uniqueJ, commonIJ, sign);

        switch (uniqueI.size()) {
            case 0: {
                double element = 0.0;

                for (const uint8_t m: commonIJ)
                    element += HSpin(m, m);

                for (uint32_t m = 0; m < commonIJ.size() - 1; m++)
                    for (uint32_t n = m + 1; n < commonIJ.size(); n++)
                        element += moTEI(commonIJ[m], commonIJ[n], commonIJ[m], commonIJ[n]);

                return element;
            }

            case 1: {
                uint32_t m = uniqueI[0], p = uniqueJ[0];

                double element = HSpin(m, p);

                for (auto &n: commonIJ)
                    element += moTEI(m, n, p, n);

                return sign * element;
            }
            case 2: {
                double element = moTEI(uniqueI[0], uniqueI[1], uniqueJ[0],
                                       uniqueJ[1]);

                return sign * element;
            }

            default: {
                return 0.0;
            }
        }
    }

    CIDeterminant::OrbitalStr ConfigurationInteraction::MixAlphaAndBeta(const CIDeterminant::OrbitalStr &alphaList,
                                                                        const CIDeterminant::OrbitalStr &betaList) {
        CIDeterminant::OrbitalStr combinedIndices(alphaList.size() + betaList.size(), (uint8_t) 0);

        size_t i = 0;

        for (const auto &m: alphaList)
            combinedIndices[i++] = m;

        for (const auto &m: betaList)
            combinedIndices[i++] = m;

        return combinedIndices;
    }

    std::vector<CIDeterminant> ConfigurationInteraction::GenerateAllDeterminants(
        uint32_t nSpinOrbs, uint32_t nElectrons) {
        auto nDeterminants = static_cast<size_t>(Integrals::Binomial(nSpinOrbs, nElectrons));
        nDeterminants = nDeterminants * nDeterminants;
        std::vector<CIDeterminant> determinants;
        determinants.reserve(nDeterminants);

        for (auto &&spinOrbCombination: iter::combinations(
                 boost::irange<uint8_t>(0, nSpinOrbs), nElectrons)) {
            CIDeterminant &det = determinants.emplace_back();

            CIDeterminant::OrbitalStr &alphaOrbs = det.GetAlphaOrbitals();
            CIDeterminant::OrbitalStr &betaOrbs = det.GetBetaOrbitals();
            alphaOrbs.reserve(nSpinOrbs);
            betaOrbs.reserve(nSpinOrbs);

            for (const auto &spinOrbIndx: spinOrbCombination) {
                // If even then alpha spin
                if (spinOrbIndx % 2 == 0) {
                    alphaOrbs.push_back(spinOrbIndx / 2);
                    continue;
                }
                // If odd then beta spin
                betaOrbs.push_back(spinOrbIndx / 2);
            }
        }

        return determinants;
    }


    void ConfigurationInteraction::PrintDeterminants(const std::vector<CIDeterminant> &determinants) {
        for (auto &det: determinants) {
            std::cout << "[Alpha, ";

            for (uint32_t i = 0; i < det.GetAlphaOrbitals().size(); i++)
                std::cout << (int) det.GetAlphaOrbitals()[i] << ", ";

            std::cout << "Beta, ";

            for (uint32_t i = 0; i < det.GetBetaOrbitals().size(); i++)
                std::cout << (int) det.GetBetaOrbitals()[i] << ", ";


            std::cout << "] \n";
        }
    }
} // SCESS
