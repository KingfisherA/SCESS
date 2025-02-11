#pragma once
#include <cstdint>
#include <utility>
#include <vector>
#include <string>

#include "GraphicsEngine/Utilities/Memory.h"
#include <Eigen/Core>
#include"Quantum/Integrals/Integrals.h"

#include  "Quantum/HF/HartreeFock.h"
#include "Viewer/Application/Application.h"

namespace GEng {
    class Thread;
}

namespace SCESS {
    struct RHFWavefunction;
}

namespace SCESS {
    class CIDeterminant {
    public:
        using OrbitalBitStr = uint64_t;

        CIDeterminant() {
        }

        CIDeterminant(const OrbitalBitStr alphaOrbs, const OrbitalBitStr betaOrbs)
            : mAlphaOrbs(alphaOrbs), mBetaOrbs(betaOrbs) {
        }

        // CIDeterminant(OrbitalBitStr &&alphaOrbs, OrbitalBitStr &&betaOrbs)
        //     : mAlphaOrbs(std::move(alphaOrbs)), mBetaOrbs(std::move(betaOrbs)) {
        // }

        CIDeterminant(const CIDeterminant &o) = default;

        void Set(OrbitalBitStr alpha, OrbitalBitStr beta) {
            mAlphaOrbs = alpha;
            mBetaOrbs = beta;
        }

        [[nodiscard]] OrbitalBitStr &GetAlphaOrbitals() {
            return mAlphaOrbs;
        }

        [[nodiscard]] const OrbitalBitStr &GetAlphaOrbitals() const {
            return mAlphaOrbs;
        }

        [[nodiscard]] OrbitalBitStr &GetBetaOrbitals() {
            return mBetaOrbs;
        }

        [[nodiscard]] const OrbitalBitStr &GetBetaOrbitals() const {
            return mBetaOrbs;
        }

        bool operator==(CIDeterminant other) const {
            return mAlphaOrbs == other.mAlphaOrbs && mBetaOrbs == other.mBetaOrbs;
        }

        CIDeterminant operator^(CIDeterminant other) const {
            return {mAlphaOrbs ^ other.mAlphaOrbs, mBetaOrbs ^ other.mBetaOrbs};
        }

        static CIDeterminant MakeReference(uint32_t nAlphaOcc, uint32_t nBetaOcc);

        static void CalculateDifferences(
            CIDeterminant I,
            CIDeterminant J,
            size_t numOfAlphaOrbitals,
            size_t numOfBetaOrbitals,
            std::vector<uint32_t> &uniqueI,
            std::vector<uint32_t> &uniqueJ,
            std::vector<uint32_t> &commonIJ,
            int8_t &sign);

        static int8_t CalculateSign(CIDeterminant I, CIDeterminant J, size_t numOfAlphaOrbitals,
                                    size_t numOfBetaOrbitals);

        void CalculateSign(CIDeterminant I, size_t numOfAlphaOrbitals,
                           size_t numOfBetaOrbitals) {
            mSign = CalculateSign(I, *this, numOfAlphaOrbitals, numOfBetaOrbitals);
        }

        [[nodiscard]] int8_t &GetSign() { return mSign; }
        [[nodiscard]] const int8_t &GetSign() const { return mSign; }

    private:
        //Store the hatree products of orbitals in bitstrings
        //i.e. 1111100 would be the ground-state hartree fock determinant of water in sto-3g
        //Because of this we can only store a maximum of 64 orbitals in each bit string
        OrbitalBitStr mAlphaOrbs = 0, mBetaOrbs = 0;
        int8_t mSign = 1;
    };

    struct CIWavefunction {
    public:
        Eigen::VectorXd Energies;
        Eigen::MatrixXd Wavefunctions;
        std::vector<CIDeterminant> Determinants;
    };

    class ConfigurationInteraction {
    public:
        /**
         * Perform Full CI on a RHF determinant
         * @param rhfWavefunction
         * @return
         */
        static uRef<CIWavefunction> FCI(const RHFWavefunction &rhfWavefunction, double Vnn,
                                        const GEng::Thread *thread = nullptr,
                                        std::stringstream *console = nullptr);


        /**
         * Generates the one electron hamiltonian for a RHF determinant with the spin
         * functions integrated out
         * @return
         */
        static Eigen::MatrixXd GenerateOneElectronHamiltonianSpin(
            const Eigen::MatrixXd& HCore,
            const Eigen::MatrixXd& C);

        static CIDeterminant::OrbitalBitStr MixAlphaAndBeta(const CIDeterminant::OrbitalBitStr alphaList,
                                                            size_t numAlpha,
                                                            const CIDeterminant::OrbitalBitStr betaList,
                                                            size_t numBeta);

        static std::vector<CIDeterminant> GenerateAllDeterminants(uint32_t nSpinOrbs, uint32_t nElectrons);

        static size_t GetNumberOfOccupiedOrbitals(CIDeterminant::OrbitalBitStr orbitalBitStr, size_t numOfOrbitals);

        //numOfOrbitals is the total size of occupied and virtual orbitals.
        //
        static double CalculateHamiltonianElement(CIDeterminant detI, CIDeterminant detJ, size_t numOfAlphaOrbitals,
                                                  size_t numOfBetaOrbitals,
                                                  const Eigen::MatrixXd &HSpin, const TEITensor &moSpinTEI);

        static void PrintDeterminants(const std::vector<CIDeterminant> &determinants,
                                      size_t numOrbitals);
    };
} // SCESS
