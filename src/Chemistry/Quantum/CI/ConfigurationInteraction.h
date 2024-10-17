#pragma once
#include <cstdint>
#include <utility>
#include <vector>
#include <string>

#include "GraphicsEngine/Utilities/Memory.h"
#include <Eigen/Core>
#include<Chemistry/Quantum/Integrals/Integrals.h>

#include "Chemistry/Quantum/HF/HartreeFock.h"

namespace SCESS {
    struct RHFResult;
}

namespace SCESS {
    class CIDeterminant {
    public:
        using OrbitalStr = std::basic_string<uint8_t>;

        CIDeterminant() {
        }

        CIDeterminant(const OrbitalStr &alphaOrbs, const OrbitalStr &betaOrbs)
            : mAlphaOrbs(alphaOrbs), mBetaOrbs(betaOrbs) {
        }

        CIDeterminant(OrbitalStr &&alphaOrbs, OrbitalStr &&betaOrbs)
            : mAlphaOrbs(std::move(alphaOrbs)), mBetaOrbs(std::move(betaOrbs)) {
        }

        CIDeterminant(const CIDeterminant &o): mAlphaOrbs(o.mAlphaOrbs), mBetaOrbs(o.mBetaOrbs) {
        }

        CIDeterminant(CIDeterminant &&o): mAlphaOrbs(std::move(o.mAlphaOrbs)), mBetaOrbs(std::move(o.mBetaOrbs)) {
        }

        void Set(const OrbitalStr &alpha, const OrbitalStr &beta) {
            mAlphaOrbs = alpha;
            mBetaOrbs = beta;
        }

        void Set(OrbitalStr &&alpha, OrbitalStr &&beta) {
            mAlphaOrbs = std::move(alpha);
            mBetaOrbs = std::move(beta);
        }

        [[nodiscard]] OrbitalStr &GetAlphaOrbitals() {
            return mAlphaOrbs;
        }

        [[nodiscard]] const OrbitalStr &GetAlphaOrbitals() const {
            return mAlphaOrbs;
        }

        [[nodiscard]] OrbitalStr &GetBetaOrbitals() {
            return mBetaOrbs;
        }

        [[nodiscard]] const OrbitalStr &GetBetaOrbitals() const {
            return mBetaOrbs;
        }

        static void CalculateDifferences(
            const CIDeterminant &I,
            const CIDeterminant &J,
            CIDeterminant::OrbitalStr &uniqueI,
            CIDeterminant::OrbitalStr &uniqueJ,
            CIDeterminant::OrbitalStr &commonIJ,
            int &sign);

    private:
        OrbitalStr mAlphaOrbs, mBetaOrbs;
    };

    struct CIResult {
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
        static uRef<CIResult> FCI(const RHFResult &rhfWavefunction, double Vnn,
                                  std::stringstream *console = nullptr);


        /**
         * Generates the one electron hamiltonian for a RHF determinant with the spin
         * functions integrated out
         * @return
         */
        static Eigen::MatrixXd GenerateOneElectronHamiltonianSpin(const RHFResult& rhfWfn);

        static CIDeterminant::OrbitalStr MixAlphaAndBeta(const CIDeterminant::OrbitalStr &alphaList,
                                                         const CIDeterminant::OrbitalStr &betaList);

        static std::vector<CIDeterminant> GenerateAllDeterminants(uint32_t nSpinOrbs, uint32_t nElectrons);

        static double CalculateHamiltonianElement(const CIDeterminant &detI, const CIDeterminant &detJ,
                                                  const Eigen::MatrixXd &HSpin, const TEITensor &moSpinTEI);

        static void PrintDeterminants(const std::vector<CIDeterminant> &determinants);
    };
} // SCESS
