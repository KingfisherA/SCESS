#pragma once
#include <cstdint>
#include <utility>
#include <vector>
#include <string>
#include <Eigen/Core>

#include "Dahlia/Utilities/Memory.h"
#include "Quantum/Integrals/Integrals.h"
#include "Quantum/HF/HF.h"


namespace DH {
    class Thread;
}

namespace SCESS::Integrals {
    struct MolInts;
}

namespace SCESS::HF {
    struct HF_Wfn;
}

namespace SCESS::CI {
    using HFWfn = HF::HF_Wfn;
    using TEITensor = Eigen::Tensor<double, 4>;
    //
    //
    using OrbBitStr = uint64_t;
    //
    class SlaterDet {
    public:
        //
        SlaterDet() {
        }
        //
        SlaterDet(const OrbBitStr alphaOrbs, const OrbBitStr betaOrbs)
            : mAlphaOrbs(alphaOrbs), mBetaOrbs(betaOrbs) {
        }
        //
        SlaterDet(const SlaterDet &o) = default;
        //
        //
        void Set(OrbBitStr alpha, OrbBitStr beta) {
            mAlphaOrbs = alpha;
            mBetaOrbs = beta;
        }
        //
        [[nodiscard]] OrbBitStr &GetAlphaOrbitals() {
            return mAlphaOrbs;
        }
        //
        [[nodiscard]] const OrbBitStr &GetAlphaOrbitals() const {
            return mAlphaOrbs;
        }
        //
        [[nodiscard]] OrbBitStr &GetBetaOrbitals() {
            return mBetaOrbs;
        }
        //
        [[nodiscard]] const OrbBitStr &GetBetaOrbitals() const {
            return mBetaOrbs;
        }
        //
        bool operator==(SlaterDet other) const {
            return mAlphaOrbs == other.mAlphaOrbs && mBetaOrbs == other.mBetaOrbs;
        }
        //
        SlaterDet operator^(SlaterDet other) const {
            return {mAlphaOrbs ^ other.mAlphaOrbs, mBetaOrbs ^ other.mBetaOrbs};
        }
        //
        static SlaterDet MakeReference(uint16_t nAlphaOcc, uint16_t nBetaOcc);
        //
        static void CalculateDifferences(SlaterDet I,
                                         SlaterDet J,
                                         uint16_t numOfAlphaOrbitals,
                                         uint16_t numOfBetaOrbitals,
                                         uint16_t *uniqueIArray,
                                         uint16_t *uniqueISize,
                                         uint16_t *uniqueJArray,
                                         uint16_t *uniqueJSize,
                                         uint16_t *commonIJArray,
                                         uint16_t *commonIJSize,
                                         int8_t *sign);
        //
        static int8_t CalculateSign(SlaterDet I, SlaterDet J,
                                    uint16_t numOfAlphaOrbitals,
                                    uint16_t numOfBetaOrbitals);
        //
        void CalculateSign(SlaterDet I,
                           uint16_t numOfAlphaOrbitals,
                           uint16_t numOfBetaOrbitals) {
            mSign = CalculateSign(I, *this, numOfAlphaOrbitals, numOfBetaOrbitals);
        }
        //
        [[nodiscard]] int8_t &GetSign() { return mSign; }
        [[nodiscard]] const int8_t &GetSign() const { return mSign; }
        //
    private:
        //Store the hatree products of orbitals in bitstrings
        //i.e. 1111100 would be the ground-state hartree fock determinant of water in sto-3g
        //Because of this we can only store a maximum of 64 orbitals in each bit string
        OrbBitStr mAlphaOrbs = 0, mBetaOrbs = 0;
        int8_t mSign = 1;
    };
    //
    //
    struct CI_Wfn {
    public:
        Eigen::VectorXd E_CI;
        Eigen::MatrixXd wfns;
    };
    //
    //
    /**
     * Perform Full CI on a RHF determinant
     * @param rhfWfn
     * @param excitorLevel Can be set to -1 for Full-CI
     * @return
     */
    CI_Wfn CI(const HFWfn &rhfWfn,
              const Integrals::MolInts &mints,
              double Vnn,
              int excitorLevel = -1,
              const DH::Thread *thread = nullptr,
              std::stringstream *console = nullptr);
    //
    /**
     * Generates the one electron hamiltonian for a RHF determinant with the spin
     * functions integrated out
     * @return
     */
    Eigen::MatrixXd GenerateHCoreSpinTransform(
        const Eigen::MatrixXd &h,
        const Eigen::MatrixXd &C);
    //
    Eigen::MatrixXd GenerateHCoreSpin(const Eigen::MatrixXd &h);
    //
    OrbBitStr MixAlphaAndBeta(OrbBitStr alphaList,
                              size_t numAlpha,
                              OrbBitStr betaList,
                              size_t numBeta);
    //
    std::vector<SlaterDet> EnumerateDet(uint32_t nSpinOrbs, uint32_t nAlphaElectrons,
                                        uint32_t nBetaElectrons, uint32_t excitorLevel);
    //Every determinant of nElectrons in nSpinOrbs spin orbitals (all M_s sectors)
    std::vector<SlaterDet> EnumerateAllDet(uint32_t nSpinOrbs, uint32_t nElectrons);
    //Every determinant with the given number of alpha and beta electrons, starting with the reference
    std::vector<SlaterDet> EnumerateAllDet(uint32_t nSpinOrbs, uint32_t nAlphaElectrons, uint32_t nBetaElectrons);
    //
    size_t GetNumOccOrbs(OrbBitStr orbitalBitStr, size_t numOfOrbitals);
    //
    //numOfOrbitals is the total size of occupied and virtual orbitals.
    //
    double CalculateHamiltonianElement(SlaterDet detI,
                                       SlaterDet detJ,
                                       size_t numOfAlphaOrbitals,
                                       size_t numOfBetaOrbitals,
                                       const Eigen::MatrixXd &HSpin,
                                       const TEITensor &spinERI);
    //
    void PrintDeterminants(const std::vector<SlaterDet> &determinants,
                           size_t numOrbitals);
    //
    //
}
