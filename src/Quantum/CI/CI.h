#pragma once
#include <array>
#include <cstdint>
#include <stdexcept>
#include <functional>
#include <unordered_map>
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
    //Hash of the alpha and beta orbital strings (the sign is not part of the determinant)
    struct SlaterDetHash {
        size_t operator()(const SlaterDet &det) const {
            //splitmix64 finaliser of the two strings
            uint64_t x = det.GetAlphaOrbitals() * 0x9E3779B97F4A7C15ULL ^ det.GetBetaOrbitals();
            x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ULL;
            x = (x ^ (x >> 27)) * 0x94D049BB133111EBULL;
            return (size_t) (x ^ (x >> 31));
        }
    };
    //
    //Index of every determinant in a list, replacing linear searches of the list
    class DeterminantIndex {
    public:
        DeterminantIndex() = default;

        explicit DeterminantIndex(const std::vector<SlaterDet> &determinants) {
            mIndices.reserve(determinants.size());
            for (uint32_t i = 0; i < determinants.size(); i++)
                mIndices.emplace(determinants[i], i);
            mSize = (uint32_t) determinants.size();
        }

        //The index of det in the list, or the size of the list if it is not there
        [[nodiscard]] uint32_t Find(SlaterDet det) const {
            const auto it = mIndices.find(det);
            return it == mIndices.end() ? mSize : it->second;
        }

        [[nodiscard]] uint32_t Size() const { return mSize; }

    private:
        std::unordered_map<SlaterDet, uint32_t, SlaterDetHash> mIndices;
        uint32_t mSize = 0;
    };
    //
    //The spin-orbital Hamiltonian over spin orbitals p = 2 i + sigma (alpha even, beta odd) of the spatial molecular
    //orbitals i, evaluated from the spatial integrals on demand so that the (2n)^4 spin-orbital tensor of
    //Integrals::Transforms::AntisymmetriseIntegrateOutSpinClosedShell is never formed
    struct SpinOrbitalIntegrals {
        //One-electron Hamiltonian over the spatial orbitals
        Eigen::MatrixXd h;
        //(pq|rs) over the spatial orbitals, chemists' notation
        TEITensor eri;

        //<pq||rs> = <pq|rs> - <pq|sr> with the spin integrated out
        [[nodiscard]] double operator()(uint32_t p, uint32_t q, uint32_t r, uint32_t s) const {
            double value = 0.0;
            if ((p & 1u) == (r & 1u) && (q & 1u) == (s & 1u))
                value += eri(p >> 1, r >> 1, q >> 1, s >> 1);
            if ((p & 1u) == (s & 1u) && (q & 1u) == (r & 1u))
                value -= eri(p >> 1, s >> 1, q >> 1, r >> 1);
            return value;
        }

        //<p|h|q>
        [[nodiscard]] double H(uint32_t p, uint32_t q) const {
            return (p & 1u) == (q & 1u) ? h(p >> 1, q >> 1) : 0.0;
        }
    };
    //
    //The integrals over the canonical RHF orbitals
    SpinOrbitalIntegrals MakeSpinOrbitalIntegrals(const HFWfn &rhfWfn, const Integrals::MolInts &mints);
    //
    //Precomputed spin-orbital integrals: HSpin from GenerateHCoreSpinTransform and spinERI from
    //AntisymmetriseIntegrateOutSpinClosedShell
    struct DenseSpinOrbitalIntegrals {
        const Eigen::MatrixXd &HSpin;
        const TEITensor &spinERI;

        [[nodiscard]] double operator()(uint32_t p, uint32_t q, uint32_t r, uint32_t s) const {
            return spinERI(p, q, r, s);
        }

        [[nodiscard]] double H(uint32_t p, uint32_t q) const { return HSpin(p, q); }
    };
    //
    //
    struct CI_Wfn {
    public:
        //Total energies (with the nuclear repulsion) of the roots found, lowest first
        Eigen::VectorXd E_CI;
        //One column of determinant coefficients per root
        Eigen::MatrixXd wfns;
    };
    //
    //
    /**
     * Perform Full CI on a RHF determinant. Spaces of up to DENSE_CI_LIMIT determinants are diagonalised in full
     * (every root). Larger spaces find the lowest root with Davidson's method: string-driven Hamiltonian products
     * for full CI, or a sparse Hamiltonian for truncated CI.
     * @param rhfWfn
     * @param excitorLevel Can be set to -1 for Full-CI
     * @return
     */
    constexpr size_t DENSE_CI_LIMIT = 2000;
    //
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
    //<detI|H|detJ> from the Slater-Condon rules (without the nuclear repulsion). Ints provides the spin-orbital
    //integrals <pq||rs> through operator() and <p|h|q> through H(): SpinOrbitalIntegrals or DenseSpinOrbitalIntegrals.
    template<class Ints>
    double CalculateHamiltonianElement(SlaterDet detI,
                                       SlaterDet detJ,
                                       size_t numOfAlphaOrbitals,
                                       size_t numOfBetaOrbitals,
                                       const Ints &ints) {
        if (numOfAlphaOrbitals != numOfBetaOrbitals)
            throw std::runtime_error(
                "Only RHF determinants with equal numbers of alpha and beta number of electrons are supported currently");

        //The arrays hold spin orbitals of both spins (at most 64 alpha and 64 beta orbitals)
        std::array<uint16_t, 128> uniqueI{}, uniqueJ{}, commonIJ{};
        uint16_t uniqueISize, uniqueJSize, commonIJSize;
        int8_t sign;
        SlaterDet::CalculateDifferences(detI, detJ,
                                        (uint16_t) numOfAlphaOrbitals, (uint16_t) numOfBetaOrbitals,
                                        uniqueI.data(), &uniqueISize,
                                        uniqueJ.data(), &uniqueJSize,
                                        commonIJ.data(), &commonIJSize,
                                        &sign);

        switch (uniqueISize) {
            case 0: {
                double element = 0.0;
                for (uint16_t m = 0; m < commonIJSize; m++)
                    element += ints.H(commonIJ[m], commonIJ[m]);
                for (uint16_t m = 0; m + 1 < commonIJSize; m++)
                    for (uint16_t n = m + 1; n < commonIJSize; n++)
                        element += ints(commonIJ[m], commonIJ[n], commonIJ[m], commonIJ[n]);
                return element;
            }
            case 1: {
                const auto m = uniqueI[0];
                const auto p = uniqueJ[0];
                double element = ints.H(m, p);
                for (uint16_t k = 0; k < commonIJSize; k++) {
                    const auto n = commonIJ[k];
                    element += ints(m, n, p, n);
                }
                return sign * element;
            }
            case 2:
                return sign * ints(uniqueI[0], uniqueI[1], uniqueJ[0], uniqueJ[1]);
            default:
                return 0.0;
        }
    }
    //
    double CalculateHamiltonianElement(SlaterDet detI,
                                       SlaterDet detJ,
                                       size_t numOfAlphaOrbitals,
                                       size_t numOfBetaOrbitals,
                                       const Eigen::MatrixXd &HSpin,
                                       const TEITensor &spinERI);
    //
    //Every determinant which differs from det by a single or double excitation conserving the numbers of alpha and
    //beta electrons, over nSpatialOrbitals orbitals of each spin
    std::vector<SlaterDet> ConnectedDeterminants(SlaterDet det, uint32_t nSpatialOrbitals);
    //
    struct DavidsonResult {
        Eigen::VectorXd eigenvalues;
        Eigen::MatrixXd eigenvectors;
        bool converged = false;
    };
    //
    /**
     * Lowest eigenpairs of a real symmetric matrix from its products with vectors (Davidson's method with the
     * diagonal preconditioner)
     * @param sigma Computes H x for a block of column vectors x
     * @param diagonal The diagonal of H
     * @param guess Initial vectors, one column per root
     * @param residualTolerance Convergence on the norm of every residual H x - e x
     */
    DavidsonResult Davidson(const std::function<Eigen::MatrixXd(const Eigen::MatrixXd &)> &sigma,
                            const Eigen::VectorXd &diagonal,
                            const Eigen::MatrixXd &guess,
                            double residualTolerance = 1e-7,
                            size_t maxIterations = 200,
                            const DH::Thread *thread = nullptr);
    //
    void PrintDeterminants(const std::vector<SlaterDet> &determinants,
                           size_t numOrbitals);
    //
    //
}
