#include <bit>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/math/special_functions/sign.hpp>
//
#include "Quantum/QMC/QMC.h"
#include "Quantum/CI/CI.h"
#include "Dahlia/Tools/Threading.h"
#include "Quantum/Integrals/IntegralTransforms.h"
#include "Tensor/TensorParallel.h"

//
//
// Linked coupled-cluster; <Di |e^{-T} H e^{T}| D0>
// Unlinked coupled-cluster; <Di | H e^{T}| D0>
namespace SCESS::QMC::CC {
    //
    //
    using namespace Eigen;
    using namespace boost::random;
    //
    //
#pragma region CCQMC Utility Functions
    //
    //
    static OrbBitStr OrbitalMask(uint32_t nSpatialOrbitals) {
        return nSpatialOrbitals >= 64 ? ~OrbBitStr(0) : (OrbBitStr(1) << nSpatialOrbitals) - 1;
    }
    //
    //Index of the n-th occupied (or unoccupied) orbital of a bit string
    static uint32_t NthOrbital(OrbBitStr bitStr, uint32_t n, bool occupied, uint32_t nSpatialOrbitals) {
        for (uint32_t i = 0; i < nSpatialOrbitals; i++)
            if ((((bitStr >> i) & 1ULL) != 0) == occupied) {
                if (n == 0)
                    return i;
                --n;
            }

        return nSpatialOrbitals;
    }
    //
    uint32_t CountNumElecs(OrbBitStr bitStr, const uint32_t nSpatialOrbitals) {
        return std::popcount(bitStr & OrbitalMask(nSpatialOrbitals));
    }
    //
    SlaterDet RandomSingleExcite(SlaterDet I,
                                 mt19937 &randGenerator,
                                 const uint32_t nSpatialOrbitals,
                                 double &pExcit) {
        //Choose an electron uniformly from every occupied spin orbital and then an unoccupied orbital of the same
        //spin uniformly. Each single excitation is generated in exactly one way.
        const uint32_t nAlpha = CountNumElecs(I.GetAlphaOrbitals(), nSpatialOrbitals);
        const uint32_t nBeta = CountNumElecs(I.GetBetaOrbitals(), nSpatialOrbitals);
        const uint32_t nOccupied = nAlpha + nBeta;

        if (nOccupied == 0)
            return I;

        uniform_int_distribution<uint32_t> electronDist(0, nOccupied - 1);
        const uint32_t electron = electronDist(randGenerator);
        const bool isAlpha = electron < nAlpha;
        OrbBitStr &bitStr = isAlpha ? I.GetAlphaOrbitals() : I.GetBetaOrbitals();

        const uint32_t nVirtual = nSpatialOrbitals - CountNumElecs(bitStr, nSpatialOrbitals);
        if (nVirtual == 0)
            return I;

        uniform_int_distribution<uint32_t> virtualDist(0, nVirtual - 1);
        const uint32_t orbitalToDelete = NthOrbital(bitStr, isAlpha ? electron : electron - nAlpha, true,
                                                    nSpatialOrbitals);
        const uint32_t orbitalToCreate = NthOrbital(bitStr, virtualDist(randGenerator), false, nSpatialOrbitals);

        bitStr ^= (1ULL << orbitalToDelete) | (1ULL << orbitalToCreate);
        pExcit *= 1.0 / (double(nOccupied) * double(nVirtual));

        return I;
    }
    //
    SlaterDet RandomDoubleExcite(SlaterDet I,
                                 mt19937 &randGenerator,
                                 const uint32_t nSpatialOrbitals,
                                 double &pExcit) {
        //Choose an unordered pair of distinct electrons uniformly, then a pair of unoccupied orbitals with the same
        //spins uniformly. Each double excitation is generated in exactly one way.
        const uint32_t nAlpha = CountNumElecs(I.GetAlphaOrbitals(), nSpatialOrbitals);
        const uint32_t nBeta = CountNumElecs(I.GetBetaOrbitals(), nSpatialOrbitals);
        const uint32_t nOccupied = nAlpha + nBeta;

        //Cannot perform a double excitation if there are less than two electrons
        if (nOccupied < 2)
            return I;

        uniform_int_distribution<uint32_t> electron1Dist(0, nOccupied - 1);
        uniform_int_distribution<uint32_t> electron2Dist(0, nOccupied - 2);
        uint32_t electron1 = electron1Dist(randGenerator);
        uint32_t electron2 = electron2Dist(randGenerator);
        if (electron2 >= electron1)
            electron2++;
        if (electron2 < electron1)
            std::swap(electron1, electron2);
        const double pElectrons = 2.0 / (double(nOccupied) * double(nOccupied - 1));

        const uint32_t nVirtAlpha = nSpatialOrbitals - nAlpha;
        const uint32_t nVirtBeta = nSpatialOrbitals - nBeta;

        if (electron2 < nAlpha || electron1 >= nAlpha) {
            //Same spin pair
            const bool isAlpha = electron2 < nAlpha;
            const uint32_t nVirtual = isAlpha ? nVirtAlpha : nVirtBeta;
            if (nVirtual < 2)
                return I;

            OrbBitStr &bitStr = isAlpha ? I.GetAlphaOrbitals() : I.GetBetaOrbitals();
            const uint32_t offset = isAlpha ? 0 : nAlpha;

            uniform_int_distribution<uint32_t> virtual1Dist(0, nVirtual - 1);
            uniform_int_distribution<uint32_t> virtual2Dist(0, nVirtual - 2);
            uint32_t virtual1 = virtual1Dist(randGenerator);
            uint32_t virtual2 = virtual2Dist(randGenerator);
            if (virtual2 >= virtual1)
                virtual2++;

            const OrbBitStr original = bitStr;
            bitStr ^= (1ULL << NthOrbital(original, electron1 - offset, true, nSpatialOrbitals)) |
                      (1ULL << NthOrbital(original, electron2 - offset, true, nSpatialOrbitals)) |
                      (1ULL << NthOrbital(original, virtual1, false, nSpatialOrbitals)) |
                      (1ULL << NthOrbital(original, virtual2, false, nSpatialOrbitals));

            pExcit *= pElectrons * 2.0 / (double(nVirtual) * double(nVirtual - 1));
        } else {
            //Opposite spin pair; electron1 is alpha and electron2 is beta
            if (nVirtAlpha == 0 || nVirtBeta == 0)
                return I;

            uniform_int_distribution<uint32_t> virtualAlphaDist(0, nVirtAlpha - 1);
            uniform_int_distribution<uint32_t> virtualBetaDist(0, nVirtBeta - 1);

            OrbBitStr &alpha = I.GetAlphaOrbitals();
            OrbBitStr &beta = I.GetBetaOrbitals();
            alpha ^= (1ULL << NthOrbital(alpha, electron1, true, nSpatialOrbitals)) |
                     (1ULL << NthOrbital(alpha, virtualAlphaDist(randGenerator), false, nSpatialOrbitals));
            beta ^= (1ULL << NthOrbital(beta, electron2 - nAlpha, true, nSpatialOrbitals)) |
                    (1ULL << NthOrbital(beta, virtualBetaDist(randGenerator), false, nSpatialOrbitals));

            pExcit *= pElectrons / (double(nVirtAlpha) * double(nVirtBeta));
        }

        return I;
    }
    //
    bool ApplyExcitor(Excitor excitor,
                      SlaterDet reference,
                      SlaterDet &det,
                      int &sign,
                      uint32_t nSpatialOrbitals) {
        const OrbBitStr mask = OrbitalMask(nSpatialOrbitals);

        //Phase of a creation or annihilation operator acting on det: (-1)^(number of occupied spin orbitals before
        //it) with every alpha spin orbital ordered before every beta spin orbital
        auto Phase = [&det, mask](bool isAlpha, uint32_t orbital) {
            const OrbBitStr below = (1ULL << orbital) - 1;
            const int count = isAlpha
                                  ? std::popcount(det.GetAlphaOrbitals() & below)
                                  : std::popcount(det.GetAlphaOrbitals() & mask) +
                                    std::popcount(det.GetBetaOrbitals() & below);
            return count % 2 == 0 ? 1 : -1;
        };

        auto Annihilate = [&](bool isAlpha, OrbBitStr orbitals) {
            OrbBitStr &bitStr = isAlpha ? det.GetAlphaOrbitals() : det.GetBetaOrbitals();
            for (uint32_t i = 0; i < nSpatialOrbitals; i++)
                if ((orbitals >> i) & 1ULL) {
                    if (!((bitStr >> i) & 1ULL))
                        return false;
                    sign *= Phase(isAlpha, i);
                    bitStr &= ~(1ULL << i);
                }
            return true;
        };

        auto Create = [&](bool isAlpha, OrbBitStr orbitals) {
            OrbBitStr &bitStr = isAlpha ? det.GetAlphaOrbitals() : det.GetBetaOrbitals();
            for (uint32_t i = nSpatialOrbitals; i-- > 0;)
                if ((orbitals >> i) & 1ULL) {
                    if ((bitStr >> i) & 1ULL)
                        return false;
                    sign *= Phase(isAlpha, i);
                    bitStr |= 1ULL << i;
                }
            return true;
        };

        const OrbBitStr alpha = excitor.GetAlphaOrbitals() & mask;
        const OrbBitStr beta = excitor.GetBetaOrbitals() & mask;

        //Annihilators are applied in ascending spin orbital order then the creators in descending order
        return Annihilate(true, alpha & reference.GetAlphaOrbitals()) &&
               Annihilate(false, beta & reference.GetBetaOrbitals()) &&
               Create(false, beta & ~reference.GetBetaOrbitals()) &&
               Create(true, alpha & ~reference.GetAlphaOrbitals());
    }
    //
    uint32_t FindExcitorIdx(Excitor ex, const std::vector<Excitor> &excitors) {
        return std::distance(excitors.begin(), std::ranges::find(excitors, ex));
    }
    //
    static std::vector<Excitor> MakeExcitors(std::vector<SlaterDet> determinants,
                                             SlaterDet reference,
                                             uint32_t nSpatialOrbitals) {
        for (Excitor &ex: determinants) {
            ex.CalculateSign(reference, nSpatialOrbitals, nSpatialOrbitals);
            ex.GetAlphaOrbitals() ^= reference.GetAlphaOrbitals();
            ex.GetBetaOrbitals() ^= reference.GetBetaOrbitals();
        }

        return determinants;
    }
    //
    std::vector<Excitor> EnumerateExcitors(SlaterDet reference, uint32_t nSpinOrbs, uint32_t nAlphaElectrons,
                                           uint32_t nBetaElectrons, uint32_t excitationLevel) {
        return MakeExcitors(CI::EnumerateDet(nSpinOrbs, nAlphaElectrons, nBetaElectrons, excitationLevel),
                            reference, nSpinOrbs / 2);
    }
    //
    std::vector<Excitor> EnumerateAllExcitors(SlaterDet reference,
                                              uint32_t nSpinOrbs,
                                              uint32_t nAlphaElectrons,
                                              uint32_t nBetaElectrons) {
        return MakeExcitors(CI::EnumerateAllDet(nSpinOrbs, nAlphaElectrons, nBetaElectrons),
                            reference, nSpinOrbs / 2);
    }
    //
    int64_t CountExcips(const std::vector<int64_t> &excitorWalkers) {
        int64_t nExcips = 0;
        for (const int64_t walkers: excitorWalkers)
            nExcips += std::abs(walkers);

        return nExcips;
    }
    //
    std::vector<double> GenerateCumulativeClusterProb(uint32_t maxClusterSize) {
        //The cumulative 1/2^(s+1) probability distribution
        std::vector<double> cumulativeClusterProb(maxClusterSize + 1, 0.0);
        cumulativeClusterProb[0] = 0.5;
        for (uint32_t s = 1; s < cumulativeClusterProb.size(); s++)
            cumulativeClusterProb[s] = cumulativeClusterProb[s - 1] + 1.0 / pow(2.0, s + 1);

        //The largest cluster takes the remaining probability so that the distribution is normalised
        cumulativeClusterProb.back() = 1.0;

        return cumulativeClusterProb;
    }
    //
    uint32_t SampleClusterSize(double randomNumber,
                               const std::vector<double> &cumulativeClusterProb,
                               double &pSize) {
        uint32_t s = 0;
        while (s + 1 < cumulativeClusterProb.size() && randomNumber >= cumulativeClusterProb[s])
            s++;

        pSize = cumulativeClusterProb[s] - (s == 0 ? 0.0 : cumulativeClusterProb[s - 1]);
        return s;
    }
    //
    void EvaluateExcitorElements(const std::vector<Excitor> &excitors,
                                 const SlaterDet reference,
                                 const CI::SpinOrbitalIntegrals &ints,
                                 const uint32_t nSpatialOrbitals,
                                 VectorXd &diagonal,
                                 VectorXd &referenceRow) {
        diagonal.resize((Index) excitors.size());
        referenceRow.resize((Index) excitors.size());
        ParallelForDynamic(excitors.size(), [&](size_t i, int) {
            const SlaterDet det = excitors[i] ^ reference;
            diagonal((Index) i) = CI::CalculateHamiltonianElement(det, det, nSpatialOrbitals, nSpatialOrbitals, ints);
            referenceRow((Index) i) = CI::CalculateHamiltonianElement(reference, det, nSpatialOrbitals,
                                                                      nSpatialOrbitals, ints);
        });
    }
    //
    double CalculateProjectedEnergy(const std::vector<Excitor> &excitors,
                                    const std::vector<int64_t> &excitorWalkers,
                                    const VectorXd &referenceRow,
                                    const CI::DeterminantIndex &excitorIndex,
                                    SlaterDet reference,
                                    uint32_t nSpatialOrbitals) {
        const double N0 = (double) excitorWalkers[0];
        if (N0 == 0.0)
            return 0.0;

        const OrbBitStr mask = OrbitalMask(nSpatialOrbitals);
        auto ExcitationLevel = [mask](Excitor ex) {
            return (std::popcount(ex.GetAlphaOrbitals() & mask) + std::popcount(ex.GetBetaOrbitals() & mask)) / 2;
        };

        double E_Proj = 0.0;
        std::vector<uint32_t> singles;

        //Connected contributions, sign_j t_j <D0|H|Dj>
        for (uint32_t I = 1; I < excitors.size(); I++) {
            const auto level = ExcitationLevel(excitors[I]);
            if (level == 1)
                singles.push_back(I);
            if (level > 2 || excitorWalkers[I] == 0)
                continue;

            SlaterDet det = reference;
            int sign = 1;
            ApplyExcitor(excitors[I], reference, det, sign, nSpatialOrbitals);
            E_Proj += referenceRow(I) * sign * double(excitorWalkers[I]) / N0;
        }

        //Disconnected contributions of the T1^2 / 2 term onto the doubles
        for (size_t i = 0; i < singles.size(); i++)
            for (size_t k = i + 1; k < singles.size(); k++) {
                if (excitorWalkers[singles[i]] == 0 || excitorWalkers[singles[k]] == 0)
                    continue;

                SlaterDet det = reference;
                int sign = 1;
                if (!ApplyExcitor(excitors[singles[i]], reference, det, sign, nSpatialOrbitals) ||
                    !ApplyExcitor(excitors[singles[k]], reference, det, sign, nSpatialOrbitals))
                    continue;

                const uint32_t J = excitorIndex.Find(det ^ reference);
                if (J == excitors.size())
                    continue;

                E_Proj += referenceRow(J) * sign * (double(excitorWalkers[singles[i]]) / N0) *
                        (double(excitorWalkers[singles[k]]) / N0);
            }

        return E_Proj;
    }
    //
    void PrintCCQMCWavefunction(const std::vector<int64_t> &excitorWalkers) {
        double D0_pop = (double) excitorWalkers[0];

        std::cout << "--- CCQMC Wavefunction ---\n";
        std::cout << std::setprecision(15);
        for (uint32_t i = 0; i < excitorWalkers.size(); ++i) {
            double amplitude = double(excitorWalkers[i]) / D0_pop;
            std::cout << "  - c" << i << " = " << amplitude << '\n';
        }
        std::cout << "\n\n";
    }
    //
    //
#pragma endregion
    //
    //
#pragma region CCQMC Implementation
    //
    //
    double CCQMC(const HF_Wfn &rhfWfn,
                 const Integrals::MolInts &mints,
                 const DH::Thread *thread,
                 const ProjEnergyCallback &onReportCycle,
                 std::optional<uint64_t> seed) {
        constexpr size_t NUM_CYCLES = 100;
        constexpr size_t NUM_REPORT_CYCLES = 2000;
        //Number of report cycles after the shift starts varying before averaging the projected energy
        constexpr size_t EQUILIBRATION_CYCLES = 50;
        //The shift is held fixed until the population reaches this number of excips
        constexpr int64_t TARGET_POPULATION = 5000;
        //Damping of the shift update
        constexpr double SHIFT_DAMPING = 0.1;
        constexpr int64_t INITIAL_REFERENCE_POPULATION = 100;

        std::cout << "--- Coupled Cluster Quantum Monte Carlo ---\n";

        //Spatial MO integrals; the spin-orbital elements are formed on demand
        const CI::SpinOrbitalIntegrals ints = CI::MakeSpinOrbitalIntegrals(rhfWfn, mints);

        //The number of doubly occupied orbitals
        //Given we are starting with a RHF determinant, the number of alpha=beta which will equal this
        //quantity.
        uint32_t nOccupied = std::distance(rhfWfn.Occ_Alpha.begin(),
                                           std::ranges::find(rhfWfn.Occ_Alpha, 0));
        uint32_t nElectrons = nOccupied * 2;
        uint32_t nSpatialOrbitals = rhfWfn.Occ_Alpha.size();
        uint32_t nSpinOrbs = nSpatialOrbitals * 2;

        mt19937 randGenerator(seed ? (uint32_t) *seed : (uint32_t) time(nullptr));
        //Draw every random number from the same engine so that no two streams are correlated
        uniform_01<double> unitDist;
        auto randomNumber = [&randGenerator, &unitDist]() { return unitDist(randGenerator); };

        //The reference determinant
        const SlaterDet D0 = SlaterDet::MakeReference(nOccupied, nOccupied);
        //Generate all possible excitors in this full CC. The first excitor is the reference (the identity).
        std::vector<Excitor> excitors = EnumerateAllExcitors(D0, nSpinOrbs, nOccupied, nOccupied);
        const uint32_t nNonRefExcitors = excitors.size() - 1;
        //The populations of walkers on each of the excitors generated above
        std::vector<int64_t> excitorWalkers(excitors.size(), 0);
        excitorWalkers[0] = INITIAL_REFERENCE_POPULATION;
        std::vector<int64_t> excitorSpawn(excitors.size(), 0);

        //Index of every excitor, and the elements <Di|H|Di> and <D0|H|Di> which are needed for every excitor; the
        //other elements are evaluated as they are needed
        const CI::DeterminantIndex excitorIndex(excitors);
        VectorXd H_Diagonal, H_ReferenceRow;
        EvaluateExcitorElements(excitors, D0, ints, nSpatialOrbitals, H_Diagonal, H_ReferenceRow);
        if (thread && thread->ShouldStop())
            return 0.0;
        //Electronic energy of the reference; the shift is a correlation energy relative to it
        const double H00 = H_Diagonal(0);

        //The cumulative 1/2^(s+1) probability distribution, the largest cluster has every electron excited
        std::vector<double> cumulativeClusterProb = GenerateCumulativeClusterProb(nElectrons);

        double shift = 0.0;
        double tau = 0.01;
        bool varyShift = false;
        size_t shiftStartCycle = NUM_REPORT_CYCLES;

        double EProjAverage = 0.0;
        size_t nAveraged = 0;
        //Current total population of excips in the system
        int64_t nExcips = CountExcips(excitorWalkers);
        int64_t nExcipsLast = nExcips;

        for (size_t reportCycle = 0; reportCycle < NUM_REPORT_CYCLES; reportCycle++) {
            for (size_t cycle = 0; cycle < NUM_CYCLES; cycle++) {
                if (thread)
                    if (thread->ShouldStop())
                        return 0.0;

                for (auto &val: excitorSpawn)
                    val = 0;

                nExcips = CountExcips(excitorWalkers);
                const double N0 = (double) excitorWalkers[0];

                //Number of cluster selection attempts this cycle
                const int64_t Na = std::max<int64_t>(nExcips, 1);

                //Reused by every attempt so that the selection does not allocate
                std::vector<uint32_t> selectedExcitorIndices;
                for (int64_t attempt = 0; attempt < Na; attempt++) {
                    // -- Cluster selection --
                    // A cluster of size s is chosen with probability pSize and s distinct excitors are
                    // chosen uniformly. The 1/s! factor of the exp(T) expansion cancels with the s! orderings
                    // of the selection.
                    double pSize = 1.0;
                    const uint32_t s = SampleClusterSize(randomNumber(), cumulativeClusterProb, pSize);
                    if (s > nNonRefExcitors)
                        continue;

                    //Probability of selecting given excitors
                    double pClust = 1.0;
                    selectedExcitorIndices.clear();
                    if (s == 0) {
                        selectedExcitorIndices.push_back(0);
                    } else {
                        uniform_int_distribution<uint32_t> excitorIdxDist(1, nNonRefExcitors);
                        bool repeated = false;
                        for (uint32_t k = 0; k < s && !repeated; k++) {
                            const uint32_t newIdx = excitorIdxDist(randGenerator);
                            //Selecting the same excitor twice ends this attempt (rather than redrawing) so the
                            //probability of every valid cluster remains s! / n^s
                            repeated = std::ranges::find(selectedExcitorIndices, newIdx) !=
                                       selectedExcitorIndices.end();
                            selectedExcitorIndices.push_back(newIdx);
                        }
                        if (repeated)
                            continue;

                        pClust = Integrals::Utils::Factorial(s) * pow(1.0 / double(nNonRefExcitors), s);
                    }

                    //Collapse the excitors onto the reference, accounting for the fermionic sign. The cluster
                    //amplitude is t_e = prod_i N_i / N_0.
                    SlaterDet Di = D0;
                    int clusterSign = 1;
                    double clusterAmplitude = 1.0;
                    bool isValidCluster = true;
                    for (uint32_t idx: selectedExcitorIndices) {
                        if (idx != 0)
                            isValidCluster &= ApplyExcitor(excitors[idx], D0, Di, clusterSign, nSpatialOrbitals);
                        clusterAmplitude *= (double) excitorWalkers[idx] / N0;
                    }

                    if (!isValidCluster || clusterAmplitude == 0.0)
                        continue;

                    const uint32_t DiIdx = excitorIndex.Find(Di ^ D0);
                    if (DiIdx == excitors.size())
                        continue;

                    //Signed contribution of this cluster to the coefficient of |Di>, in excips
                    const double clusterWeight = clusterSign * N0 * clusterAmplitude;
                    const double pSelect = pSize * pClust;

                    // -- Spawning --
                    // A random single or double excitation of |Di> is performed to generate |Dj>. The excitor
                    // of |Dj> changes by dN_j = -tau sign_j H_ji weight, where sign_j is the sign of the excitor
                    // acting on the reference
                    double pExcit = 0.5;
                    SlaterDet Dj = randomNumber() < 0.5
                                       ? RandomSingleExcite(Di, randGenerator, nSpatialOrbitals, pExcit)
                                       : RandomDoubleExcite(Di, randGenerator, nSpatialOrbitals, pExcit);
                    const uint32_t DjIdx = excitorIndex.Find(Dj ^ D0);

                    //Cannot spawn onto itself and check J is valid
                    if (DiIdx != DjIdx && DjIdx != excitors.size()) {
                        const double Hji = CI::CalculateHamiltonianElement(Dj, Di, nSpatialOrbitals,
                                                                           nSpatialOrbitals, ints);

                        double pSpawn = tau * abs(clusterWeight * Hji) / (Na * pSelect * pExcit);
                        auto nSpawn = static_cast<int64_t>(pSpawn);
                        pSpawn = pSpawn - (double) nSpawn;
                        if (pSpawn > randomNumber())
                            nSpawn++;

                        const int spawnSign = -boost::math::sign(Hji) * boost::math::sign(clusterWeight) *
                                              excitors[DjIdx].GetSign();
                        excitorSpawn[DjIdx] += spawnSign * nSpawn;
                    }

                    // -- Death --
                    // dN_i = -tau sign_i (H_ii - H_00 - S) weight
                    const double Hii = H_Diagonal(DiIdx) - H00;

                    double pCloneDeath = tau * abs(clusterWeight) * (Hii - shift) / (Na * pSelect);
                    auto nKill = static_cast<int64_t>(pCloneDeath);
                    pCloneDeath = pCloneDeath - (double) nKill;
                    if (abs(pCloneDeath) > randomNumber()) {
                        if (pCloneDeath > 0.0)
                            nKill++;
                        else
                            nKill--;
                    }

                    excitorSpawn[DiIdx] -= nKill * boost::math::sign(clusterWeight) * excitors[DiIdx].GetSign();
                }

                // -- Annihilation --
                // Pairs of particles with opposite signs on the same excitor are removed
                for (uint32_t I = 0; I < excitorWalkers.size(); I++)
                    excitorWalkers[I] += excitorSpawn[I];

                if (excitorWalkers[0] == 0)
                    excitorWalkers[0] = INITIAL_REFERENCE_POPULATION;
            }

            nExcips = CountExcips(excitorWalkers);

            double EProj = rhfWfn.E_HF + CalculateProjectedEnergy(excitors, excitorWalkers, H_ReferenceRow,
                                                                  excitorIndex, D0, nSpatialOrbitals);

            if (onReportCycle)
                onReportCycle(reportCycle, EProj);

            if (!varyShift && nExcips >= TARGET_POPULATION) {
                varyShift = true;
                shiftStartCycle = reportCycle;
            }

            if (varyShift && reportCycle > shiftStartCycle && nExcipsLast != 0 && nExcips != 0)
                shift -= SHIFT_DAMPING / (tau * NUM_CYCLES) * log((double) nExcips / (double) nExcipsLast);

            if (varyShift && reportCycle >= shiftStartCycle + EQUILIBRATION_CYCLES) {
                EProjAverage += EProj;
                nAveraged++;
            }

            std::cout << "  Cycle(" << reportCycle << "); E(Proj) = " << EProj <<
                    " Shift = " << shift <<
                    " Population = " << nExcips;
            if (nAveraged != 0)
                std::cout << " E(Proj) Average = " << EProjAverage / double(nAveraged);
            std::cout << '\n';

            std::cout << "  Pop D0 = " << excitorWalkers[0] << "  Pop Dj = " << excitorWalkers.back() << '\n';

            nExcipsLast = nExcips;
        }

        if (nAveraged == 0) {
            std::cout << "!!! CCQMC population never reached the target; returning the final projected energy !!!\n";
            return rhfWfn.E_HF + CalculateProjectedEnergy(excitors, excitorWalkers, H_ReferenceRow, excitorIndex, D0,
                                                         nSpatialOrbitals);
        }

        return EProjAverage / double(nAveraged);
    }
    //
    //
#pragma endregion
    //
    //
}
