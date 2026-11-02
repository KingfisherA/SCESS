#include <Quantum/Q_PCH.h>
//
#include <bit>
#include <fmt/format.h>
//
#include <boost/math/special_functions/sign.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/random/uniform_int.hpp>
//
#include "F12_QMC.h"
#include "Quantum/F12/F12.h"
#include "Quantum/QMC/QMC.h"
#include "Quantum/CI/CI.h"
#include "Quantum/Integrals/IntegralsF12.h"
#include "Quantum/Integrals/IntegralTransforms.h"
//
#include "Tensor/TensorOperations.h"
#include "Tensor/TensorSlice.h"
//
#include <fmt/color.h>

#include "Dahlia/Tools/Threading.h"
//


//
//
namespace SCESS::F12_QMC {
    //
    //
    using namespace Eigen;
    using namespace boost::random;
    //
    //
#pragma region CCSD_QMC_MP2_F12 Utilities
    //
    //
    double EvaluateGeminalElement(SlaterDet Dj,
                                  ElectronPair gemPair,
                                  const Tensor4 &V_pqxy,
                                  const Tensor4 &h_pqxy,
                                  uint32_t nAlphaElectrons,
                                  uint32_t nBetaElectrons,
                                  uint32_t numOfAlphaOrbitals,
                                  uint32_t numOfBetaOrbitals) {
        const SlaterDet D0 = SlaterDet::MakeReference(nAlphaElectrons, nBetaElectrons);

        //Geminal electron 1 (alpha) and 2 (beta)
        auto [g1, g2] = gemPair;

        const OrbBitStr alphaRemoved = D0.GetAlphaOrbitals() & ~Dj.GetAlphaOrbitals();
        const OrbBitStr alphaAdded = Dj.GetAlphaOrbitals() & ~D0.GetAlphaOrbitals();
        const OrbBitStr betaRemoved = D0.GetBetaOrbitals() & ~Dj.GetBetaOrbitals();
        const OrbBitStr betaAdded = Dj.GetBetaOrbitals() & ~D0.GetBetaOrbitals();

        //Dj = Phi_{g1 g2}^{pq}; it may only differ from the reference by replacing the alpha orbital g1 and/or the
        //beta orbital g2
        if (std::popcount(alphaRemoved) > 1 || std::popcount(betaRemoved) > 1)
            return 0.0;
        if ((alphaRemoved != 0 && alphaRemoved != (1ULL << g1)) ||
            (betaRemoved != 0 && betaRemoved != (1ULL << g2)))
            return 0.0;

        const int p = alphaRemoved != 0 ? std::countr_zero(alphaAdded) : g1;
        const int q = betaRemoved != 0 ? std::countr_zero(betaAdded) : g2;

        //Phase of Dj relative to the excitor a+_p a+_q a_g2 a_g1 acting on the reference
        const int8_t sign = SlaterDet::CalculateSign(D0, Dj, numOfAlphaOrbitals, numOfBetaOrbitals);

        //<pq| H Q12 S f12 |g1 g2> with the fixed amplitude (SP) geminal S f12|g1 g2> = 3/8 f12|g1 g2> + 1/8 f12|g2 g1>
        return sign * ((3.0 / 8.0) * (V_pqxy(p, q, g1, g2) + h_pqxy(p, q, g1, g2))
                       + (1.0 / 8.0) * (V_pqxy(p, q, g2, g1) + h_pqxy(p, q, g2, g1)));
    }
    //
    ElectronPair GetGeminalElectronPair(Excitor gemExcitor,
                                        uint32_t nAlphaOrbitals,
                                        uint32_t nBetaOrbitals) {
        ElectronPair gemPair = {-1, -1};

        for (int i = 0; i < nAlphaOrbitals; i++)
            if ((1ULL << i) & gemExcitor.GetAlphaOrbitals()) {
                if (gemPair.first == -1)
                    gemPair.first = i;
                else {
                    gemPair.second = i;
                    return gemPair;
                }
            }

        for (int i = 0; i < nBetaOrbitals; i++)
            if ((1ULL << i) & gemExcitor.GetBetaOrbitals()) {
                if (gemPair.first == -1)
                    gemPair.first = i;
                else {
                    gemPair.second = i;
                    return gemPair;
                }
            }

        return gemPair;
    }
    //
    //
#pragma endregion
    //
    //
#pragma region CCSD_QMC_MP2_F12 Implementation
    //
    //
    double CCSD_QMC_MP2_F12(const HF::HF_Wfn &refWfn,
                            const Integrals::MolInts &mints,
                            const Integrals::F12::F12Integrals &f12Ints,
                            F12::F12Idxs f12Idxs) {
        using namespace QMC::CC;
        enum {
            //All OBS orbitals
            P, Q, R, S
        };
        TensorSliceIndices idxs = {
            //All OBS indices
            {P, {0, f12Idxs.nOBS}},
            {Q, {0, f12Idxs.nOBS}},
            {R, {0, f12Idxs.nOBS}},
            {S, {0, f12Idxs.nOBS}}
        };

        fmt::print("\n\n--- CCSD Quantum Monte Carlo MP2 F12; CCSD[MP2-F12] ---\n\n");


        // -- Parameters for the CCSD_QMC_MP2_F12 procedure --
        //Number of cycles
        constexpr size_t NUM_CYCLES = 100;
        //Number of report cycles
        constexpr size_t NUM_REPORT_CYCLES = 400;
        //Report cycle to start averaging the projected energy (once the shift is varying)
        constexpr size_t START_AVERAGE = 200;
        //Target excip population at which the shift starts to vary
        constexpr int64_t TARGET_POPULATION = 5000;
        //Damping of the shift update
        constexpr double SHIFT_DAMPING = 0.1;
        constexpr int64_t INITIAL_REFERENCE_POPULATION = 100;
        //CCSD truncation level of the excitors
        constexpr uint32_t TRUNCATION_LEVEL = 2;


        //The integrals from TransformF12Basis are in physicists' notation but the antisymmetrisation expects
        //chemists' notation
        TEITensor spinERI = Integrals::Transforms::AntisymmetriseIntegrateOutSpinClosedShell(
            Integrals::Transforms::ConvertPhysicistsToChemists(
                Slice<Indices<P, Q, R, S> >(idxs, mints.eri)));

        MatrixXd h = mints.T.block(0, 0, f12Idxs.nOBS, f12Idxs.nOBS) +
                     mints.V.block(0, 0, f12Idxs.nOBS, f12Idxs.nOBS);
        MatrixXd HSpin = CI::GenerateHCoreSpin(h);

        fmt::print("Evaluating F12 T2'(1) Intermediates...\n");
        Tensor2 F_pq;
        TEITensor V_ijxy, C_xyab, X_xywz, B_xywz;
        F12::MP2_F12_Fixed(refWfn, mints, f12Ints, f12Idxs, &F_pq, &V_ijxy, &C_xyab, &X_xywz, &B_xywz);
        TEITensor V_pqxy, h_pqxy;
        F12::EvaluateAdditionalMP2Intermediates(mints, f12Ints, f12Idxs,
                                                &V_pqxy, &h_pqxy);
        fmt::print(fg(fmt::color::green), "Done!\n");

        //The number of doubly occupied orbitals
        //Given we are starting with a RHF determinant, the number of alpha=beta which will equal this
        //quantity.
        uint32_t nOccupied = std::distance(refWfn.Occ_Alpha.begin(),
                                           std::ranges::find(refWfn.Occ_Alpha, 0));
        uint32_t nAlphaElectrons = nOccupied;
        uint32_t nBetaElectrons = nOccupied;
        uint32_t nElectrons = nAlphaElectrons + nBetaElectrons;
        uint32_t nSpatialOrbitals = refWfn.Occ_Alpha.size();
        uint32_t nSpinOrbs = nSpatialOrbitals * 2;

        mt19937 randGenerator(time(nullptr));
        //Draw every random number from the same engine so that no two streams are correlated
        uniform_01<double> unitDist;
        auto randomNumber = [&randGenerator, &unitDist]() { return unitDist(randGenerator); };

        //The reference determinant
        const SlaterDet D0 = SlaterDet::MakeReference(nOccupied, nOccupied);
        //Generate all single and double excitors for CCSD. The first excitor is the reference (the identity).
        std::vector<Excitor> excitors = EnumerateExcitors(D0, nSpinOrbs, nOccupied, nOccupied, TRUNCATION_LEVEL);
        //Number of conventional excitors including the reference
        const uint32_t nConventional = excitors.size();

        //Make the Hamiltonian matrix
        MatrixXd H_CI = MatrixXd::Zero(excitors.size(), excitors.size());
        for (uint32_t i = 0; i < H_CI.rows(); i++)
            for (uint32_t j = i; j < H_CI.cols(); j++)
                H_CI(i, j) = H_CI(j, i) = CI::CalculateHamiltonianElement(
                                 excitors[i] ^ D0,
                                 excitors[j] ^ D0,
                                 nSpatialOrbitals,
                                 nSpatialOrbitals,
                                 HSpin,
                                 spinERI);
        //Electronic energy of the reference; the shift is a correlation energy relative to it
        const double H00 = H_CI(0, 0);

        //Add the geminal (MEG) excitors, one for each pair of an alpha and a beta occupied orbital
        uint32_t nGeminalExcitors = 0;
        for (uint32_t i = 0; i < nAlphaElectrons; i++)
            for (uint32_t j = 0; j < nBetaElectrons; j++) {
                Excitor gemExcitor = MakeGeminalExcitor();
                gemExcitor.GetAlphaOrbitals() |= 1ULL << i;
                gemExcitor.GetBetaOrbitals() |= 1ULL << j;

                excitors.push_back(gemExcitor);
                nGeminalExcitors++;
            }

        //The populations of walkers on each of the excitors generated above. The MEG amplitudes are fixed so
        //their populations are never used.
        std::vector<int64_t> excitorWalkers(excitors.size(), 0);
        excitorWalkers[0] = INITIAL_REFERENCE_POPULATION;
        std::vector<int64_t> excitorSpawn(excitors.size(), 0);

        //The cumulative 1/2^(s+1) probability distribution. Clusters larger than the truncation level + 2 cannot
        //couple to the excitor space through the Hamiltonian
        std::vector<double> cumulativeClusterProb =
                GenerateCumulativeClusterProb(std::min(nElectrons, TRUNCATION_LEVEL + 2));

        //Energy shift
        double shift = 0.0;
        double tau = 0.01;
        bool varyShift = false;

        double EProjCycle = 0.0;
        double EProjAverage = 0.0;
        double EProjCycle_F12 = 0.0;
        double EProjAverage_F12 = 0.0;
        size_t nAveraged = 0;
        //Current total population of excips in the system
        int64_t nExcips = CountExcips(excitorWalkers);
        int64_t nExcipsLast = nExcips;

        fmt::print("\n\n| Cycle No.      | E(Proj_CCSD)/Eh| E(Proj_F12)/Eh | E(Proj Tot)/Eh "
            "| E(S)/Eh        | n(Excips)      | n(Excips D0)   |\n");

        for (size_t reportCycle = 0; reportCycle < NUM_REPORT_CYCLES; reportCycle++) {
            for (size_t cycle = 0; cycle < NUM_CYCLES; cycle++) {
                for (auto &val: excitorSpawn)
                    val = 0;

                nExcips = CountExcips(excitorWalkers);
                const double N0 = (double) excitorWalkers[0];

                double EProjCluster = 0.0;
                double EProjCluster_F12 = 0.0;
                //Number of cluster selection attempts this cycle
                int64_t Na = std::max<int64_t>(nExcips, 1);
                if (Na > 100000)
                    Na = 5000;

                for (int64_t attempt = 0; attempt < Na; attempt++) {
                    // -- Cluster selection --
                    // A cluster of size s is chosen with probability pSize. A MEG excitor may only be selected in a
                    // cluster of size one; larger clusters are formed from distinct conventional excitors.
                    double pSize = 1.0;
                    const uint32_t s = SampleClusterSize(randomNumber(), cumulativeClusterProb, pSize);

                    //Probability of selecting given excitors
                    double pClust = 1.0;
                    std::vector<uint32_t> selectedExcitorIndices;
                    if (s == 0) {
                        selectedExcitorIndices.push_back(0);
                    } else if (s == 1) {
                        uniform_int_distribution<uint32_t> excitorIdxDist(1, excitors.size() - 1);
                        selectedExcitorIndices.push_back(excitorIdxDist(randGenerator));
                        pClust = 1.0 / (double) (excitors.size() - 1);
                    } else {
                        if (s > nConventional - 1)
                            continue;

                        uniform_int_distribution<uint32_t> excitorIdxDist(1, nConventional - 1);
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

                        pClust = Integrals::Utils::Factorial(s) *
                                 pow(1.0 / (double) (nConventional - 1), s);
                    }

                    const double pSelect = pSize * pClust;

                    if (IsGeminalExcitor(excitors[selectedExcitorIndices[0]])) {
                        // -- MEG cluster --
                        // The MEG amplitude is fixed at one so the cluster carries N_0 excips
                        const ElectronPair gemPair = GetGeminalElectronPair(
                            excitors[selectedExcitorIndices[0]], nAlphaElectrons, nBetaElectrons);
                        const double clusterWeight = N0;

                        //Energy contribution of the MEG from the fixed amplitude MP2-F12 pair energy
                        EProjCluster_F12 += F12::CalculateMP2_F12_FixedEnergyElectronPair(
                            gemPair.first, gemPair.second,
                            F_pq, V_ijxy, X_xywz, B_xywz) / pSelect;

                        //Spawn onto a uniformly chosen conventional (non-reference) excitor
                        const double pExcit = 1.0 / double(nConventional - 1);
                        uniform_int_distribution<uint32_t> excitorIdxDist(1, nConventional - 1);
                        const uint32_t DjIdx = excitorIdxDist(randGenerator);

                        const double Hj_MEG = EvaluateGeminalElement(excitors[DjIdx] ^ D0,
                                                                     gemPair,
                                                                     V_pqxy,
                                                                     h_pqxy,
                                                                     nAlphaElectrons,
                                                                     nBetaElectrons,
                                                                     nSpatialOrbitals,
                                                                     nSpatialOrbitals);

                        double pSpawn = tau * abs(clusterWeight * Hj_MEG) / (Na * pSelect * pExcit);
                        int64_t nSpawn = static_cast<int64_t>(pSpawn);
                        pSpawn = pSpawn - (double) nSpawn;
                        if (pSpawn > randomNumber())
                            nSpawn++;

                        const int spawnSign = -boost::math::sign(Hj_MEG) * boost::math::sign(clusterWeight) *
                                              excitors[DjIdx].GetSign();
                        excitorSpawn[DjIdx] += spawnSign * nSpawn;

                        //Death/cloning is skipped for a MEG as its amplitude is fixed
                        continue;
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

                    //Clusters beyond the truncation level have no excitor but may still spawn into the excitor space
                    const uint32_t DiIdx = FindExcitorIdx(Di ^ D0, excitors);
                    const bool isInExcitorSpace = DiIdx < nConventional;

                    //Signed contribution of this cluster to the coefficient of |Di>, in excips
                    const double clusterWeight = clusterSign * N0 * clusterAmplitude;

                    //Update the projected energy
                    if (isInExcitorSpace && DiIdx != 0)
                        EProjCluster += clusterSign * clusterAmplitude * H_CI(0, DiIdx) / pSelect;

                    // -- Spawning --
                    // A random single or double excitation of |Di> is performed to generate |Dj>. The excitor
                    // of |Dj> changes by dN_j = -tau sign_j H_ji weight, where sign_j is the sign of the excitor
                    // acting on the reference
                    double pExcit = 0.5;
                    SlaterDet Dj = randomNumber() < 0.5
                                       ? RandomSingleExcite(Di, randGenerator, nSpatialOrbitals, pExcit)
                                       : RandomDoubleExcite(Di, randGenerator, nSpatialOrbitals, pExcit);
                    const uint32_t DjIdx = FindExcitorIdx(Dj ^ D0, excitors);

                    //Cannot spawn onto itself and check J is valid
                    if (!(Dj == Di) && DjIdx < nConventional) {
                        const double Hji = isInExcitorSpace
                                               ? H_CI(DjIdx, DiIdx)
                                               : CI::CalculateHamiltonianElement(Dj, Di,
                                                   nSpatialOrbitals, nSpatialOrbitals, HSpin, spinERI);

                        double pSpawn = tau * abs(clusterWeight * Hji) / (Na * pSelect * pExcit);
                        int64_t nSpawn = static_cast<int64_t>(pSpawn);
                        pSpawn = pSpawn - (double) nSpawn;
                        if (pSpawn > randomNumber())
                            nSpawn++;

                        const int spawnSign = -boost::math::sign(Hji) * boost::math::sign(clusterWeight) *
                                              excitors[DjIdx].GetSign();
                        excitorSpawn[DjIdx] += spawnSign * nSpawn;
                    }

                    if (!isInExcitorSpace)
                        continue;

                    // -- Death --
                    // dN_i = -tau sign_i (H_ii - H_00 - S) weight
                    const double Hii = H_CI(DiIdx, DiIdx) - H00;

                    double pCloneDeath = tau * abs(clusterWeight) * (Hii - shift) / (Na * pSelect);
                    int64_t nKill = static_cast<int64_t>(pCloneDeath);
                    pCloneDeath = pCloneDeath - static_cast<double>(nKill);
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
                for (uint32_t I = 0; I < nConventional; I++)
                    excitorWalkers[I] += excitorSpawn[I];

                if (excitorWalkers[0] == 0)
                    excitorWalkers[0] = INITIAL_REFERENCE_POPULATION;

                //Update the cycle averages for the projected energy
                EProjCycle += EProjCluster / static_cast<double>(Na);
                EProjCycle_F12 += EProjCluster_F12 / static_cast<double>(Na);
            }

            nExcips = CountExcips(excitorWalkers);

            double EProj = (EProjCycle / NUM_CYCLES) + refWfn.E_HF;
            double EProj_F12 = EProjCycle_F12 / NUM_CYCLES;

            EProjCycle = 0.0;
            EProjCycle_F12 = 0.0;

            if (!varyShift && nExcips >= TARGET_POPULATION)
                varyShift = true;
            else if (varyShift && nExcipsLast != 0 && nExcips != 0)
                shift -= SHIFT_DAMPING / (tau * NUM_CYCLES) * log((double) nExcips / (double) nExcipsLast);

            if (varyShift && reportCycle >= START_AVERAGE) {
                EProjAverage += EProj;
                EProjAverage_F12 += EProj_F12;
                nAveraged++;
            }

            fmt::print(" {:<16} {:<16.8} {:<16.8} {:<16.8} {:<16.8} {:<16} {:<16} \n",
                       reportCycle,
                       EProj,
                       EProj_F12,
                       EProj + EProj_F12,
                       shift,
                       nExcips,
                       excitorWalkers[0]);

            nExcipsLast = nExcips;
        }

        fmt::print("\n\n");
        PrintCCQMCWavefunction(std::vector<int64_t>(excitorWalkers.begin(),
                                                    excitorWalkers.begin() + nConventional));

        if (nAveraged == 0) {
            fmt::print(fg(fmt::color::crimson),
                       "!!! No projected energies were averaged (population never reached the target) !!!\n");
            return std::numeric_limits<double>::quiet_NaN();
        }

        return (EProjAverage + EProjAverage_F12) / double(nAveraged);
    }
    //
    //
#pragma endregion
    //
    //
}
