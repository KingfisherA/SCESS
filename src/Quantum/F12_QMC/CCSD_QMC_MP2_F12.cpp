#include <Quantum/Q_PCH.h>
//
#include <fmt/format.h>
//
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
                                  const Tensor4 &h_apxy,
                                  uint32_t nAlphaElectrons,
                                  uint32_t nBetaElectrons,
                                  uint32_t nOcc,
                                  uint32_t numOfAlphaOrbitals,
                                  uint32_t numOfBetaOrbitals) {
        auto *uniqueI = reinterpret_cast<uint16_t *>(
            alloca(std::max(numOfAlphaOrbitals, numOfBetaOrbitals) * sizeof(uint16_t)));
        auto *uniqueJ = reinterpret_cast<uint16_t *>(
            alloca(std::max(numOfAlphaOrbitals, numOfBetaOrbitals) * sizeof(uint16_t)));
        auto *commonIJ = reinterpret_cast<uint16_t *>(
            alloca(std::max(numOfAlphaOrbitals, numOfBetaOrbitals) * sizeof(uint16_t)));
        uint16_t uniqueISize, uniqueJSize, commonIJSize;
        int8_t sign;

        SlaterDet D0 = SlaterDet::MakeReference(nAlphaElectrons, nBetaElectrons);

        SlaterDet::CalculateDifferences(D0, Dj,
                                        numOfAlphaOrbitals, numOfBetaOrbitals,
                                        uniqueI, &uniqueISize,
                                        uniqueJ, &uniqueJSize,
                                        commonIJ, &commonIJSize,
                                        &sign);

        //Geminal electron 1 and 2
        auto [g1, g2] = gemPair;

        switch (uniqueISize) {
            case 0: {
                double sum = 0.0;

                for (uint32_t n = 0; n < commonIJSize / 2; ++n) {
                    uint32_t i = n * 2 + 0;
                    uint32_t j = n * 2 + 1;
                    sum += (3.0 / 8.0) * V_pqxy(commonIJ[i], commonIJ[j], g1, g2)
                            + (1.0 / 8.0) * V_pqxy(commonIJ[i], commonIJ[j], g1, g2);
                }

                return sum;
            }

            case 1: {
                double sum = 0.0;
                for (uint32_t i = 0; i < uniqueJSize / 2; ++i) {
                    int a = int(uniqueJ[i * 2 + 0]) / 2;
                    int b = int(uniqueJ[i * 2 + 1]) / 2;

                    sum += (3.0 / 8.0) * (V_pqxy(a, b, g1, g2) + h_apxy(a - nOcc, b, g1, g2))
                            + (1.0 / 8.0) * (V_pqxy(a, b, g1, g2) + h_apxy(a - nOcc, b, g1, g2));
                }

                return sum;
            }
            case 2: {
                double sum = 0.0;
                for (uint32_t i = 0; i < uniqueJSize / 2; ++i) {
                    int a = int(uniqueJ[i * 2 + 0]) / 2;
                    int b = int(uniqueJ[i * 2 + 1]) / 2;

                    sum += (3.0 / 8.0) * (V_pqxy(a, b, g1, g2) + h_apxy(a - nOcc, b, g1, g2))
                            + (1.0 / 8.0) * (V_pqxy(a, b, g1, g2) + h_apxy(a - nOcc, b, g1, g2));
                }

                return sum;
            }

            default: {
                return 0.0;
            }
        }
    }
    //
    std::pair<int, int> GetGeminalElectronPair(Excitor gemExcitor,
                                               uint32_t nAlphaOrbitals,
                                               uint32_t nBetaOrbitals) {
        std::pair<uint32_t, uint32_t> gemPair = {-1, -1};

        for (int i = 0; i < nAlphaOrbitals; i++)
            if ((0x1 << i) & gemExcitor.GetAlphaOrbitals()) {
                if (gemPair.first == -1)
                    gemPair.first = i;
                else {
                    gemPair.second = i;
                    return gemPair;
                }
            }

        for (int i = 0; i < nBetaOrbitals; i++)
            if ((0x1 << i) & gemExcitor.GetBetaOrbitals()) {
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
        //Number of attempts
        int64_t Na = 100;
        //Number of report cycles
        constexpr size_t NUM_REPORT_CYCLES = 400;
        //Cycle to start the population control shift
        constexpr int START_SHIFT = 5;
        //Cycle to start averaging the projected energy
        constexpr int START_AVERAGE = 200;
        //Target excip population
        constexpr int TARGET_POPULATION = 5000;


        TEITensor spinERI = Integrals::Transforms::AntisymmetriseIntegrateOutSpinClosedShell(
            Slice<Indices<P, Q, R, S> >(idxs,
                                        mints.eri));

        MatrixXd h = mints.T.block(0, 0, f12Idxs.nOBS, f12Idxs.nOBS) +
                     mints.V.block(0, 0, f12Idxs.nOBS, f12Idxs.nOBS);
        MatrixXd HSpin = CI::GenerateHCoreSpin(h);

        //Tensor4 W_ = EvaluateW(f12Ints.f12);

        fmt::print("Evaluating F12 T2'(1) Intermediates...\n");
        Tensor2 F_pq;
        TEITensor V_ijxy, C_xyab, X_xywz, B_xywz;
        F12::MP2_F12_Fixed(refWfn, mints, f12Ints, f12Idxs, &F_pq, &V_ijxy, &C_xyab, &X_xywz, &B_xywz);
        TEITensor V_pqxy, h_apxy;
        F12::EvaluateAdditionalMP2Intermediates(mints, f12Ints, f12Idxs,
                                                &V_pqxy, &h_apxy);
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

        //The first nAnnihilationBits number of bits indicate an annihilation of an occupied orbital
        uint32_t nAnnihilationBits = nOccupied;
        //The remaining bits in the excitor bit string are the creation operators for the virtual orbitals
        uint32_t nCreationBits = nSpatialOrbitals - nOccupied;

        mt19937 randGenerator(time(nullptr));
        uniform_01<mt19937> randomNumber(randGenerator);

        //The reference determinant
        const SlaterDet D0 = SlaterDet::MakeReference(nOccupied, nOccupied);
        //Generate all single and double excitors for CCSD
        //The first "nAnnihilationBits" bits are the annihilators and the
        //remaining "nCreationBits" are the creation operators
        //nElectrons
        std::vector<Excitor> excitors = EnumerateExcitors(D0, nSpinOrbs, nOccupied, nOccupied, 2);
        //The populations of walkers on each of the excitors generated above
        std::vector<int64_t> excitorWalkers(excitors.size(), 0);
        excitorWalkers[0] = 100;
        std::vector<int64_t> cumulativeWalkers(excitors.size(), 0);
        std::vector<int64_t> excitorSpawn(excitors.size(), 0);

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


        uint32_t nGeminalExcitors = 0;
        //Add the geminal (MEG) excitors
        for (uint32_t i = 0; i < nAlphaElectrons; i++)
            for (uint32_t j = 0; j < nBetaElectrons; j++) {
                Excitor gemExcitor = MakeGeminalExcitor();
                gemExcitor.GetAlphaOrbitals() |= 1 << i;
                gemExcitor.GetBetaOrbitals() |= 1 << j;

                excitors.push_back(gemExcitor);
                nGeminalExcitors++;
            }


        //The cumulative 1/2^(s+1) probability distribution
        std::vector<double> cumulativeClusterProb = GenerateCumulativeClusterProb(nElectrons);

        //Energy shift
        double shift = 0.0;
        double tau = 0.01;

        // Excitors; These are operators which take the reference determinant onto another
        // determinant such as |Di> = ai^\dag |D0> where ai is the excitor

        double EProjCycle = 0.0;
        double EProjAverage = 0.0;
        double EProjCycle_F12 = 0.0;
        double EProjAverage_F12 = 0.0;
        //Current total population of excips in the system
        int64_t nExcips = 0;
        int64_t nExcipsLast = 0;

        uint32_t n_E_F12_Cycle = 0;
        uint32_t n_E_Cycle = 0;

        int64_t D0_Population = excitorWalkers[0];

        fmt::print("\n\n| Cycle No.      | E(Proj_CCSD)/Eh| E(Proj_F12)/Eh | E(Proj Tot)/Eh "
            "| n(E)           | n(E_F12)       | E(S)/Eh        | n(Excips)      | n(Excips D0)   |\n");

        for (size_t reportCycle = 0; reportCycle < NUM_REPORT_CYCLES; reportCycle++) {
            for (size_t cycle = 0; cycle < NUM_CYCLES; cycle++) {
                for (auto &val: excitorSpawn)
                    val = 0;
                UpdateCumulativePop(cumulativeWalkers, excitorWalkers, nExcips);


                double EProjCluster = 0.0;
                double EProjCluster_F12 = 0.0;
                size_t n_E_F12_Cluster = 0;
                size_t n_E_Cluster = 0;
                Na = nExcips;
                if (abs(static_cast<double>(Na)) > 1E5)
                    Na = 5000;

                for (int64_t cluster = 0; cluster < Na; cluster++) {
                    // -- Cluster selection --
                    // Randomly choose a number of excitors to sample between 0 and the truncation level (l) + 2
                    // These excitators will be combined with each other and used to generate a given
                    // excited determinant |Di>
                    // This step is equivalent to truncating the exp(T) Taylor expansion at the l+2th power
                    // The 1/n! factor from the Taylor expansion is combined into the weight of the cluster

                    //Maximum excitation is l+2 where l is truncation of the exciation operator
                    //or the number of electrons
                    //uint32_t maxExcitation = std::min(nElectrons, truncationLevel);
                    uniform_int_distribution<uint32_t> clusterSizeDist(0, nElectrons);
                    //A random cluster is chosen of size s
                    uint32_t s = 0;

                    //Reject this cluster with probability 1.0 - 1/2^(s+1)
                    do {
                        s = clusterSizeDist(randGenerator);
                    } while (randomNumber() > cumulativeClusterProb[s]);

                    double pSize = 1.0 / pow(2.0, s + 1); //1.0 - cumulativeClusterProb[s - 1];

                    //Probability of selecting given excitors
                    double pClust = 1.0;
                    //Need to be able to select the same thing twice
                    std::unordered_set<uint32_t> selectedExcitorIndices;
                    if (s == 0) {
                        selectedExcitorIndices.insert(0);
                        pClust = 1.0;
                    } else if (s == 1) {
                        uniform_int_distribution<uint32_t> excitorIdxDist(0, excitors.size() - 1);
                        selectedExcitorIndices.insert(excitorIdxDist(randGenerator));
                        pClust = 1.0 / (double) (excitors.size() - 1);
                    } else {
                        uniform_int_distribution<uint32_t> excitorIdxDist(0, excitors.size() - 1 - nGeminalExcitors);
                        //uniform_01<mt19937> geminalExcitor(randGenerator);
                        selectedExcitorIndices.reserve(s);
                        while (selectedExcitorIndices.size() != s) {
                            uint32_t newIndx = excitorIdxDist(randGenerator);
                            if (selectedExcitorIndices.contains(newIndx)) {
                                //This attempt has failed given we can't choose the same index twice but we cannot just
                                //choose another index to ensure we can calculate the probability of selecting this
                                //precise cluster
                                //goto CLUSTER_LOOP_BEGIN;
                                cluster--;
                                continue;
                            }
                            selectedExcitorIndices.insert(newIndx);
                        }

                        pClust = Integrals::Utils::Factorial(s) *
                                 pow(1.0 / (double) (excitors.size() - nGeminalExcitors - 1), s);
                    }


                    //Collapse the excitors together taking care of sign
                    Excitor collapsedExcitor;
                    collapsedExcitor.GetSign() = 1.0;
                    double excitorAmplitude = 1.0;
                    bool isGeminal = false;
                    ElectronPair gemPair;

                    for (uint32_t idx: selectedExcitorIndices) {
                        if (IsGeminalExcitor(excitors[idx])) {
                            collapsedExcitor = excitors[idx];
                            gemPair = GetGeminalElectronPair(excitors[idx], nAlphaElectrons, nBetaElectrons);
                            isGeminal = true;
                            break;
                        }

                        collapsedExcitor.GetAlphaOrbitals() |= excitors[idx].GetAlphaOrbitals();
                        collapsedExcitor.GetBetaOrbitals() |= excitors[idx].GetBetaOrbitals();
                        collapsedExcitor.GetSign() *= excitors[idx].GetSign();

                        excitorAmplitude *= ((double) excitorWalkers[idx] / (double) D0_Population);
                    }

                    //Generate determinant by applying new excitor to the HF determinant
                    SlaterDet Di = D0;

                    if (!isGeminal) {
                        Di.GetAlphaOrbitals() = collapsedExcitor.GetAlphaOrbitals() ^ D0.GetAlphaOrbitals();
                        Di.GetBetaOrbitals() = collapsedExcitor.GetBetaOrbitals() ^ D0.GetBetaOrbitals();
                        Di.GetSign() = collapsedExcitor.GetSign();
                        // Di.CalculateSign(D0, nSpatialOrbitals, nSpatialOrbitals);

                        uint32_t numAlphaElecs = CountNumElecs(Di.GetAlphaOrbitals(), nSpatialOrbitals);
                        uint32_t numBetaElecs = CountNumElecs(Di.GetBetaOrbitals(), nSpatialOrbitals);

                        //We have formed an invalid  cluster
                        if (numAlphaElecs + numBetaElecs > nElectrons)
                            continue;
                    }

                    uint32_t DiIdx = FindExcitorIdx(Di ^ D0, excitors);
                    if (DiIdx == excitors.size()) continue;
                    uint32_t DjIdx = 0;

                    //Now randomly produce a determinant that is 1 or 2 excitations away from Dj
                    SlaterDet Dj;
                    //Probability of generating j given we have i
                    //Half chance of selecting a single excitation and half chance of selecting a double
                    double pExcit = 1.0; //0.5;
                    if (isGeminal) {
                        pExcit = 1.0 / double(excitors.size() - nGeminalExcitors);
                        uniform_int_distribution<uint32_t> excitorIdxDist(1, excitors.size() - nGeminalExcitors - 1);
                        DjIdx = excitorIdxDist(randGenerator);
                        Dj = excitors[DjIdx] ^ D0;
                    } else {
                        if (randomNumber() < 0.5) {
                            //Perform a single excitation/de-excitation
                            Dj = RandomSingleExcite(Di, randGenerator,
                                                    nSpatialOrbitals, nSpinOrbs, nElectrons,
                                                    nAnnihilationBits, nCreationBits, pExcit);
                        } else {
                            //Perform a double excitation/de-excitation
                            Dj = RandomDoubleExcite(Di, randGenerator,
                                                    nSpatialOrbitals, nSpinOrbs, nElectrons,
                                                    nAnnihilationBits, nCreationBits, pExcit);
                        }
                        DjIdx = FindExcitorIdx(Dj ^ D0, excitors);
                    }

                    //Update the projected energy
                    if (isGeminal) {
                        //Add contributions for each pair of electrons in the determinant
                        double H0i = F12::CalculateMP2_F12_FixedEnergyElectronPair(
                            gemPair.first, gemPair.second,
                            F_pq, V_ijxy, X_xywz, B_xywz);
                        EProjCluster_F12 += 1.0 * H0i;
                        n_E_F12_Cluster++;
                    } else if (DiIdx != 0) {
                        double H0i = H_CI(0, DiIdx);
                        EProjCluster += excitorAmplitude * H0i / pClust;
                        n_E_Cluster++;
                    }

                    // -- Spawning --
                    // A random single or double excitation of |Di> is performed to generate |Dj>
                    // We now spawn a particle onto |Dj> with probability proportional to Hij which corresponds to
                    // <Di | H | Dj>. The sign of the particle is equal to Hij

                    //Cannot spawn onto itself and check J is valid
                    if ((isGeminal || DiIdx != DjIdx) && DjIdx != excitors.size()) {
                        Dj.CalculateSign(D0, nSpatialOrbitals, nSpatialOrbitals);

                        double Hij = 0.0;
                        //int64_t IJ_Sign = 1;
                        if (isGeminal) {
                            Hij = 0.6 * EvaluateGeminalElement(Dj,
                                                               gemPair,
                                                               V_pqxy,
                                                               h_apxy,
                                                               nAlphaElectrons,
                                                               nBetaElectrons,
                                                               nOccupied,
                                                               nSpatialOrbitals,
                                                               nSpatialOrbitals);
                        } else {
                            Hij = H_CI(DiIdx, DjIdx);
                        }

                        double pSpawn = tau * abs(excitorAmplitude * Hij) / (Na * pClust * pExcit * pSize);
                        int64_t nSpawn = static_cast<int64_t>(pSpawn);
                        pSpawn = pSpawn - (double) nSpawn;
                        if (pSpawn > randomNumber())
                            nSpawn++;

                        if (nSpawn != 0) {
                            if (isGeminal) {
                                int gemSign = boost::math::sign(Hij) * (int) Dj.GetSign();

                                if (Hij > 0.0) {
                                    //Opposite sign to parent
                                    if (gemSign < 0)
                                        excitorSpawn[DjIdx] += nSpawn;
                                    else
                                        excitorSpawn[DjIdx] -= nSpawn;
                                } else {
                                    //Same sign as parent
                                    if (gemSign > 0)
                                        excitorSpawn[DjIdx] += nSpawn;
                                    else
                                        excitorSpawn[DjIdx] -= nSpawn;
                                }
                            } else {
                                if (Hij > 0.0) {
                                    //Opposite sign to parent
                                    if (excitorWalkers[DiIdx] < 0)
                                        excitorSpawn[DjIdx] += nSpawn;
                                    else
                                        excitorSpawn[DjIdx] -= nSpawn;
                                } else {
                                    //Same sign as parent
                                    if (excitorWalkers[DiIdx] > 0)
                                        excitorSpawn[DjIdx] += nSpawn;
                                    else
                                        excitorSpawn[DjIdx] -= nSpawn;
                                }
                            }
                        }
                    }

                    //If a geminal (MEG) was selected, skip the cloning/death stage as we are using
                    //fixed amplitudes and thus not sampling MEG amplitudes
                    if (isGeminal)
                        continue;

                    double Hii = H_CI(DiIdx, DiIdx);

                    // -- Death --
                    // A new particle is spawned back onto |Di> with the sign of |Hii-S| with sign opposite
                    // to the parent
                    double pCloneDeath = tau * abs(excitorAmplitude) * (Hii - shift) / (Na * pClust * pSize);
                    int64_t nKill = static_cast<int64_t>(pCloneDeath);
                    pCloneDeath = pCloneDeath - static_cast<double>(nKill);
                    if (abs(pCloneDeath) > randomNumber()) {
                        if (pCloneDeath > 0.0)
                            nKill++;
                        else
                            nKill--;
                    }

                    //nKill *= static_cast<int64_t>(Di.GetSign());

                    if (nKill != 0) {
                        if (excitorWalkers[DiIdx] > 0) {
                            excitorSpawn[DiIdx] -= nKill;
                        } else {
                            excitorSpawn[DiIdx] += nKill;
                        }
                    }
                }

                // -- Annihilation --
                // Pairs of particles with opposite signs on the same excitor are removed
                for (uint32_t I = 0; I < excitorWalkers.size(); I++)
                    excitorWalkers[I] += excitorSpawn[I];

                D0_Population += excitorSpawn[0];
                if (D0_Population == 0) {
                    D0_Population = 100;
                    excitorWalkers[0] = 100;
                }

                //Update the cycle averages for the projected energy
                if (n_E_Cluster != 0)
                    EProjCluster /= static_cast<double>(n_E_Cluster);

                EProjCycle += EProjCluster;
                EProjCycle_F12 += EProjCluster_F12;
                n_E_Cycle += n_E_Cluster;
                n_E_F12_Cycle += n_E_F12_Cluster;
            }

            double EProj = (EProjCycle / NUM_CYCLES) + refWfn.E_HF;
            double EProj_F12 = EProjCycle_F12 / (
                                   static_cast<double>(n_E_F12_Cycle) / static_cast<double>(nGeminalExcitors));

            EProjCycle = 0.0;
            EProjCycle_F12 = 0.0;

            if (reportCycle >= START_AVERAGE) {
                EProjAverage += EProj;
                EProjAverage_F12 += EProj_F12;
            }

            if (reportCycle == START_SHIFT)
                shift = refWfn.E_HF;

            if (reportCycle > START_SHIFT && nExcipsLast != 0 && nExcips != 0) {
                //(100.0 / (tau * 1000.0))
                shift -= (1.0 * tau) * log((double) nExcips / (double) nExcipsLast);
                shift -= (500.0 * tau) * log((double) nExcips / (double) TARGET_POPULATION);
            }

            fmt::print(" {:<16} {:<16.8} {:<16.8} {:<16.8} {:<16} {:<16} {:<16.8} {:<16} {:<16} \n",
                       reportCycle,
                       EProj,
                       EProj_F12,
                       EProj + EProj_F12,
                       n_E_Cycle,
                       n_E_F12_Cycle,
                       shift,
                       nExcips,
                       excitorWalkers[0]);


            n_E_Cycle = 0;
            n_E_F12_Cycle = 0;

            nExcipsLast = nExcips;
        }
        EProjAverage /= double(NUM_REPORT_CYCLES - START_AVERAGE);
        EProjAverage_F12 /= double(NUM_REPORT_CYCLES - START_AVERAGE);

        fmt::print("\n\n");
        PrintCCQMCWavefunction(excitorWalkers);

        return EProjAverage + EProjAverage_F12;
    }
    //
    //
#pragma endregion
    //
    //
}
