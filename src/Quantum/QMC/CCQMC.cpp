#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/math/special_functions/sign.hpp>
//
#include "Quantum/QMC/QMC.h"
#include "Quantum/CI/CI.h"
#include "GraphicsEngine/Tools/Threading.h"
#include "Quantum/Integrals/IntegralTransforms.h"

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
    //Return the number of annihilation and creation operators in a given bitstring
    void CountOperations(Excitor excitor,
                         const uint32_t nAnnihilationBits,
                         const uint32_t nCreationBits,
                         uint32_t &numAnnihAlpha,
                         uint32_t &numCreationAlpha,
                         uint32_t &numAnnihBeta,
                         uint32_t &numCreationBeta) {
        for (uint32_t i = 0; i < nAnnihilationBits + nCreationBits; ++i) {
            if (excitor.GetAlphaOrbitals() & (0x1 << i)) {
                if (i < nAnnihilationBits) {
                    numAnnihAlpha++;
                } else {
                    numCreationAlpha++;
                }
            }
        }

        for (uint32_t i = 0; i < nAnnihilationBits + nCreationBits; ++i) {
            if (excitor.GetBetaOrbitals() & (0x1 << i)) {
                if (i < nAnnihilationBits) {
                    numAnnihBeta++;
                } else {
                    numCreationBeta++;
                }
            }
        }
    }
    //
    uint32_t CountNumElecs(OrbBitStr bitStr, const uint32_t nSpatialOrbitals) {
        uint32_t numOrbs = 0;

        for (uint32_t i = 0; i < nSpatialOrbitals; i++) {
            if (bitStr & (0x1) << i)numOrbs++;
        }
        return numOrbs;
    }
    //
    void DestroyCreateElec(OrbBitStr &sourceBitStr,
                           OrbBitStr &destBitStr,
                           uint32_t orbitalToDelete,
                           uint32_t orbitalToCreate,
                           const uint32_t nAnnihilationBits,
                           const uint32_t nCreationBits) {
        uint32_t bitCounter = 0;

        OrbBitStr srcString = sourceBitStr;
        OrbBitStr destString = destBitStr;

        for (uint32_t i = 0; i < nAnnihilationBits + nCreationBits; i++) {
            if (srcString & (0x1 << i)) {
                if (orbitalToDelete == bitCounter) {
                    sourceBitStr &= ~(0x1 << i);
                    break;
                }

                bitCounter++;
            }
        }

        bitCounter = 0;
        for (uint32_t i = 0; i < nAnnihilationBits + nCreationBits; i++) {
            if ((destString & (0x1 << i)) == 0) {
                if (orbitalToCreate == bitCounter) {
                    destBitStr |= (0x1 << i);
                    break;
                }

                bitCounter++;
            }
        }
    }
    //
    SlaterDet RandomSingleExcite(SlaterDet I,
                                 mt19937 &randGenerator,
                                 const uint32_t nSpatialOrbitals,
                                 const uint32_t nSpinOrbs,
                                 const uint32_t nElectrons,
                                 const uint32_t nAnnihilationBits,
                                 const uint32_t nCreationBits,
                                 double &pExcit) {
        uniform_01<mt19937> random(randGenerator);

        //Randomly select whether to excite the alpha or beta electron
        OrbBitStr *sourceBitStr = nullptr;
        OrbBitStr *destBitStr = nullptr;

        if (random() > 0.5) {
            sourceBitStr = &I.GetAlphaOrbitals();
            destBitStr = &I.GetAlphaOrbitals();
        } else {
            sourceBitStr = &I.GetBetaOrbitals();
            destBitStr = &I.GetBetaOrbitals();
        }

        //pExcit *= 0.5;

        //Number of occupied spaces
        uint32_t numSrc = CountNumElecs(*sourceBitStr, nSpatialOrbitals);
        //Number of free spaces
        uint32_t numDest = nSpatialOrbitals - CountNumElecs(*destBitStr, nSpatialOrbitals);

        uniform_int_distribution<uint32_t> deleteDist(0, numSrc - 1);
        uniform_int_distribution<uint32_t> createDist(0, numDest - 1);

        if (numSrc == 0 || numDest == 0)
            return I;

        //Index of the orbital to destroy from the occupied set
        uint32_t orbitalToDelete = deleteDist(randGenerator);
        //Index of the orbital to create from the free set
        uint32_t orbitalToCreate = createDist(randGenerator);

        if (orbitalToCreate != orbitalToDelete)
            DestroyCreateElec(*sourceBitStr, *destBitStr,
                              orbitalToDelete, orbitalToCreate,
                              nAnnihilationBits, nCreationBits);

        //pExcit *= (1.0 / numSrc) * (1.0 / numDest);

        pExcit *= (1.0 / (nSpinOrbs - nElectrons));

        return I;
    }
    SlaterDet RandomDoubleExcite(SlaterDet I,
                                 mt19937 &randGenerator,
                                 const uint32_t nSpatialOrbitals,
                                 const uint32_t nSpinOrbs,
                                 const uint32_t nElectrons,
                                 const uint32_t nAnnihilationBits,
                                 const uint32_t nCreationBits,
                                 double &pExcit) {
        //Randomly select whether to excite the alpha or beta electron
        uniform_01<mt19937> random(randGenerator);

        //Cannot perform a double excitation if there are less than two electrons
        if (nElectrons < 2)
            return I;

        //Number of occupied spaces
        uint32_t numAlpha = CountNumElecs(I.GetAlphaOrbitals(), nSpatialOrbitals);
        //Number of free spaces
        uint32_t numBeta = CountNumElecs(I.GetBetaOrbitals(), nSpatialOrbitals);

        OrbBitStr *sourceBitStr1 = nullptr;
        OrbBitStr *destBitStr1 = nullptr;

        if (random() < 0.5) {
            sourceBitStr1 = &I.GetAlphaOrbitals();
            destBitStr1 = &I.GetAlphaOrbitals();
        } else {
            sourceBitStr1 = &I.GetBetaOrbitals();
            destBitStr1 = &I.GetBetaOrbitals();
        }

        OrbBitStr *sourceBitStr2 = nullptr;
        OrbBitStr *destBitStr2 = nullptr;

        if (random() < 0.5) {
            sourceBitStr2 = &I.GetAlphaOrbitals();
            destBitStr2 = &I.GetAlphaOrbitals();
        } else {
            sourceBitStr2 = &I.GetBetaOrbitals();
            destBitStr2 = &I.GetBetaOrbitals();
        }

        //0.5*0.5
        //pExcit *= 0.25;
        pExcit *= 2.0 * (1.0 / ((nSpinOrbs - nElectrons) / 2));

        //Number of occupied spaces
        uint32_t numSrc1 = CountNumElecs(*sourceBitStr1, nSpatialOrbitals);
        uint32_t numSrc2 = CountNumElecs(*sourceBitStr2, nSpatialOrbitals);
        //Number of free spaces
        uint32_t numDest1 = nSpatialOrbitals - CountNumElecs(*destBitStr1, nSpatialOrbitals);
        uint32_t numDest2 = nSpatialOrbitals - CountNumElecs(*destBitStr2, nSpatialOrbitals);

        if (numSrc1 == 0 || numSrc2 == 0 || numDest1 == 0 || numDest2 == 0)
            return I;

        uniform_int_distribution<uint32_t> deleteDist1(0, numSrc1 - 1);
        uniform_int_distribution<uint32_t> deleteDist2(0, numSrc2 - 1);
        uniform_int_distribution<uint32_t> createDist1(0, numDest1 - 1);
        uniform_int_distribution<uint32_t> createDist2(0, numDest2 - 1);

        //Index of the orbital to destroy from the occupied set
        uint32_t orbitalToDelete1 = deleteDist1(randGenerator);
        uint32_t orbitalToDelete2 = deleteDist2(randGenerator);
        //Index of the orbital to create from the free set
        uint32_t orbitalToCreate1 = createDist1(randGenerator);
        uint32_t orbitalToCreate2 = createDist2(randGenerator);

        DestroyCreateElec(*sourceBitStr1, *destBitStr1,
                          orbitalToDelete1, orbitalToCreate1,
                          nAnnihilationBits, nCreationBits);
        DestroyCreateElec(*sourceBitStr2, *destBitStr2,
                          orbitalToDelete2, orbitalToCreate2,
                          nAnnihilationBits, nCreationBits);

        //(1.0 / numSrc1) * (1.0 / numDest1) * (1.0 / numSrc2) * (1.0 / numDest2);

        return I;
    }
    //
    //
    uint32_t FindExcitorIdx(Excitor ex, const std::vector<Excitor> &excitors) {
        return std::distance(excitors.begin(), std::ranges::find(excitors, ex));
    }
    //
    std::vector<Excitor> EnumerateExcitors(SlaterDet reference, uint32_t nSpinOrbs, uint32_t nAlphaElectrons,
                                           uint32_t nBetaElectrons, uint32_t excitationLevel) {
        std::vector<Excitor> excitors = CI::EnumerateDet(nSpinOrbs, nAlphaElectrons, nBetaElectrons, excitationLevel);

        for (Excitor &ex: excitors) {
            ex.CalculateSign(reference, nSpinOrbs / 2, nSpinOrbs / 2);
            ex.GetAlphaOrbitals() ^= reference.GetAlphaOrbitals();
            ex.GetBetaOrbitals() ^= reference.GetBetaOrbitals();
        }

        return excitors;
    }
    //
    std::vector<Excitor> EnumerateAllExcitors(SlaterDet reference,
                                              uint32_t nSpinOrbs,
                                              uint32_t nElectrons) {
        std::vector<Excitor> excitors = CI::EnumerateAllDet(nSpinOrbs, nElectrons);

        for (Excitor &ex: excitors) {
            ex.CalculateSign(reference, nSpinOrbs / 2, nSpinOrbs / 2);
            ex.GetAlphaOrbitals() ^= reference.GetAlphaOrbitals();
            ex.GetBetaOrbitals() ^= reference.GetBetaOrbitals();
        }

        return excitors;
    }
    //
    void UpdateCumulativePop(std::vector<int64_t> &cumulativeWalkers,
                             std::vector<int64_t> &excitorWalkers,
                             int64_t &nExcips) {
        cumulativeWalkers[0] = excitorWalkers[0];
        nExcips = std::abs(excitorWalkers[0]);
        for (size_t i = 1; i < cumulativeWalkers.size(); i++) {
            cumulativeWalkers[i] = excitorWalkers[i] + cumulativeWalkers[i - 1];
            nExcips += std::abs(excitorWalkers[i]);
        }
    }
    //
    std::vector<double> GenerateCumulativeClusterProb(uint32_t nElectrons) {
        //The cumulative 1/2^(s+1) probability distribution
        std::vector<double> cumulativeClusterProb(nElectrons, 0.0);
        cumulativeClusterProb[0] = 0.5;
        for (uint32_t i = 1; i < cumulativeClusterProb.size(); i++) {
            cumulativeClusterProb[i] = cumulativeClusterProb[i - 1] + 1.0 / pow(2.0, i + 1);
        }
        cumulativeClusterProb.back() -= 1.0 - cumulativeClusterProb.back();

        return cumulativeClusterProb;
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
                 const GEng::Thread *thread,
                 boost::circular_buffer<glm::vec2> *projEnergyGraph) {
        constexpr size_t NUM_CYCLES = 100;
        //Number of attempts
        uint32_t Na = 200;
        constexpr size_t NUM_REPORT_CYCLES = 2000;

        std::cout << "--- Coupled Cluster Quantum Monte Carlo ---\n";

        TEITensor spinERI = Integrals::Transforms::AntisymmetriseIntegrateOutSpinClosedShell(
            Integrals::Transforms::TransformTEI(mints.eri, rhfWfn.C_Alpha));

        auto h = mints.T + mints.V;
        MatrixXd hSpin = CI::GenerateHCoreSpinTransform(h, rhfWfn.C_Alpha);

        //The number of doubly occupied orbitals
        //Given we are starting with a RHF determinant, the number of alpha=beta which will equal this
        //quantity.
        uint32_t nOccupied = std::distance(rhfWfn.Occ_Alpha.begin(),
                                           std::ranges::find(rhfWfn.Occ_Alpha, 0));
        uint32_t nElectrons = nOccupied * 2;
        uint32_t nSpatialOrbitals = rhfWfn.Occ_Alpha.size();
        uint32_t nSpinOrbs = nSpatialOrbitals * 2;

        //The first nAnnihilationBits number of bits indicate an annihilation of an occupied orbital
        uint32_t nAnnihilationBits = nOccupied;
        //The remaining bits in the excitor bit string are the creation operators for the virtual orbitals
        uint32_t nCreationBits = nSpatialOrbitals - nOccupied;

        mt19937 randGenerator(time(nullptr));
        uniform_01<mt19937> randomNumber(randGenerator);

        //The reference determinant
        const SlaterDet D0 = SlaterDet::MakeReference(nOccupied, nOccupied);
        //Generate all possible excitors in this full CC
        //The first "nAnnihilationBits" bits are the annihilators and the
        //remaning "nCreationBits" are the creation operators
        std::vector<Excitor> excitors = EnumerateAllExcitors(D0, nSpinOrbs, nElectrons);
        //The populations of walkers on each of the excitors generated above
        std::vector<int64_t> excitorWalkers(excitors.size(), 0);
        excitorWalkers[0] = 100;
        std::vector<int64_t> cumulativeWalkers(excitors.size(), 0);
        std::vector<int64_t> excitorSpawn(excitors.size(), 0);

        //Make the FCI Hamiltonian matrix
        MatrixXd H_FCI = MatrixXd::Zero(excitors.size(), excitors.size());
        for (uint32_t i = 0; i < H_FCI.rows(); i++)
            for (uint32_t j = i; j < H_FCI.cols(); j++)
                H_FCI(i, j) = H_FCI(j, i) =
                              CI::CalculateHamiltonianElement(
                                  excitors[i] ^ D0, excitors[j] ^ D0, nSpatialOrbitals,
                                  nSpatialOrbitals, hSpin, spinERI);

        //The cumulative 1/2^(s+1) probability distribution
        std::vector<double> cumulativeClusterProb = GenerateCumulativeClusterProb(nElectrons);

        double shift = 0.0;
        double tau = 0.1;

        // Excitors; These are operators which take the reference determinant onto another
        // determinant such as |Di> = ai^\dag |D0> where ai is the excitor

        double EProjCycle = 0.0;
        double EProjAverage = 0.0;
        //Current total population of excips in the system
        int64_t nExcips = 0;
        int64_t nExcipsLast = 0;

        int64_t D0_Population = excitorWalkers[0];

        for (size_t reportCycle = 0; reportCycle < NUM_REPORT_CYCLES; reportCycle++) {
            EProjCycle = 0.0;

            for (size_t cycle = 0; cycle < NUM_CYCLES; cycle++) {
                if (thread)
                    if (thread->ShouldStop())
                        return 0.0;

                for (auto &val: excitorSpawn)
                    val = 0;

                UpdateCumulativePop(cumulativeWalkers, excitorWalkers, nExcips);

                double EProjCluster = 0.0;
                Na = nExcips;

                for (size_t cluster = 0; cluster < Na; cluster++) {
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
                    if (s > 0) {
                        //uniform_int_distribution<uint32_t> excitorPopDist(0, totalPopulation);
                        uniform_int_distribution<uint32_t> excitorIdxDist(1, excitors.size() - 1);
                        selectedExcitorIndices.reserve(s);
                        while (selectedExcitorIndices.size() != s) {
                            uint32_t newIndx = excitorIdxDist(randGenerator);
                            if (selectedExcitorIndices.contains(newIndx)) {
                                //This attempt has failed given we can't choose the same index twice but we cannot just
                                //choose another index to ensure we can calculate the probability of selecting this
                                //precise cluster
                                //goto CLUSTER_LOOP_BEGIN;
                                continue;
                            }
                            selectedExcitorIndices.insert(newIndx);
                        }
                        pClust = Integrals::Utils::Factorial(s) * pow(1.0 / excitors.size(), s);
                    } else
                        selectedExcitorIndices.insert(0);

                    //Collapse the excitors together taking care of sign
                    Excitor collapsedExcitor;
                    collapsedExcitor.GetSign() = 1.0;
                    double excitorAmplitude = 1.0;

                    for (uint32_t idx: selectedExcitorIndices) {
                        collapsedExcitor.GetAlphaOrbitals() |= excitors[idx].GetAlphaOrbitals();
                        collapsedExcitor.GetBetaOrbitals() |= excitors[idx].GetBetaOrbitals();
                        collapsedExcitor.GetSign() *= excitors[idx].GetSign();

                        excitorAmplitude *= ((double) excitorWalkers[idx] / (double) D0_Population);
                    }

                    //Generate determinant by applying new excitor to the HF determinant
                    SlaterDet Di;
                    Di.GetAlphaOrbitals() = collapsedExcitor.GetAlphaOrbitals() ^ D0.GetAlphaOrbitals();
                    Di.GetBetaOrbitals() = collapsedExcitor.GetBetaOrbitals() ^ D0.GetBetaOrbitals();
                    Di.GetSign() = collapsedExcitor.GetSign();
                    // Di.CalculateSign(D0, nSpatialOrbitals, nSpatialOrbitals);

                    uint32_t numAlphaElecs = CountNumElecs(Di.GetAlphaOrbitals(), nSpatialOrbitals);
                    uint32_t numBetaElecs = CountNumElecs(Di.GetBetaOrbitals(), nSpatialOrbitals);

                    //We have formed a cluster
                    if (numAlphaElecs + numBetaElecs > nElectrons)
                        continue;
                    //goto CLUSTER_LOOP_BEGIN;

                    uint32_t DiIdx = FindExcitorIdx(Di ^ D0, excitors);
                    if (DiIdx == excitors.size()) continue;
                    uint32_t DjIdx = 0;

                    //Now randomly produce a determinant that is 1 or 2 excitations away from Dj
                    SlaterDet Dj;
                    //Probability of generating j given we have i
                    //Half chance of selecting a single excitation and half chance of selecting a double
                    double pExcit = 1.0; //0.5;
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

                    //Update the projected energy
                    if (DiIdx != 0) {
                        double H0i = H_FCI(0, DiIdx);
                        //static_cast<double>(Di.GetSign()) *
                        EProjCluster +=
                                excitorAmplitude * H0i / (pClust);
                    }

                    // -- Spawning --
                    // A random single or double excitation of |Di> is performed to generate |Dj>
                    // We now spawn a particle onto |Dj> with probability proportional to Hij which corresponds to
                    // <Di | H | Dj>. The sign of the particle is equal to Hij

                    //Cannot spawn onto itself and check J is valid
                    if (DiIdx != DjIdx && DjIdx != excitors.size()) {
                        Dj.CalculateSign(D0, nSpatialOrbitals, nSpatialOrbitals);

                        double Hij = H_FCI(DiIdx, DjIdx);

                        //1.2 *
                        double pSpawn = 4.0 * tau * abs(excitorAmplitude * Hij) / (Na * pClust * pExcit * pSize); //
                        auto nSpawn = static_cast<int64_t>(pSpawn);
                        pSpawn = pSpawn - (double) nSpawn;
                        if (pSpawn > randomNumber())
                            nSpawn++;

                        nSpawn *= static_cast<int64_t>(Dj.GetSign() * Di.GetSign());
                        if (-Hij > 0.0) {
                            excitorSpawn[DjIdx] += nSpawn;
                        } else {
                            excitorSpawn[DjIdx] -= nSpawn;
                        }
                    }

                    double Hii = H_FCI(DiIdx, DiIdx);

                    // -- Death --
                    // A new particle is spawned back onto |Di> with the sign of |Hii-S| with sign opposite
                    // to the parent
                    double pCloneDeath = tau * abs(excitorAmplitude) * (Hii - shift) / (Na * pClust);
                    //*pSize
                    auto nKill = static_cast<int64_t>(pCloneDeath);
                    pCloneDeath = pCloneDeath - (double) nKill;
                    if (abs(pCloneDeath) > randomNumber()) {
                        if (pCloneDeath > 0.0)
                            nKill++;
                        else
                            nKill--;
                    }

                    nKill *= static_cast<int64_t>(Di.GetSign());
                    if (excitorWalkers[DiIdx] > 0)
                        excitorSpawn[DiIdx] -= nKill;
                    else
                        excitorSpawn[DiIdx] += nKill;
                }

                EProjCluster /= (double) (Na);
                EProjCycle += EProjCluster;
                // -- Annihilation --
                // Pairs of particles with opposite signs on the same excitor are removed
                for (uint32_t I = 0; I < excitorWalkers.size(); I++)
                    excitorWalkers[I] += excitorSpawn[I];

                D0_Population += excitorSpawn[0];
                if (D0_Population == 0) {
                    D0_Population = 100;
                    excitorWalkers[0] = 100;
                }
            }

            //rhfWfn.Energy +
            //double EProj = (EProjCycle / NUM_CYCLES) + rhfWfn.Energy;

            double EProj = rhfWfn.E_HF;

            for (uint32_t I = 1; I < excitorWalkers.size(); I++) {
                EProj += H_FCI(0, I) * double(excitorWalkers[I]) / double(excitorWalkers[0]);
            }

            if (projEnergyGraph)
                projEnergyGraph->push_back(glm::vec2{(float) reportCycle, (float) EProj});

            if (reportCycle >= 50) {
                EProjAverage += EProj;
            }


            if (reportCycle == 50)
                shift = rhfWfn.E_HF;

            if (reportCycle > 50 && nExcipsLast != 0 && nExcips != 0)
                shift = shift - (100.0 / (tau * 1000.0)) * log((double) nExcips / (double) nExcipsLast);

            std::cout << "  Cycle(" << reportCycle << "); E(Proj) = " << EProj <<
                    " Shift = " << shift <<
                    " Population = " <<
                    nExcips
                    << " E(Proj) Average = " << EProjAverage / double(reportCycle - 50) << '\n';

            std::cout << "  Pop D0 = " << excitorWalkers[0] << "  Pop Dj = " << excitorWalkers.back() << '\n';

            nExcipsLast = nExcips;
        }
        EProjAverage /= double(NUM_REPORT_CYCLES - 50);


        return EProjAverage;
    }
    //
    //
#pragma endregion
    //
    //
}
