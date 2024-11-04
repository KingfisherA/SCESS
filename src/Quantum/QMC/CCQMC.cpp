#include  "Quantum/QMC/QMC.h"
#include  "Quantum/CI/ConfigurationInteraction.h"

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/random/uniform_int.hpp>

// Linked coupled-cluster; <Di |e^{-T} H e^{T}| D0>
// Unlinked coupled-cluster; <Di | H e^{T}| D0>

namespace SCESS {
    using namespace Eigen;
    using namespace boost::random;

    double QuantumMonteCarlo::CCQMC(const RHFWavefunction &rhfWfn, double Vnn) {
        constexpr size_t NUM_CYCLES = 20;
        constexpr size_t NUM_REPORT_CYCLES = 500000;

        std::cout << "--- Coupled Cluster Quantum Monte Carlo ---\n";

        TEITensor moSpinTEI = Integrals::Transforms::AntisymmetriseIntegrateOutSpinClosedShell(
            Integrals::Transforms::ToMolecularTEI(rhfWfn.Integrals.tei, rhfWfn.Coeffs));

        MatrixXd HSpin = ConfigurationInteraction::GenerateOneElectronHamiltonianSpin(rhfWfn);

        //The number of doubly occupied orbitals
        //Given we are starting with a RHF determinant, the number of alpha=beta which will equal this
        //quantity.
        uint32_t nOccupied = std::distance(rhfWfn.Occupation.begin(),
                                           std::ranges::find(rhfWfn.Occupation, 0));
        uint32_t nElectrons = nOccupied * 2;
        uint32_t nSpatialOrbitals = rhfWfn.Occupation.size();
        uint32_t nSpinOrbs = nSpatialOrbitals * 2;

        //The first nAnnihilationBits number of bits indicate an annihilation of an occupied orbital
        uint32_t nAnnihilationBits = nOccupied;
        //The remaining bits in the excitor bit string are the creation operators for the virtual orbitals
        uint32_t nCreationBits = nSpatialOrbitals - nOccupied;

        mt19937 randGenerator(time(nullptr));
        uniform_01<mt19937> randomNumber(randGenerator);

        //The reference determinant
        const CIDeterminant D0 = CIDeterminant::MakeReference(nOccupied, nOccupied);
        //Generate all possible excitors in this full CC
        //The first "nAnnihilationBits" bits are the annihilators and the
        //remaning "nCreationBits" are the creation operators
        std::vector<Excitor> excitors = GenerateAllExcitors(D0, nSpinOrbs, nElectrons);
        //The populations of walkers on each of the excitors generated above
        std::vector<int64_t> excitorWalkers(excitors.size(), 0);
        excitorWalkers[0] = 100;
        std::vector<int64_t> cumulativeWalkers(excitors.size(), 0);
        std::vector<int64_t> excitorSpawn(excitors.size(), 0);

        //Make the Hamiltonian matrix
        MatrixXd H_FCI = MatrixXd::Zero(excitors.size(), excitors.size());
        for (uint32_t i = 0; i < H_FCI.rows(); i++)
            for (uint32_t j = i; j < H_FCI.cols(); j++)
                H_FCI(i, j) = H_FCI(j, i) =
                              ConfigurationInteraction::CalculateHamiltonianElement(
                                  excitors[i] ^ D0, excitors[j] ^ D0, nSpatialOrbitals,
                                  nSpatialOrbitals, HSpin, moSpinTEI);


        //The cumulative 1/2^(s+1) probability distribution
        std::vector<double> cumulativeClusterProb(nElectrons + 3, 0.0);
        cumulativeClusterProb[0] = 0.5;
        for (uint32_t i = 1; i < cumulativeClusterProb.size(); i++)
            cumulativeClusterProb[i] = cumulativeClusterProb[i - 1] + 1.0 / pow(2.0, i + 1);

        //Return the number of annihilation and creation operators in a given bitstring
        auto CountOperations = [&nAnnihilationBits, &nCreationBits](Excitor excitor,
                                                                    uint32_t &numAnnihAlpha,
                                                                    uint32_t &numCreationAlpha,
                                                                    uint32_t &numAnnihBeta,
                                                                    uint32_t &numCreationBeta) {
            for (uint32_t i = 0; i < nAnnihilationBits + nCreationBits; ++i) {
                if (excitor.GetAlphaOrbitals() & (0x1 << i)) {
                    if (i < nAnnihilationBits) {
                        numAnnihAlpha++;
                    }
                    else {
                        numCreationAlpha++;
                    }
                }
            }

            for (uint32_t i = 0; i < nAnnihilationBits + nCreationBits; ++i) {
                if (excitor.GetBetaOrbitals() & (0x1 << i)) {
                    if (i < nAnnihilationBits) {
                        numAnnihBeta++;
                    }
                    else {
                        numCreationBeta++;
                    }
                }
            }
        };

        auto RandomSingleExcite = [&nAnnihilationBits, &nCreationBits, &randGenerator, &randomNumber](CIDeterminant I) {
            // uniform_int_distribution<uint32_t> dist(0, nAnnihilationBits + nCreationBits - 1);
            bool isAlpha = false;

            if (I.GetAlphaOrbitals() > 0 && I.GetBetaOrbitals() > 0) {
                isAlpha = randomNumber() > 0.5;
            } else if (I.GetAlphaOrbitals() > 0) {
                isAlpha = true;
            } else {
                isAlpha = false;
            }

            //Randomly select whether to excite the alpha or beta electron
            Excitor::OrbitalBitStr &bitStr = isAlpha ? I.GetAlphaOrbitals() : I.GetBetaOrbitals();
            // uint32_t bitToAnnihilate = 0, bitToCreate = 0, attempts = 0;
            //
            // do {
            //     bitToAnnihilate = dist(randGenerator);
            //     bitToCreate = dist(randGenerator);
            //
            //     if (attempts > 15) {
            //         bitToAnnihilate = 0;
            //         bitToCreate = 0;
            //         break;
            //     }
            //
            //     attempts++;
            // } while (!(bitStr & (0x1 << bitToAnnihilate)
            //            && !(bitStr & (0x1 << bitToCreate))));
            //
            // bitStr ^= (0x1 << bitToAnnihilate);
            // bitStr |= (0x1 << bitToCreate);

            if (randomNumber() < 0.5) {
                bitStr = bitStr << 1;
                if (bitStr & (0x1 << (nAnnihilationBits + nCreationBits))) {
                    bitStr |= 0x1;
                    bitStr ^= 0x1 << (nAnnihilationBits + nCreationBits);
                }
            } else {
                bool hasShiftedOffEnd = bitStr & 0x1;
                bitStr = bitStr >> 1;
                if (hasShiftedOffEnd)
                    bitStr |= 0x1 << (nAnnihilationBits + nCreationBits - 1);
            }

            I.GetSign() *= -1;
            return I;
        };

        auto FindExcitorIdx = [&excitors](Excitor ex)-> uint32_t {
            return (uint32_t) std::distance(excitors.begin(),
                                            std::find(excitors.begin(), excitors.end(), ex));
        };


        double shift = 0.0;
        double tau = 0.001;

        // Excitors; These are operators which take the reference determinant onto another
        // determinant such as |Di> = ai^\dag |D0> where ai is the excitor

        double EProjCycle = 0.0;
        double EProjAverage = 0.0;
        int64_t totalPopulation = 0;
        int64_t totalPopLast = 0;

        auto UpdateCumulativePop = [&cumulativeWalkers, &excitorWalkers, &totalPopulation]() {
            cumulativeWalkers[0] = excitorWalkers[0];
            for (size_t i = 1; i < cumulativeWalkers.size(); i++)
                cumulativeWalkers[i] = excitorWalkers[i] + cumulativeWalkers[i - 1];

            totalPopulation = cumulativeWalkers.back();
        };

        for (size_t reportCycle = 0; reportCycle < NUM_REPORT_CYCLES; reportCycle++) {
            for (size_t cycle = 0; cycle < NUM_CYCLES; cycle++) {
                for (auto &val: excitorSpawn)
                    val = 0;
                UpdateCumulativePop();

                double EProjCluster = 0.0;

                for (size_t cluster = 0; cluster < abs(totalPopulation); cluster++) {
                CLUSTER_LOOP_BEGIN:
                    // -- Cluster selection --
                    // Randomly choose a number of excitors to sample between 0 and the truncation level (l) + 2
                    // These excitators will be combined with each other and used to generate a given
                    // excited determinant |Di>
                    // This step is equivalent to truncating the exp(T) Taylor expansion at the l+2th power
                    // The 1/n! factor from the Taylor expansion is combined into the weight of the cluster

                    //A random cluster is chosen of size s
                    uniform_int_distribution<uint32_t> clusterSizeDist(0, nElectrons);
                    uint32_t s = clusterSizeDist(randGenerator);

                    //Reject this cluster with probability 1.0 - 1/2^(s+1)
                    while (randomNumber() > cumulativeClusterProb[s]) {
                        s = clusterSizeDist(randGenerator);
                    }

                    //Need to be able to select the same thing twice
                    std::unordered_set<uint32_t> selectedExcitorIndices;
                    if (s != 0) {
                        //uniform_int_distribution<uint32_t> excitorPopDist(0, totalPopulation);
                        uniform_int_distribution<uint32_t> excitorIdxDist(0, excitors.size() - 1);
                        selectedExcitorIndices.reserve(s);
                        while (selectedExcitorIndices.size() != s) {
                            uint32_t newIndx = excitorIdxDist(randGenerator);
                            if (selectedExcitorIndices.contains(newIndx)) {
                                //This attempt has failed given we can't choose the same index twice but we cannot just
                                //choose another index to ensure we can calculate the probability of selecting this
                                //precise cluster
                                goto CLUSTER_LOOP_BEGIN;
                            }
                            selectedExcitorIndices.insert(newIndx);
                        }
                    } else
                        selectedExcitorIndices.insert(0);

                    //Collapse the excitors together taking care of sign
                    Excitor collapsedExcitor;
                    uint32_t numAnniAlpha = 0,
                            numCreationAlpha = 0,
                            numAnnihBeta = 0,
                            numCreationBeta = 0;
                    double excitorAmplitude = 1.0;

                    for (uint32_t idx: selectedExcitorIndices) {
                        collapsedExcitor.GetAlphaOrbitals() |= excitors[idx].GetAlphaOrbitals();
                        collapsedExcitor.GetBetaOrbitals() |= excitors[idx].GetBetaOrbitals();
                        collapsedExcitor.GetSign() *= excitors[idx].GetSign();

                        excitorAmplitude *= ((double) excitorWalkers[idx] / (double) excitorWalkers[0]);
                    }

                    CountOperations(collapsedExcitor, numAnniAlpha,
                                    numCreationAlpha,
                                    numAnnihBeta, numCreationBeta);
                    //We cannot excite or annihilate more than the number of electrons we have
                    if ((numAnniAlpha + numAnnihBeta) / 2 > nElectrons) continue;

                    //Generate determinant by applying new excitor to the HF determinant
                    CIDeterminant Di;
                    Di.GetAlphaOrbitals() = collapsedExcitor.GetAlphaOrbitals() ^ D0.GetAlphaOrbitals();
                    Di.GetBetaOrbitals() = collapsedExcitor.GetBetaOrbitals() ^ D0.GetBetaOrbitals();
                    Di.GetSign() = collapsedExcitor.GetSign();

                    uint32_t DiIdx = FindExcitorIdx(Di ^ D0);
                    if (DiIdx == excitors.size()) continue;
                    uint32_t DjIdx = 0;

                    //Update the projected energy
                    double H0i = H_FCI(0, DiIdx);
                    EProjCluster += excitorAmplitude * H0i;

                    //Now randomly produce a determinant that is 1 or 2 excitations away from Dj
                    CIDeterminant Dj;
                    if (randomNumber() < 0.5)
                    //Perform a single excitation/de-excitation
                        Dj = RandomSingleExcite(Di);
                    else {
                        //Perform a double excitation/de-excitation
                        Dj = RandomSingleExcite(Di);
                        Dj = RandomSingleExcite(Dj);
                    }
                    DjIdx = FindExcitorIdx(Dj ^ D0);

                    // -- Spawning --
                    // A random single or double excitation of |Di> is performed to generate |Dj>
                    // We now spawn a particle onto |Dj> with probability proportional to Hij which corresponds to
                    // <Di | H | Dj>. The sign of the particle is equal to Hij
                    double Hij = H_FCI(DiIdx, DjIdx);
                    double Hjj = H_FCI(DjIdx, DjIdx);

                    //Cannot spawn onto itself
                    if (DiIdx == DjIdx) continue;

                    // = Amplitude/(nA*pSelect)
                    //double xe = abs(excitorCoeff) / (NUM_ATTEMPTS * pSelect);
                    double xe = 1.0;
                    //((double) excitorWalkers[DiIdx] / (double) excitorWalkers[0]) / (totalPopulation * 1.0);

                    if (s > 2)
                        xe = ((pow(2.0, s + 1)) / boost::math::factorial<double>(s)) *
                             pow((double) totalPopulation / (double) excitorWalkers[0], s - 1)
                             * (1.0 / (1.0 + (double) excitorWalkers[0] / (double) totalPopulation));

                    //Probability of spawning onto j from i
                    double pSpawn = tau * abs(Hij * xe);
                    auto nSpawn = int64_t(pSpawn);
                    pSpawn = pSpawn - (double) nSpawn;
                    //nSpawn = 0;
                    if (pSpawn > randomNumber())
                        nSpawn++;

                    if (-(Hij * xe) > 0.0) {
                        excitorSpawn[DjIdx] += nSpawn;
                    } else {
                        excitorSpawn[DjIdx] -= nSpawn;
                    }


                    // -- Death --
                    // A new particle is spawned back onto |Di> with the sign of |Hii-S| with sign opposite
                    // to the parent
                    double pCloneDeath = tau * abs(xe * (Hjj - shift));
                    auto nKill = static_cast<int64_t>(pSpawn);
                    pCloneDeath = pCloneDeath - (double) nKill;
                    //nKill = 0;
                    if (abs(pCloneDeath) > randomNumber())
                        nKill++;

                    if ((-Hjj + shift) * xe > 0.0) {
                        excitorWalkers[DjIdx] += nKill;
                    } else {
                        excitorWalkers[DjIdx] -= nKill;
                    }
                }

                EProjCluster /= abs(totalPopulation);
                EProjCycle += EProjCluster;
                // -- Annihilation --
                // Pairs of particles with opposite signs on the same excitor are removed
                for (uint32_t I = 0; I < excitorWalkers.size(); I++)
                    excitorWalkers[I] += excitorSpawn[I];
            }

            //double EProj = rhfWfn.Energy;
            double EProj = (EProjCycle / (double) (NUM_CYCLES)) + Vnn;

            // for (uint32_t i = 1; i < excitorWalkers.size(); i++)
            //     EProj += H_FCI(0, i) * abs((double) excitorWalkers[i] / (double) excitorWalkers[0]);

            if (reportCycle >= 5000) {
                EProjAverage += EProj;
            }

            if (reportCycle >= 50)
                shift = shift - (1.0 / (tau * 100.0)) * log((double) totalPopulation / (double) totalPopLast);

            std::cout << "  Cycle(" << reportCycle << "); E(Proj) = " << EProj <<
                    " Shift = " << shift <<
                    " Population = " <<
                    totalPopulation
                    << " E(Proj) Average = " << EProjAverage / double(reportCycle + 1 - 5000) << '\n';

            EProjCycle = 0.0;
            totalPopLast = totalPopulation;
        }


        return EProjAverage / double(NUM_REPORT_CYCLES + 1 - 5000);
    }

    std::vector<Excitor> QuantumMonteCarlo::GenerateAllExcitors(CIDeterminant reference, uint32_t nSpinOrbs,
                                                                uint32_t nElectrons) {
        std::vector<Excitor> excitors =
                ConfigurationInteraction::GenerateAllDeterminants(nSpinOrbs, nElectrons);

        for (Excitor &ex: excitors) {
            ex.GetAlphaOrbitals() ^= reference.GetAlphaOrbitals();
            ex.GetBetaOrbitals() ^= reference.GetBetaOrbitals();

            uint32_t flips = 0;
            for (uint32_t i = 0; i < nSpinOrbs / 2; i++)
                if ((0x1 << i) & ex.GetAlphaOrbitals())
                    flips++;

            ex.GetSign() = 1;
            if (flips % 2 == 0)
                ex.GetSign() = -1;
        }

        return excitors;
    }
}
