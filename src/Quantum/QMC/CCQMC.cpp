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
        constexpr size_t NUM_CYCLES = 100;
        //Number of attempts
        uint32_t Na = 500;
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

        //To check the FCI hamiltonian matrix has been constructed correctly
        // SelfAdjointEigenSolver<MatrixXd> eigensolver(H_FCI);
        // VectorXd fciEnergies = eigensolver.eigenvalues().real();
        // MatrixXd fciWavefuncs = eigensolver.eigenvectors().real();

        //The cumulative 1/2^(s+1) probability distribution
        std::vector<double> cumulativeClusterProb(nElectrons, 0.0);
        cumulativeClusterProb[0] = 0.5;
        for (uint32_t i = 1; i < cumulativeClusterProb.size(); i++) {
            cumulativeClusterProb[i] = cumulativeClusterProb[i - 1] + 1.0 / pow(2.0, i + 1);
        }
        cumulativeClusterProb.back() -= 1.0 - cumulativeClusterProb.back();


        auto CountNumElecs = [&nSpatialOrbitals](CIDeterminant::OrbitalBitStr bitStr) {
            uint32_t numOrbs = 0;

            for (uint32_t i = 0; i < nSpatialOrbitals; i++) {
                if (bitStr & (0x1) << i)numOrbs++;
            }
            return numOrbs;
        };

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
        };

        auto RandomSingleExcite = [&nAnnihilationBits, &nCreationBits,&nOccupied,
                    &nSpatialOrbitals, &nElectrons, &randomNumber,&randGenerator, &CountNumElecs](CIDeterminant I) {
            //Randomly select whether to excite the alpha or beta electron

            //Number of occupied spaces
            uint32_t numAlpha = CountNumElecs(I.GetAlphaOrbitals());
            //Number of free spaces
            uint32_t numBeta = CountNumElecs(I.GetBetaOrbitals());
            Excitor::OrbitalBitStr *sourceBitStr = nullptr;
            Excitor::OrbitalBitStr *destBitStr = nullptr;

            if (numAlpha == nElectrons) {
                sourceBitStr = &I.GetAlphaOrbitals();
                destBitStr = &I.GetBetaOrbitals();
            } else if (numBeta == nElectrons) {
                sourceBitStr = &I.GetBetaOrbitals();
                destBitStr = &I.GetAlphaOrbitals();
            } else {
                sourceBitStr = randomNumber() > 0.5 ? &I.GetAlphaOrbitals() : &I.GetBetaOrbitals();
                destBitStr = randomNumber() > 0.5 ? &I.GetAlphaOrbitals() : &I.GetBetaOrbitals();
            }

            //Number of occupied spaces
            uint32_t numSrc = CountNumElecs(*sourceBitStr);
            //Number of free spaces
            uint32_t numDest = nSpatialOrbitals - CountNumElecs(*destBitStr);

            if (numSrc == 0) {
                numDest = nSpatialOrbitals;
                numSrc = nElectrons;
                std::swap(sourceBitStr, destBitStr);
            }

            uniform_int_distribution<uint32_t> deleteDist(0, numSrc - 1);
            uniform_int_distribution<uint32_t> createDist(0, numDest - 1);

            //Index of the orbital to destroy from the occupied set
            uint32_t orbitalToDelete = deleteDist(randGenerator);
            //Index of the orbital to create from the free set
            uint32_t orbitalToCreate = createDist(randGenerator);

            uint32_t bitCounter = 0;
            for (uint32_t i = 0; i < nAnnihilationBits + nCreationBits; i++) {
                if (*sourceBitStr & (0x1 << i)) {
                    if (orbitalToDelete == bitCounter) {
                        *sourceBitStr &= ~(0x1 << i);
                        break;
                    }

                    bitCounter++;
                }
            }

            bitCounter = 0;
            for (uint32_t i = 0; i < nAnnihilationBits + nCreationBits; i++) {
                if ((*destBitStr & (0x1 << i)) == 0) {
                    if (orbitalToCreate == bitCounter) {
                        *destBitStr |= (0x1 << i);
                        break;
                    }

                    bitCounter++;
                }
            }

            return I;
        };

        auto FindExcitorIdx = [&excitors](Excitor ex)-> uint32_t {
            return (uint32_t) std::distance(excitors.begin(),
                                            std::find(excitors.begin(), excitors.end(), ex));
        };


        double shift = 0.0;
        double tau = 0.1;

        // Excitors; These are operators which take the reference determinant onto another
        // determinant such as |Di> = ai^\dag |D0> where ai is the excitor

        double EProjCycle = 0.0;
        double EProjAverage = 0.0;
        //Current total population of excips in the system
        int64_t nExcips = 0;
        int64_t nExcipsLast = 0;

        auto UpdateCumulativePop = [&cumulativeWalkers, &excitorWalkers, &nExcips]() {
            cumulativeWalkers[0] = excitorWalkers[0];
            for (size_t i = 1; i < cumulativeWalkers.size(); i++)
                cumulativeWalkers[i] = excitorWalkers[i] + cumulativeWalkers[i - 1];

            nExcips = cumulativeWalkers.back();
        };

        for (size_t reportCycle = 0; reportCycle < NUM_REPORT_CYCLES; reportCycle++) {
            EProjCycle = 0.0;
            for (size_t cycle = 0; cycle < NUM_CYCLES; cycle++) {
                for (auto &val: excitorSpawn)
                    val = 0;
                UpdateCumulativePop();

                double EProjCluster = 0.0;

                for (size_t cluster = 0; cluster < Na; cluster++) {
                CLUSTER_LOOP_BEGIN:
                    // -- Cluster selection --
                    // Randomly choose a number of excitors to sample between 0 and the truncation level (l) + 2
                    // These excitators will be combined with each other and used to generate a given
                    // excited determinant |Di>
                    // This step is equivalent to truncating the exp(T) Taylor expansion at the l+2th power
                    // The 1/n! factor from the Taylor expansion is combined into the weight of the cluster

                    //Maximum excitation is l+2 where l is truncation of the exciation operator
                    //or the number of electrons
                    //uint32_t maxExcitation = std::min(nElectrons, truncationLevel);
                    uniform_int_distribution<uint32_t> clusterSizeDist(1, nElectrons);
                    //A random cluster is chosen of size s
                    uint32_t s = 0;

                    //Reject this cluster with probability 1.0 - 1/2^(s+1)
                    do {
                        s = clusterSizeDist(randGenerator);
                    } while (randomNumber() > cumulativeClusterProb[s - 1]);

                    //Need to be able to select the same thing twice
                    std::unordered_set<uint32_t> selectedExcitorIndices;
                    if (s > 1) {
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

                    if (nExcips == 0) {
                        excitorWalkers[0] = 100;
                        break;
                    }

                    double pSelect = boost::math::factorial<double>(s) * pow(1.0 / excitors.size(), s);

                    //Collapse the excitors together taking care of sign
                    Excitor collapsedExcitor;
                    // uint32_t numAnniAlpha = 0,
                    //         numCreationAlpha = 0,
                    //         numAnnihBeta = 0,
                    //         numCreationBeta = 0;
                    double excitorAmplitude = 1.0;

                    for (uint32_t idx: selectedExcitorIndices) {
                        collapsedExcitor.GetAlphaOrbitals() |= excitors[idx].GetAlphaOrbitals();
                        collapsedExcitor.GetBetaOrbitals() |= excitors[idx].GetBetaOrbitals();
                        collapsedExcitor.GetSign() *= excitors[idx].GetSign();

                        excitorAmplitude *= ((double) excitorWalkers[idx] / (double) excitorWalkers[0]);
                    }
                    excitorAmplitude *= collapsedExcitor.GetSign();

                    // CountOperations(collapsedExcitor, numAnniAlpha,
                    //                 numCreationAlpha,
                    //                 numAnnihBeta, numCreationBeta);
                    // //We cannot excite or annihilate more than the number of electrons we have
                    // if ((numAnniAlpha + numAnnihBeta) / 2 > nElectrons) goto CLUSTER_LOOP_BEGIN;

                    //Generate determinant by applying new excitor to the HF determinant
                    CIDeterminant Di;
                    Di.GetAlphaOrbitals() = collapsedExcitor.GetAlphaOrbitals() ^ D0.GetAlphaOrbitals();
                    Di.GetBetaOrbitals() = collapsedExcitor.GetBetaOrbitals() ^ D0.GetBetaOrbitals();
                    Di.GetSign() = collapsedExcitor.GetSign();

                    uint32_t DiIdx = FindExcitorIdx(Di ^ D0);
                    if (DiIdx == excitors.size()) goto CLUSTER_LOOP_BEGIN;
                    uint32_t DjIdx = 0;

                    //Now randomly produce a determinant that is 1 or 2 excitations away from Dj
                    CIDeterminant Dj;
                    //Probability of generating j given we have i
                    double pGen = 0.25 / nSpatialOrbitals;
                    if (randomNumber() < 0.5) {
                        //Perform a single excitation/de-excitation
                        Dj = RandomSingleExcite(Di);
                    } else {
                        //Perform a double excitation/de-excitation
                        Dj = RandomSingleExcite(Di);
                        Dj = RandomSingleExcite(Dj);
                        pGen = pGen * pGen;
                    }
                    DjIdx = FindExcitorIdx(Dj ^ D0);

                    //Update the projected energy
                    double H0i = H_FCI(0, DiIdx);
                    EProjCluster += excitorAmplitude * H0i / pSelect;

                    // -- Spawning --
                    // A random single or double excitation of |Di> is performed to generate |Dj>
                    // We now spawn a particle onto |Dj> with probability proportional to Hij which corresponds to
                    // <Di | H | Dj>. The sign of the particle is equal to Hij
                    double Hij = H_FCI(DiIdx, DjIdx);
                    double Hjj = H_FCI(DjIdx, DjIdx);

                    //Cannot spawn onto itself
                    if (DiIdx == DjIdx) goto CLUSTER_LOOP_BEGIN;

                    //Probability of spawning onto j from i
                    double pSpawn = tau * abs(excitorAmplitude) * abs(Hij) / (Na * pSelect * pGen);
                    auto nSpawn = int64_t(pSpawn);
                    pSpawn = pSpawn - (double) nSpawn;
                    //nSpawn = 0;
                    if (pSpawn > randomNumber())
                        nSpawn++;

                    if (excitorAmplitude > 0.0) {
                        excitorSpawn[DjIdx] -= nSpawn;
                    } else {
                        excitorSpawn[DjIdx] += nSpawn;
                    }


                    // -- Death --
                    // A new particle is spawned back onto |Di> with the sign of |Hii-S| with sign opposite
                    // to the parent
                    double pCloneDeath = tau * abs(excitorAmplitude) * abs(Hjj - shift) / (Na * pSelect);
                    auto nKill = static_cast<int64_t>(pSpawn);
                    pCloneDeath = pCloneDeath - (double) nKill;
                    //nKill = 0;
                    if (pCloneDeath > randomNumber())
                        nKill++;

                    if (excitorAmplitude < 0.0) {
                        excitorSpawn[DjIdx] += nKill;
                    } else {
                        excitorSpawn[DjIdx] -= nKill;
                    }
                }

                EProjCluster /= (double) Na;
                EProjCycle += EProjCluster;
                // -- Annihilation --
                // Pairs of particles with opposite signs on the same excitor are removed
                for (uint32_t I = 0; I < excitorWalkers.size(); I++)
                    excitorWalkers[I] += excitorSpawn[I];
            }

            //double EProj = rhfWfn.Energy;
            double EProj = (EProjCycle / NUM_CYCLES)*0.1;// + Vnn;

            // for (uint32_t i = 1; i < excitorWalkers.size(); i++)
            //     EProj += H_FCI(0, i) * abs((double) excitorWalkers[i] / (double) excitorWalkers[0]);

            if (reportCycle >= 50) {
                EProjAverage += EProj;
            }

            if (reportCycle >= 50 && nExcipsLast != 0 && nExcips != 0)
                shift = shift - (1.0 / (tau * 100.0)) * log(abs((double) nExcips / (double) nExcipsLast));

            std::cout << "  Cycle(" << reportCycle << "); E(Proj) = " << EProj <<
                    " Shift = " << shift <<
                    " Population = " <<
                    nExcips
                    << " E(Proj) Average = " << EProjAverage / double(reportCycle + 1 - 5000) << '\n';

            nExcipsLast = nExcips;
        }


        return EProjAverage / double(NUM_REPORT_CYCLES + 1 - 5000);
    }

    std::vector<Excitor> QuantumMonteCarlo::GenerateAllExcitors(CIDeterminant reference, uint32_t nSpinOrbs,
                                                                uint32_t nElectrons) {
        std::vector<Excitor> excitors =
                ConfigurationInteraction::GenerateAllDeterminants(nSpinOrbs, nElectrons);


        using OrbitalBitStr = CIDeterminant::OrbitalBitStr;

        std::vector<uint32_t> uniqueI;
        std::vector<uint32_t> uniqueJ;
        std::vector<uint32_t> commonIJ;
        int sign = 1;

        for (Excitor &ex: excitors) {
            CIDeterminant::CalculateDifferences(reference, ex,
                                                nSpinOrbs / 2,
                                                nSpinOrbs / 2, uniqueI, uniqueJ, commonIJ,
                                                sign);
            ex.GetSign() = sign;

            ex.GetAlphaOrbitals() ^= reference.GetAlphaOrbitals();
            ex.GetBetaOrbitals() ^= reference.GetBetaOrbitals();
        }

        return excitors;
    }
}
