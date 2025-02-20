#include <Quantum/Q_PCH.h>
#include "F12_QMC.h"
#include "Quantum/F12/F12.h"
//
#include "Tensor/TensorOperations.h"
#include "Tensor/TensorSlice.h"
#include "Tensor/TensorEinsum.h"
#include "Quantum/QMC/QMC.h"
#include "Quantum/CI/CI.h"
//
#include "GraphicsEngine/Tools/Threading.h"
//
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/random/uniform_int.hpp>
//
#include "Quantum/Integrals/IntegralsF12.h"
#include "Quantum/Integrals/IntegralTransforms.h"
#include "Quantum/F12/F12.h"
#include "Utilities/Utilities.h"


//
//
namespace SCESS::F12_QMC {
    //
    //
    using namespace Eigen;
    using namespace boost::random;
    //
    //
#pragma region CCQMC(MP2-F12) Utilities
    //
    //
    std::vector<Excitor> MakeGeminalExcitors(uint32_t nSpinOrbs) {
        return CI::EnumerateDet(nSpinOrbs, 2);
    }
    //
    std::pair<int, int> GetGeminalExcitorElecs(Excitor D_gem, uint32_t nSpatialOrbs) {
        //The index of the electron pair the geminal electrons were put back into
        int x = -1, y = -1;
        for (int i = 0; i < nSpatialOrbs; i++)
            if (D_gem.GetAlphaOrbitals() & (1 << i))
                if (x == -1)
                    x = i;
                else
                    y = i;

        if (y == -1)
            for (int i = 0; i < nSpatialOrbs; i++)
                if (D_gem.GetBetaOrbitals() & (1 << i)) {
                    y = i;
                    break;
                }

        return {x, y};
    }
    //
    double EvaluateGeminalElement(Excitor Di_gem,
                                  SlaterDet Dj,
                                  const TEITensor &V_pqrs,
                                  const TEITensor &C_pqrs,
                                  uint32_t nOcc) {
        size_t numOfAlphaOrbitals = 2, numOfBetaOrbitals = 2, nSpatialOrbs = 2;

        // if (numOfAlphaOrbitals != numOfBetaOrbitals)
        //     throw std::runtime_error(
        //         "Only RHF determinants with equal numbers of alpha and beta number of electrons are supported currently");


        std::vector<uint32_t> uniqueI, uniqueJ, commonIJ;
        uniqueI.reserve(numOfAlphaOrbitals);
        uniqueJ.reserve(numOfAlphaOrbitals);
        commonIJ.reserve(numOfAlphaOrbitals);

        SlaterDet D0 = SlaterDet::MakeReference(1, 1);

        int8_t sign = 1;
        SlaterDet::CalculateDifferences(D0, Dj, numOfAlphaOrbitals, numOfBetaOrbitals,
                                        uniqueI, uniqueJ, commonIJ, sign);


        auto [x,y] = GetGeminalExcitorElecs(Di_gem, nSpatialOrbs);

        if (uniqueI.size() == 1)
            return 0.0;

        if (uniqueI.size() == 0) {
            return V_pqrs(0, 0, x, y) + C_pqrs(0, 0, x, y);
        }

        uint32_t a = uniqueJ[0], b = uniqueJ[1];
        return V_pqrs(a, b, x, y) + C_pqrs(a, b, x, y);
    }
    //
    double EvaluateF12PairEnergy(
        long p, long q,
        const Tensor2 &F_pq,
        const TEITensor &V_pqrs,
        const TEITensor &B_pqrs,
        const TEITensor &X_pqrs) {
        double E_F12_V = (5.0 / 8.0) * V_pqrs(p, q, p, q);
        E_F12_V -= (1.0 / 8.0) * V_pqrs(q, p, p, q);

        double E_F12_X = (F_pq(p, p) + F_pq(q, q)) * (7.0 / 32.0) * X_pqrs(p, q, p, q);
        E_F12_X += (F_pq(p, p) + F_pq(q, q)) * (1.0 / 32.0) * X_pqrs(p, q, q, p);

        double E_F12_B = (7.0 / 32.0) * B_pqrs(p, q, p, q);
        E_F12_B += (1.0 / 32.0) * B_pqrs(p, q, q, p);

        return 2.0 * E_F12_V + E_F12_B - E_F12_X;
    }
    //
    //
#pragma endregion
    //
    //
#pragma region CCQMC(MP2-F12) Implementation
    //
    //
    double CCQMC_MP2_F12(const HF::HF_Wfn &refWfn,
                         const Integrals::MolInts &mints,
                         const Integrals::F12::F12Integrals &f12Ints,
                         F12::F12Idxs f12Idxs) {
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

        using namespace QMC::CC;
        //
        //Number of cycles
        constexpr size_t NUM_CYCLES = 100;
        //Number of attempts
        uint32_t Na = 200;
        constexpr size_t NUM_REPORT_CYCLES = 400;

        std::cout << "--- Coupled Cluster Quantum Monte Carlo F12 ---\n";

        TEITensor spinERI = Integrals::Transforms::AntisymmetriseIntegrateOutSpinClosedShell(
            Slice<Indices<P, Q, R, S> >(idxs,
                                        mints.eri));

        MatrixXd h = mints.T.block(0, 0, 2, 2) +
                     mints.V.block(0, 0, 2, 2);
        MatrixXd HSpin = CI::GenerateHCoreSpin(h);

        //Tensor4 W_ = EvaluateW(f12Ints.f12);

        // Tensor2 F_;
        // TEITensor Vij_, X_, B_;
        // F12::MP2_F12_Fixed(refWfn, mints, f12Ints, f12Idxs, &F_, &Vij_, &X_, &B_);
        // TEITensor Vab_ = F12::EvaluateVab(f12Ints.f12Repulsion, f12Ints.f12, mints.eri, f12Idxs);

        Tensor2 F_pq;
        TEITensor V_pqrs, C_pqrs, B_pqrs, X_pqrs;
        F12::EvaluatePerturbative_CC_F12Intermediates(refWfn, mints, f12Ints, f12Idxs,
                                                      &F_pq, &V_pqrs, &C_pqrs, &B_pqrs, &X_pqrs);

        double val = EvaluateF12PairEnergy(0, 0, F_pq, V_pqrs, B_pqrs, X_pqrs);

        //The number of doubly occupied orbitals
        //Given we are starting with a RHF determinant, the number of alpha=beta which will equal this
        //quantity.
        uint32_t nOccupied = std::distance(refWfn.Occ_Alpha.begin(),
                                           std::ranges::find(refWfn.Occ_Alpha, 0));
        uint32_t nElectrons = nOccupied * 2;
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
        //Generate all possible excitors in this full CC
        //The first "nAnnihilationBits" bits are the annihilators and the
        //remaining "nCreationBits" are the creation operators
        std::vector<Excitor> excitors = EnumerateExcitors(D0, nSpinOrbs, nElectrons);
        //The populations of walkers on each of the excitors generated above
        std::vector<int64_t> excitorWalkers(excitors.size(), 0);
        excitorWalkers[0] = 100;
        std::vector<int64_t> cumulativeWalkers(excitors.size(), 0);
        std::vector<int64_t> excitorSpawn(excitors.size(), 0);

        std::vector<Excitor> geminalExcitors = MakeGeminalExcitors(nSpinOrbs);

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

        //The cumulative 1/2^(s+1) probability distribution
        std::vector<double> cumulativeClusterProb = GenerateCumulativeClusterProb(nElectrons + geminalExcitors.size());

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
                    //Negative indices indicate geminal excitors
                    std::unordered_set<int> selectedExcitorIndices;
                    if (s > 0) {
                        uniform_int_distribution<int>
                                excitorIdxDist(-(geminalExcitors.size() - 1), excitors.size() - 1);
                        selectedExcitorIndices.reserve(s);
                        while (selectedExcitorIndices.size() != s) {
                            int newIndx = excitorIdxDist(randGenerator);
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
                                 pow(1.0 / (excitors.size() + geminalExcitors.size()), s);
                    } else
                        selectedExcitorIndices.insert(0);

                    //Collapse the excitors together taking care of sign
                    Excitor collapsedExcitor;
                    collapsedExcitor.GetSign() = 1.0;
                    double excitorAmplitude = 1.0;
                    bool isGeminal = false;
                    Excitor geminalExcitor;

                    for (uint32_t i = 0; i < nSpatialOrbitals; i++) {
                        geminalExcitor.GetAlphaOrbitals() |= 1 << i;
                        geminalExcitor.GetBetaOrbitals() |= 1 << i;
                    }

                    for (int idx: selectedExcitorIndices) {
                        if (idx < 0) {
                            isGeminal = true;
                            geminalExcitor.GetAlphaOrbitals() =
                                    geminalExcitors[abs(idx)].GetAlphaOrbitals() & geminalExcitor.GetAlphaOrbitals();
                            geminalExcitor.GetBetaOrbitals() =
                                    geminalExcitors[abs(idx)].GetBetaOrbitals() & geminalExcitor.GetBetaOrbitals();
                            continue;
                        }

                        collapsedExcitor.GetAlphaOrbitals() |= excitors[idx].GetAlphaOrbitals();
                        collapsedExcitor.GetBetaOrbitals() |= excitors[idx].GetBetaOrbitals();
                        collapsedExcitor.GetSign() *= excitors[idx].GetSign();

                        excitorAmplitude *= ((double) excitorWalkers[idx] / (double) D0_Population);
                    }

                    // if (geminalExcitor.GetAlphaOrbitals() == 0 && geminalExcitor.GetBetaOrbitals() == 0)
                    //     continue;

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
                    if (isGeminal) {
                        // auto orbIdxs = FindOrbitalIdxs(Di);
                        //Add contributions for each pair of electrons in the determinant

                        auto [x,y] = GetGeminalExcitorElecs(geminalExcitor, nSpatialOrbitals);

                        double H0i = EvaluateF12PairEnergy(0,0, F_pq, V_pqrs, B_pqrs, X_pqrs);
                        EProjCluster += 1.0 * H0i / (pClust*Na);
                    }
                    else if (DiIdx != 0) {
                        double H0i = H_CI(0, DiIdx);
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

                        double Hij = 0.0;
                        //  if (isGeminal)
                        //     Hij = EvaluateGeminalElement(geminalExcitor, Dj, V_pqrs, C_pqrs, nOccupied);
                        // else
                        Hij = H_CI(DiIdx, DjIdx);

                        //1.2 *
                        double pSpawn = tau * abs(excitorAmplitude * Hij) / (Na * pClust * pExcit * pSize * 10.0); //*2.0
                        auto nSpawn = static_cast<int64_t>(pSpawn);
                        pSpawn = pSpawn - (double) nSpawn;
                        if (pSpawn > randomNumber())
                            nSpawn++;

                        nSpawn *= static_cast<int64_t>(Dj.GetSign() * Di.GetSign());
                        if (nSpawn != 0) {
                            if (Hij > 0.0) {
                                excitorSpawn[DjIdx] -= nSpawn;
                            } else {
                                excitorSpawn[DjIdx] += nSpawn;
                            }
                        }
                    }

                    double Hii = H_CI(DiIdx, DiIdx);

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
                    if (nKill != 0) {
                        if (excitorWalkers[DiIdx] > 0) {
                            excitorSpawn[DiIdx] -= nKill;
                        } else {
                            excitorSpawn[DiIdx] += nKill;
                        }
                    }
                }

                EProjCluster /= (double) (Na);

                if (reportCycle < 2)
                    EProjCluster = 0.0;

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

            double EProj = (EProjCycle / NUM_CYCLES) + refWfn.E_HF;

            // double EProj = refWfn.Energy;
            //
            // for (uint32_t I = 1; I < excitorWalkers.size(); I++) {
            //     EProj += H_FCI(0, I) * double(excitorWalkers[I]) / double(excitorWalkers[0]);
            // }


            if (reportCycle >= 50) {
                EProjAverage += EProj;
            }


            if (reportCycle == 50)
                shift = refWfn.E_HF;

            if (reportCycle > 15 && nExcipsLast != 0 && nExcips != 0)
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
    // auto FindOrbitalIdxs = [&nSpatialOrbitals, &nElectrons](SlaterDet det) {
    //     std::vector<uint32_t> orbitalIndices;
    //     orbitalIndices.reserve(nElectrons);
    //
    //     for (int i = 0; i < nSpatialOrbitals; ++i)
    //         if (det.GetAlphaOrbitals() & (1 << i))
    //             orbitalIndices.push_back(i);
    //
    //     for (int i = 0; i < nSpatialOrbitals; ++i)
    //         if (det.GetBetaOrbitals() & (1 << i))
    //             orbitalIndices.push_back(i);
    //
    //     return orbitalIndices;
    // };
    //
    //
}
