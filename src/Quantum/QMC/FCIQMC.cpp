#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/range/irange.hpp>
#include <boost/math/special_functions/binomial.hpp>
#include <unordered_map>
#include  <Dahlia/Tools/Threading.h>
//
#include "Quantum/QMC/QMC.h"
#include "Quantum/CI/CI.h"
#include "Utilities/Utilities.h"
#include "Quantum/Integrals/Integrals.h"
#include "Quantum/Integrals/IntegralTransforms.h"
#include "Quantum/HF/HF.h"
#include "Dahlia/Tools/Threading.h"


//
//
namespace SCESS::QMC::FCI {
    using namespace Eigen;
    //
    //
    double FCIQMC(const HF_Wfn &rhfWfn,
                  const Integrals::MolInts &mints,
                  double Vnn,
                  const DH::Thread *thread,
                  std::stringstream *console,
                  const ProjEnergyCallback &onReportCycle,
                  std::optional<uint64_t> seed) {
        using CI::SlaterDetHash;
        if (console) (*console) << "--- Full CI Quantum Monte Carlo ---\n";

        const uint32_t nOccupied = std::distance(rhfWfn.Occ_Alpha.begin(), std::ranges::find(rhfWfn.Occ_Alpha, 0)),
                nSpatialOrbitals = rhfWfn.Occ_Alpha.size();

        //Spatial MO integrals; the spin-orbital elements are formed on demand
        const CI::SpinOrbitalIntegrals ints = CI::MakeSpinOrbitalIntegrals(rhfWfn, mints);
        auto Element = [&](SlaterDet I, SlaterDet J) {
            return CI::CalculateHamiltonianElement(I, J, nSpatialOrbitals, nSpatialOrbitals, ints);
        };

        //The exact energy for comparison when the M_s = 0 space is small enough to diagonalise quickly
        constexpr double EXACT_FCI_LIMIT = 1e5;
        const double nDeterminants = std::pow(boost::math::binomial_coefficient<double>(nSpatialOrbitals, nOccupied),
                                              2.0);
        std::optional<double> E_FCI;
        if (nDeterminants <= EXACT_FCI_LIMIT) {
            E_FCI = CI::CI(rhfWfn, mints, Vnn, -1, thread).E_CI(0);
            if (console) (*console) << " --- FCI Energy: " << *E_FCI << " hartrees---\n\n\n\n";
        }
        if (thread && thread->ShouldStop())
            return 0.0;

        //Implementation of the QMC algorithm

        constexpr size_t NUM_CYCLES = 1000;
        constexpr size_t NUM_REPORT_CYCLES = 2000;

        using namespace boost::random;
        mt19937 generator(seed ? (uint32_t) *seed : (uint32_t) time(nullptr));
        uniform_01<double> unitDist;
        auto randomNumber = [&generator, &unitDist]() { return unitDist(generator); };

        //The occupied determinants with their walker populations and diagonal elements
        struct Walkers {
            int64_t population = 0;
            double Hii = 0.0;
        };
        std::unordered_map<SlaterDet, Walkers, SlaterDetHash> walkers;
        std::unordered_map<SlaterDet, int64_t, SlaterDetHash> spawned;

        //The reference Hartree-Fock determinant starts with a population of 100
        const SlaterDet D0 = SlaterDet::MakeReference(nOccupied, nOccupied);
        //Electronic energy of the reference; the shift is a correlation energy relative to it
        const double H00 = Element(D0, D0);
        walkers[D0] = {100, H00};

        double runningTotalProjEnergy = 0.0, averageEnergy = 0.0;
        //Imaginary time step
        double tau = 0.0001;
        // The shift factor in the projection energy operator
        double shift = 0.0;
        int64_t nParticles = 0;
        int64_t nParticlesLast = 0;

        for (size_t reportCycle = 0; reportCycle < NUM_REPORT_CYCLES; reportCycle++) {
            for (size_t cycle = 0; cycle < NUM_CYCLES; cycle++) {
                if (thread)
                    if (thread->ShouldStop())
                        return 0.0;

                spawned.clear();

                for (const auto &[Di, w]: walkers) {
                    if (w.population == 0)
                        continue;
                    const int parentSign = w.population > 0 ? 1 : -1;

                    //-- Spawning --
                    //Every walker attempts one spawn onto a random single or double excitation Dj of Di with the
                    //probability tau |H_ij| / p_gen
                    for (int64_t k = 0; k < std::abs(w.population); k++) {
                        double pGen = 0.5;
                        const SlaterDet Dj = randomNumber() < 0.5
                                                 ? CC::RandomSingleExcite(Di, generator, nSpatialOrbitals, pGen)
                                                 : CC::RandomDoubleExcite(Di, generator, nSpatialOrbitals, pGen);
                        if (Dj == Di)
                            continue;
                        const double Hij = Element(Dj, Di);
                        if (Hij == 0.0)
                            continue;

                        double pSpawn = tau * std::abs(Hij) / pGen;
                        auto nSpawn = static_cast<int64_t>(pSpawn);
                        pSpawn -= (double) nSpawn;
                        if (pSpawn > randomNumber())
                            nSpawn++;
                        if (nSpawn == 0)
                            continue;

                        //Opposite sign to the parent for a positive element, the same sign otherwise
                        spawned[Dj] += (Hij >= 0.0 ? -parentSign : parentSign) * nSpawn;
                    }

                    //-- Death and cloning --
                    double pCloneDeath = std::abs(w.population) * tau * (w.Hii - H00 - shift);
                    auto nKill = static_cast<int64_t>(pCloneDeath);
                    pCloneDeath -= (double) nKill;
                    if (std::abs(pCloneDeath) > randomNumber())
                        nKill += pCloneDeath > 0.0 ? 1 : -1;

                    if (nKill > std::abs(w.population))
                        throw std::runtime_error("Cannot create anti-particles");

                    spawned[Di] -= parentSign * nKill;
                }

                //-- Annihilation --
                for (const auto &[Dj, n]: spawned) {
                    if (n == 0)
                        continue;
                    auto it = walkers.find(Dj);
                    if (it == walkers.end())
                        walkers.emplace(Dj, Walkers{n, Element(Dj, Dj)});
                    else
                        it->second.population += n;
                }
                std::erase_if(walkers, [&D0](const auto &entry) {
                    return entry.second.population == 0 && !(entry.first == D0);
                });
            }

            const int64_t N0 = walkers.at(D0).population;
            double projEnergy = rhfWfn.E_HF;
            if (N0 != 0)
                for (const auto &[Dj, w]: walkers)
                    if (!(Dj == D0))
                        projEnergy += Element(D0, Dj) * double(w.population) / double(N0);

            if (onReportCycle)
                onReportCycle(reportCycle, projEnergy);

            if (reportCycle + 1 > NUM_REPORT_CYCLES - 100) {
                runningTotalProjEnergy += projEnergy;
                averageEnergy = runningTotalProjEnergy / double(reportCycle + 1 - (NUM_REPORT_CYCLES - 100));
            }

            if (console)
                (*console) << "  Cycle(" << reportCycle + 1 << "); E(Projected) = " << projEnergy << ", E(average) =" <<
                        averageEnergy << '\n';

            //Total number of walkers of either sign
            nParticles = 0;
            for (const auto &[Dj, w]: walkers)
                nParticles += std::abs(w.population);

            if (reportCycle + 1 > 2 && nParticlesLast != 0 && nParticles != 0)
                shift = shift - (1.0 / (tau * 1000.0)) * log((double) nParticles / (double) nParticlesLast);

            if (console) (*console) << "      Shift: " << shift << '\n';
            if (console) (*console) << "      Number of particles: " << nParticles << '\n';
            if (console) (*console) << "      Occupied determinants: " << walkers.size() << '\n';

            nParticlesLast = nParticles;
        }

        if (console) (*console) << "E(Projected) = " << averageEnergy << " hartrees\n";
        if (console && E_FCI) (*console) << "E(FCI) = " << *E_FCI << " hartrees\n";

        return averageEnergy;
    }
    //
    //
} // SCESS
