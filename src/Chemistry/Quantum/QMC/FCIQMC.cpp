#include "Application/Application.h"
#include "Chemistry/Quantum/QMC/QMC.h"
#include "Chemistry/Quantum/CI/ConfigurationInteraction.h"
#include "Utilities/Utilities.h"

#include<boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/range/irange.hpp>

namespace SCESS {
    double QuantumMonteCarlo::SimpleFCIQMC(const RHFResult &rhfWavefunction, double Vnn, std::stringstream *console) {
        using namespace Eigen;

        //Evaluate the Full CI  Hamiltonian matrix
        uint32_t nElectrons = std::distance(rhfWavefunction.Occupation.begin(),
                                            std::ranges::find(rhfWavefunction.Occupation, 0)) * 2,
                nSpatialOrbitals = rhfWavefunction.Occupation.size(),
                nSpinOrbs = nSpatialOrbitals * 2;

        std::vector<CIDeterminant> determinants = ConfigurationInteraction::GenerateAllDeterminants(
            nSpinOrbs, nElectrons);

        MatrixXd H_FCI = MatrixXd::Zero((Eigen::Index) determinants.size(),
                                        (Eigen::Index) determinants.size()); {
            MatrixXd HSpin = ConfigurationInteraction::GenerateOneElectronHamiltonianSpin(rhfWavefunction);
            TEITensor moSpinTEI = Integrals::Transforms::AntisymmetriseIntegrateOutSpinClosedShell(
                Integrals::Transforms::ToMolecularTEI(rhfWavefunction.Integrals.tei, rhfWavefunction.Coeffs));

            for (uint32_t i = 0; i < H_FCI.rows(); i++)
                for (uint32_t j = i; j < H_FCI.cols(); j++)
                    H_FCI(i, j) = H_FCI(j, i) =
                                  ConfigurationInteraction::CalculateHamiltonianElement(determinants[i],
                                      determinants[j], HSpin, moSpinTEI);
        }


        SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(H_FCI);
        VectorXd Energies = eigensolver.eigenvalues().real();
        MatrixXd Wavefunctions = eigensolver.eigenvectors().real();

        for (double &energy: Energies)
            energy += Vnn;

        VectorXd WaveFunc0 = Wavefunctions.col(0);


        //Implementation of the QMC algorithm

        constexpr size_t NUM_CYCLES = 200,
                NUM_REPORT_CYCLES = 10;
        constexpr double EPSILON = 0.001;

        using namespace boost::random;
        mt19937 generator(time(0));
        uniform_01<mt19937, double> randomNumber(generator);

        //The key is the determinant the given particle is on
        std::unordered_map<uint32_t, Particle> walkers;
        std::unordered_map<uint32_t, Particle> spawn;

        // Zero is the reference determinant
        // Give it an initial population of 200
        walkers.emplace(0, Particle{200});

        //Reference determinant in the Hartree-Fock determinant which is given by the first index 0
        double runningTotalProjEnergy = 0.0;
        double H00 = H_FCI(0, 0);
        //Imaginary time step
        double tau = 0.0001;
        // The shift factor in the projection energy operator
        double shift = 0.0;
        double nParticles = 0.0;
        double nParticlesLast = 0.0;

        // I: Index of the determinant which we are attempting to spawn on to
        // H_row: Row of the Hamiltonian matrix connected to this determinant
        // tau: Imaginary timestep
        // populationParent: Population of the parent determinant
        auto AttemptSpawn = [&randomNumber, &spawn](uint32_t I,
                                                    int populationParent, double tau, const VectorXd &H_row) {
            for (uint32_t J = 0; J < H_row.size(); J++) {
                //Cannot spawn a determinant onto itself
                if (I == J) continue;

                //  Cannot spawn if the matrix element between this determinant and another one is too small
                //  This indicates that the determinants are disconnected
                if (abs(H_row(J)) < EPSILON) continue;

                //Attempt to spawn particle
                double rate = abs(tau * H_row(J));
                uint32_t nSpawn = static_cast<uint32_t>(rate);
                rate = rate - (double) nSpawn;
                //If the rate is greater than the random number, spawn a new particle
                if (rate > randomNumber()) nSpawn++;

                if (!spawn.contains(J))
                    spawn.emplace(J, Particle{0});

                Particle &JDetParticle = spawn.at(J);

                if (H_row(J) > 0.0) {
                    if (populationParent < 0)
                        JDetParticle.population += nSpawn;
                    else
                        JDetParticle.population -= nSpawn;
                } else {
                    if (populationParent > 0)
                        JDetParticle.population += nSpawn;
                    else
                        JDetParticle.population -= nSpawn;
                }
            }
        };


        auto SimpleDeathClone = [&randomNumber](double tau, double Hii, double H00, double shift, int &population) {
            double rate = abs(population) * tau * (Hii - H00 - shift);
            int numKill = static_cast<int>(rate);
            rate = rate - (double) numKill;

            if (abs(rate) > randomNumber()) {
                if (rate > 0.0)
                    numKill++;
                else
                    numKill--;
            }

            //Particle death takes population closer to zero
            if (population > 0)
                population -= numKill;
            else
                population += numKill;
        };

        auto SimpleAnnihilation = [&spawn, &walkers]() {
            for (auto &[I, walker]: walkers) {
                if (spawn.contains(I))
                    walker.population += spawn.at(I).population;
            }
        };

        auto UpdateShift = [](double &shift, size_t nParticles, size_t nParticlesLast) -> void {
        };

        for (size_t reportCycle = 0; reportCycle < NUM_REPORT_CYCLES; reportCycle++) {
            runningTotalProjEnergy = 0.0;

            for (size_t cycle = 0; cycle < NUM_CYCLES; cycle++) {
                spawn.clear();
                size_t numSpawningAttempts = std::accumulate(walkers.begin(), walkers.end(), (size_t) 0,
                                                             [](size_t &&acc, const std::pair<uint32_t, Particle> &p) {
                                                                 return acc + p.second.population;
                                                             });

                //  Loop over all determinants in the Hilbert space
                //  These will all be considered to be "walkers"
                for (uint32_t I = 0; I < determinants.size(); I++) {
                    double H0i = H_FCI(0, I);
                    double Hii = H_FCI(I, I);

                    if (!walkers.contains(I))
                        walkers.emplace(I, Particle{0});

                    Particle &currentParticle = walkers.at(I);

                    // Only update the projected energy if this particle is not on the
                    // reference determinant
                    if (I != 0) {
                        //Update the total projected energy
                        runningTotalProjEnergy += H0i * currentParticle.population;
                    }

                    //Attempt to spawn more determinants
                    for (uint32_t j = 0; j < currentParticle.population; j++) {
                        AttemptSpawn(I, currentParticle.population, tau, H_FCI.row(I));
                    }

                    SimpleDeathClone(tau, Hii, H00, shift, currentParticle.population);
                }

                SimpleAnnihilation();

                if (console)
                    *console << "       Cycle(" << cycle + 1 << "; E(Projected) = " << runningTotalProjEnergy / double(
                                cycle + 1)
                            << '\n';
            }

            if (console)
                *console << "   Report Cycle(" << reportCycle + 1 << "; E(Projected) = " << runningTotalProjEnergy /
                        NUM_CYCLES
                        << '\n';


            nParticles = std::accumulate(walkers.begin(), walkers.end(), (size_t) 0,
                                         [](size_t &&acc, const std::pair<uint32_t, Particle> &p) {
                                             return acc + p.second.population;
                                         });

            nParticlesLast = nParticles;
        }

        runningTotalProjEnergy /= NUM_CYCLES;

        if (console)
            *console << "Final total projected energy: " << runningTotalProjEnergy << '\n';

        return runningTotalProjEnergy;
    }
} // SCESS
