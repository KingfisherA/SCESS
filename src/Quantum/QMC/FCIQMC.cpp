#include "Quantum/QMC/QMC.h"
#include "Quantum/CI/ConfigurationInteraction.h"
#include "Utilities/Utilities.h"

#include "GraphicsEngine/Utilities/Threading.h"

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/range/irange.hpp>
#include <GraphicsEngine/Utilities/Threading.h>

#include "Quantum/Integrals/IntegralTransforms.h"

namespace SCESS {
    double QuantumMonteCarlo::FCIQMC(const RHFWavefunction &rhfWavefunction, double Vnn,
                                     const GEng::Thread *thread,
                                     std::stringstream *console,
                                     boost::circular_buffer<glm::vec2> *projEnergyGraph) {
        using namespace Eigen;

        if (console) (*console) << "--- Full CI Quantum Monte Carlo ---\n";

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
                Integrals::Transforms::ToMolecularTEI(rhfWavefunction.Integrals.eri, rhfWavefunction.C));

            for (uint32_t i = 0; i < H_FCI.rows(); i++)
                for (uint32_t j = i; j < H_FCI.cols(); j++)
                    H_FCI(i, j) = H_FCI(j, i) =
                                  ConfigurationInteraction::CalculateHamiltonianElement(determinants[i],
                                      determinants[j], nSpatialOrbitals, nSpatialOrbitals, HSpin, moSpinTEI);
        }


        SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(H_FCI);
        VectorXd Energies = eigensolver.eigenvalues().real();
        MatrixXd Wavefunctions = eigensolver.eigenvectors().real();

        for (double &energy: Energies)
            energy += Vnn;

        VectorXd WaveFunc0 = Wavefunctions.col(0);

        if (console) (*console) << " --- FCI Energy: " << Energies(0) << " hartrees---\n\n\n\n";

        //Implementation of the QMC algorithm

        constexpr size_t NUM_CYCLES = 1000;
        constexpr size_t NUM_REPORT_CYCLES = 2000;
        constexpr double EPSILON = 0.01;

        using namespace boost::random;
        mt19937 generator(time(nullptr));
        uniform_01<mt19937, double> randomNumber(generator);

        using Particle = int;
        //The populations of the walker and spawned particles 
        std::vector<Particle> walkerPopulations(determinants.size(), 0);
        std::vector<Particle> spawnPopulations(determinants.size(), 0);

        // Zero is the reference Hartree-Fock determinant
        // Give it an initial population of 100
        walkerPopulations[0] = 100;

        //Reference determinant in the Hartree-Fock determinant which is given by the first index 0
        double runningTotalProjEnergy = 0.0, averageEnergy = 0.0;
        double H00 = H_FCI(0, 0);
        //Imaginary time step
        double tau = 0.0001;
        // The shift factor in the projection energy operator
        double shift = 0.0; //Energies(0) - Vnn;
        int nParticles = 0.0;
        int nParticlesLast = 0.0;

        for (size_t reportCycle = 0; reportCycle < NUM_REPORT_CYCLES; reportCycle++) {
            for (size_t cycle = 0; cycle < NUM_CYCLES; cycle++) {

                if(thread)
                    if(thread->ShouldStop())
                        return 0.0;

                for (Particle &p: spawnPopulations)
                    p = 0;

                //  Loop over all determinants in the Hilbert space
                //  These will all be considered to be "walkers"
                for (uint32_t I = 0; I < determinants.size(); I++) {
                    double Hii = H_FCI(I, I);

                    Particle &currentWalker = walkerPopulations[I];

                    //Attempt to spawn determinants from I onto other determinants
                    for (uint32_t k = 0; k < abs(currentWalker); k++)
                        for (uint32_t J = 0; J < determinants.size(); J++) {
                            if (I == J) continue;

                            if (abs(H_FCI(I, J)) < EPSILON) continue;

                            double pSpawn = abs(tau * H_FCI(I, J));
                            int nSpawn = int(pSpawn);
                            pSpawn = pSpawn - nSpawn;
                            if (pSpawn > randomNumber()) nSpawn++;

                            Particle &JSpawn = spawnPopulations[J];
                            if (H_FCI(I, J) >= 0.0) {
                                //Opposite sign to parent
                                if (currentWalker < 0)
                                    JSpawn += nSpawn;
                                else
                                    JSpawn -= nSpawn;
                            } else {
                                //Same sign as parent
                                if (currentWalker > 0)
                                    JSpawn += nSpawn;
                                else
                                    JSpawn -= nSpawn;
                            }
                        }

                    //Death cloning stage
                    double pCloneDeath = abs(currentWalker) * tau * (Hii - rhfWavefunction.Energy - shift);
                    int nKill = int(pCloneDeath);
                    pCloneDeath = pCloneDeath - nKill;

                    if (abs(pCloneDeath) > randomNumber()) {
                        if (pCloneDeath > 0.0)
                            nKill++;
                        else
                            nKill--;
                    }

                    if (nKill > abs(currentWalker))
                        throw std::runtime_error("Cannot create anti-particles");

                    if (currentWalker > 0)
                        spawnPopulations[I] -= nKill;
                    else
                        spawnPopulations[I] += nKill;
                }

                //Annihilation stage
                for (uint32_t I = 0; I < walkerPopulations.size(); I++)
                    walkerPopulations[I] += spawnPopulations[I];
            }

            double projEnergy = rhfWavefunction.Energy;

            if (walkerPopulations[0] != 0)
                for (uint32_t I = 1; I < walkerPopulations.size(); I++) {
                    projEnergy += H_FCI(0, I) * double(walkerPopulations[I]) / double(walkerPopulations[0]);
                }

            if (projEnergyGraph)
                projEnergyGraph->push_back(glm::vec2{(float) reportCycle, (float) projEnergy});

            if (reportCycle + 1 > NUM_REPORT_CYCLES - 100) {
                runningTotalProjEnergy += projEnergy;
                averageEnergy = runningTotalProjEnergy / double(reportCycle + 1 - (NUM_REPORT_CYCLES - 100));
            }


            if (console)
                (*console) << "  Cycle(" << reportCycle + 1 << "); E(Projected) = " << projEnergy << ", E(average) =" <<
                        averageEnergy
                        << '\n';

            nParticles = std::accumulate(walkerPopulations.begin(), walkerPopulations.end(), 0,
                                         [](int &&acc, const Particle &p)-> int {
                                             return acc + p;
                                         });


            if (reportCycle + 1 > 2 && nParticlesLast != 0)
                shift = shift - (1.0 / (tau * 1000.0)) * log((double) nParticles / (double) nParticlesLast);

            if (console) (*console) << "      Shift: " << shift << '\n';
            if (console) (*console) << "      Number of particles: " << nParticles << '\n';

            nParticlesLast = nParticles;
        }

        if (console) (*console) << "E(Projected) = " << averageEnergy << " hartrees\n";
        if (console) (*console) << "E(FCI) = " << Energies(0) << " hartrees\n";

        return averageEnergy;
    }
} // SCESS
