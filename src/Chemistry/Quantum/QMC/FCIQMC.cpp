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

        //The key is the determinant the given particle is on
        std::unordered_map<uint32_t, Particle> particles;
        std::vector<int> spawn;

        constexpr size_t NUM_CYCLES = 20;
        constexpr double EPSILON = 0.001;

        using namespace boost::random;
        mt19937 generator;
        uniform_01<mt19937, double> randomNumber(generator);

        //Testing
        for (auto i: boost::irange(0, 10))
            std::cout << randomNumber() << '\n';


        // I: Index of the determinant which we are attempting to spawn on to
        // H_row: Row of the Hamiltonian matrix connected to this determinant
        // tau: Imaginary timestep
        auto AttemptSpawn = [&randomNumber](uint32_t I, double tau, VectorXd &H_row) {
            for (uint32_t J = 0; J < H_row.size(); J++) {
                //Cannot spawn a determinant onto itself
                if (I == J) continue;

                //  Cannot spawn if the matrix element between this determinant and another one is too small
                //  This indicates that the determinants are disconnected
                if (abs(H_row(J)) < EPSILON) continue;

                //Attempt to spawn particle
                randomNumber();
                double rate = abs(tau * H_row(J));
                uint32_t nSpawn = static_cast<uint32_t>(rate);
                rate = rate - (double) nSpawn;
                if (rate > randomNumber()) nSpawn++;
            }
        };

        //Reference determinant in the Hartree-Fock determinant which is given by the first index 0

        size_t D0_population = 0;
        double runningTotalProjEnergy = 0.0;

        for (size_t cycle = 0; cycle < NUM_CYCLES; cycle++) {
            spawn.clear();
            double cycleProjEnergy = 0.0;
            size_t numSpawningAttempts = particles.size();

            //Loop over all determinants in the Hilbert space
            for (size_t I = 0; I < determinants.size(); I++) {
                double H0i = H_FCI(0, I);
                double Hii = H_FCI(I, I);

                Particle &currentParticle = particles.at(I);

                // Only update the projected energy if this particle is not on the
                // reference determinant
                if (I != 0) {
                    //Update the the projected energy
                    runningTotalProjEnergy += H0i * currentParticle.population;
                } else {
                    D0_population += currentParticle.population;
                }

                //Attempt to spawn more determinants
                for (uint32_t j = 0; j < currentParticle.population; j++) {
                }
            }
        }

        runningTotalProjEnergy /= NUM_CYCLES;


        return runningTotalProjEnergy;
    }
} // SCESS
