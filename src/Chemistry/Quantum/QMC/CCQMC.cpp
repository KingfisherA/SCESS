#include "Chemistry/Quantum/QMC/QMC.h"
#include "Chemistry/Quantum/CI/ConfigurationInteraction.h"

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_01.hpp>

// Linked coupled-cluster; <Di |e^{-T} H e^{T}| D0>
// Unlinked coupled-cluster; <Di | H e^{T}| D0>

namespace SCESS {
    using namespace Eigen;
    using namespace boost::random;

    double QuantumMonteCarlo::CCQMC(const RHFWavefunction &rhfWfn) {
        constexpr size_t NUM_CYCLES = 100;
        constexpr size_t NUM_REPORT_CYCLES = 100;

        constexpr uint32_t MAX_EXCITATION_LEVEL = 2;

        TEITensor moTEI = Integrals::Transforms::AntisymmetriseIntegrateOutSpinClosedShell(
            Integrals::Transforms::ToMolecularTEI(rhfWfn.Integrals.tei, rhfWfn.Coeffs));

        MatrixXd HSpin = ConfigurationInteraction::GenerateOneElectronHamiltonianSpin(rhfWfn);


        //nElectrons = Number of spin orbitals occupied
        uint32_t nElectrons = std::distance(rhfWfn.Occupation.begin(),
                                            std::ranges::find(rhfWfn.Occupation, 0)) * 2,
                nSpatialOrbitals = rhfWfn.Occupation.size(),
                nSpinOrbs = nSpatialOrbitals * 2;

        //Generate all the excitors in this full CC
        std::vector<Excitor> excitors;
        excitors.reserve(static_cast<size_t>(Integrals::Binomial(nSpinOrbs, nElectrons)));

        for (uint32_t i = 0; i < nElectrons; i++) {
        }


        mt19937 randGenerator(time(nullptr));
        uniform_01<mt19937> randomNumber(randGenerator);

        double shift = 0.0;

        // Excitors; These are operators which take the reference determinant onto another
        // determinant such as |Di> = ai^\dag |D0> where ai is the excitor


        for (size_t reportCycle = 0; reportCycle < NUM_REPORT_CYCLES; reportCycle++) {
            for (size_t cycle = 0; cycle < NUM_CYCLES; cycle++) {
                uint32_t maxClusterSize = std::min<uint32_t>(
                    std::min<uint32_t>(nElectrons, MAX_EXCITATION_LEVEL + 2),
                    0);

                //Update the projected energy by sampling the clusters randomly


                // -- Cluster selection --
                // Randomly choose a number of excitors to sample between 0 and the truncation level (l) + 2
                // These excitators will be combined with each other and used to generate a given
                // excited determinant |Di>
                // This step is equivalent to truncating the exp(T) Taylor expansion at the l+2th power
                // The 1/n! factor from the Taylor expansion is combined into the weight of the cluster
                int numExcitors = 0;


                // -- Spawning --
                // A random single or double excitation of |Di> is performed to generate |Dj>
                // We now spawn a particle onto |Dj> with probability proportional to Hij which corresponds to
                // <Di | H | Dj>. The sign of the particle is equal to Hij


                // -- Death --
                // A new particle is spawned back onto |Di> with the sign of |Hii-S| with sign opposite
                // to the parent


                // -- Annihilation --
                // Pairs of particles with opposite signs on the same excitor are removed
            }
        }


        return 0.0;
    }
}
