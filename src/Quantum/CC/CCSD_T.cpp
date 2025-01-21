#include <Quantum/QuantumPCH.h>

#include <Eigen/Core>
#include <unsupported/Eigen/CXX11/Tensor>

#include "GraphicsEngine/Utilities/Threading.h"

#include "Quantum/CC/CoupledCluster.h"

namespace SCESS {
    using namespace Eigen;

    uRef<CCSD_TWavefunction> CoupledCluster::CCSD_T(const RHFWavefunction &rhfResult, const GEng::Thread *thread,
                                                    std::stringstream *console) {
        if (console)
            *console << "--- Coupled Cluster Singles & Doubles with perturbative Triples (CCSD(T)) ---\n";

        //Perform CCSD with DIIS acceleration and get t1 and t2 amplitudes
        uRef<CCSDWavefunction> ccsd = CoupledCluster::CCSD_DIIS(rhfResult, thread, console);

        if (ccsd == nullptr) {
            if (console)
                *console << "Cannot calculate triples correction as the CCSD t amplitudes did not converge.\n";

            return nullptr;
        }

        if (console)
            *console << "Calculating triples correction to the energy ... \n";

        //Number of occupied spin orbitals
        Eigen::Index nOcc = std::distance(rhfResult.Occupation.begin(),
                                          std::ranges::find(rhfResult.Occupation, 0));
        //For two electrons in each orbital
        nOcc *= 2;
        //Number of virtual spin orbitals
        Eigen::Index nVirt = rhfResult.Occupation.size() * 2 - nOcc;

        //Get the Fock matrix in canonical MO basis set
        FockTensor F = TransformFockToMO(rhfResult);

        if (thread)
            if (thread->ShouldStop())
                return nullptr;

        //Calculate the triples perturbation correction energy
        double triplesCorrection = CalculateTriplesCorrection(*ccsd, F, nOcc, nVirt, thread);

        if (console) {
            *console << "E(T) = " << triplesCorrection << "\n";
            *console << "E(Total) = E(Hartree Fock) + E(CCSD) + E((T)) = " << rhfResult.Energy + ccsd->
                    ccsdCorrelationEnergy + triplesCorrection << '\n';
        }

        return std::make_unique<CCSD_TWavefunction>(*ccsd, triplesCorrection);
    }

    //Computes the triple correction energy with the spin integrated out
    double CoupledCluster::CalculateTriplesCorrection(const CCSDWavefunction &ccsd,
                                                      const FockTensor &F, const Eigen::Index nOcc,
                                                      const Eigen::Index nVirt,
                                                      const GEng::Thread *thread) {
        using IdxP = Eigen::IndexPair<int>;
        const TEITensor &moTEI = ccsd.moTEI;


        //Construct W; the connected triples tensor
        auto MakeWijkabc = [&moTEI, &ccsd, nOcc, nVirt]() {
            //jkae,eibc -> jkaibc -> ijkabc
            Tensor<double, 6> Wijkabc = ccsd.t2Amps.contract(moTEI.slice(array<Eigen::Index, 4>{nOcc, 0, nOcc, nOcc},
                                                                         array<Eigen::Index, 4>{
                                                                             nVirt, nOcc, nVirt, nVirt
                                                                         }),
                                                             array<IdxP, 1>{
                                                                 IdxP{3, 0}
                                                             }).shuffle(array<int, 6>{3, 0, 1, 2, 4, 5});

            //imbc,majk-> ibcajk -> ijkabc
            Wijkabc -= ccsd.t2Amps.contract(moTEI.slice(array<Eigen::Index, 4>{0, nOcc, 0, 0},
                                                        array<Eigen::Index, 4>{nOcc, nVirt, nOcc, nOcc}),
                                            array<IdxP, 1>{
                                                IdxP{1, 0}
                                            }).shuffle(array<int, 6>{0, 4, 5, 3, 1, 2});


            //P(i/jk) * P(a/bc)
            // (ijk - jik - kji) * (abc - bac - cba)
            //(ijkabc - ijkbac - ijkcba - jikabc + jikbac + jikcba - kjiabc + kjibac + kjicba)
            //Perform the Pijkabc permutation operation
            Wijkabc =
                    //ijkabc
                    Wijkabc.shuffle(array<int, 6>{0, 1, 2, 3, 4, 5})
                    //- ijkbac
                    - Wijkabc.shuffle(array<int, 6>{0, 1, 2, 4, 3, 5})
                    //- ijkcba
                    - Wijkabc.shuffle(array<int, 6>{0, 1, 2, 5, 4, 3})
                    //- jikabc
                    - Wijkabc.shuffle(array<int, 6>{1, 0, 2, 3, 4, 5})
                    //jikbac
                    + Wijkabc.shuffle(array<int, 6>{1, 0, 2, 4, 3, 5})
                    //jikcba
                    + Wijkabc.shuffle(array<int, 6>{1, 0, 2, 5, 4, 3})
                    //- kjiabc
                    - Wijkabc.shuffle(array<int, 6>{2, 1, 0, 3, 4, 5})
                    //+ kjibac
                    + Wijkabc.shuffle(array<int, 6>{2, 1, 0, 4, 3, 5})
                    //+ kjicba
                    + Wijkabc.shuffle(array<int, 6>{2, 1, 0, 5, 4, 3});

            return Wijkabc;
        };

        //Construct the V tensor; the disconnected triples tensor
        auto MakeVijkabc = [&ccsd, &moTEI, nOcc, nVirt]() {
            const T1Tensor &t1 = ccsd.t1Amps;

            //ia,jkbc -> iajkbc -> ijkabc
            Tensor<double, 6> Vijkabc = t1.contract(moTEI.slice(array<Eigen::Index, 4>{0, 0, nOcc, nOcc},
                                                                array<Eigen::Index, 4>{nOcc, nOcc, nVirt, nVirt}),
                                                    array<IdxP, 0>{}).shuffle(array<int, 6>{0, 2, 3, 1, 4, 5});

            //P(i/jk) * P(a/bc)
            // (ijk - jik - kji) * (abc - bac - cba)
            //(ijkabc - ijkbac - ijkcba - jikabc + jikbac + jikcba - kjiabc + kjibac + kjicba)
            //Perform the Pijkabc permutation operation
            Vijkabc =
                    //ijkabc
                    Vijkabc.shuffle(array<int, 6>{0, 1, 2, 3, 4, 5})
                    //- ijkbac
                    - Vijkabc.shuffle(array<int, 6>{0, 1, 2, 4, 3, 5})
                    //- ijkcba
                    - Vijkabc.shuffle(array<int, 6>{0, 1, 2, 5, 4, 3})
                    //- jikabc
                    - Vijkabc.shuffle(array<int, 6>{1, 0, 2, 3, 4, 5})
                    //jikbac
                    + Vijkabc.shuffle(array<int, 6>{1, 0, 2, 4, 3, 5})
                    //jikcba
                    + Vijkabc.shuffle(array<int, 6>{1, 0, 2, 5, 4, 3})
                    //- kjiabc
                    - Vijkabc.shuffle(array<int, 6>{2, 1, 0, 3, 4, 5})
                    //- kjibac
                    + Vijkabc.shuffle(array<int, 6>{2, 1, 0, 4, 3, 5})
                    //- kjicba
                    + Vijkabc.shuffle(array<int, 6>{2, 1, 0, 5, 4, 3});

            return Vijkabc;
        };

        //Construct D; The denominator tensor composed of the diagonal elements of the
        //Fock matrix
        auto MakeDijkabc = [&F, nOcc, nVirt]() {
            Tensor<double, 6> Dijkabc(nOcc, nOcc, nOcc, nVirt, nVirt, nVirt);

            //Slices of the Fock matricies for occupied and virtual orbitals
            Tensor<double, 2> Focc = F.slice(array<Eigen::Index, 2>{0, 0},
                                             array<Eigen::Index, 2>{nOcc, nOcc});
            Tensor<double, 2> Fvirt = F.slice(array<Eigen::Index, 2>{nOcc, nOcc},
                                              array<Eigen::Index, 2>{nVirt, nVirt});

            for (Eigen::Index i = 0; i < nOcc; i++)
                for (Eigen::Index j = 0; j < nOcc; j++)
                    for (Eigen::Index k = 0; k < nOcc; k++)
                        for (Eigen::Index a = 0; a < nVirt; a++)
                            for (Eigen::Index b = 0; b < nVirt; b++)
                                for (Eigen::Index c = 0; c < nVirt; c++) {
                                    Dijkabc(i, j, k, a, b, c) =
                                            Focc(i, i)
                                            + Focc(j, j)
                                            + Focc(k, k)
                                            - Fvirt(a, a)
                                            - Fvirt(b, b)
                                            - Fvirt(c, c);
                                }

            return Dijkabc;
        };

        //Evaluate the triples perturbed correction energy
        double tCorrectionEnergy = 0.0; {
            Tensor<double, 6> Wijkabc = MakeWijkabc();

            if(thread)
                if(thread->ShouldStop())
                    return 0.0;

            Tensor<double, 6> sumTensor = (Wijkabc + MakeVijkabc()) / MakeDijkabc();

            if(thread)
                if(thread->ShouldStop())
                    return 0.0;

            Tensor<double, 0> correctionETensor = Wijkabc.contract(sumTensor, array<IdxP, 6>{
                                                                       IdxP{0, 0},
                                                                       IdxP{1, 1},
                                                                       IdxP{2, 2},
                                                                       IdxP{3, 3},
                                                                       IdxP{4, 4},
                                                                       IdxP{5, 5}
                                                                   });
            tCorrectionEnergy = correctionETensor(0);
        }

        return (1.0 / 36.0) * tCorrectionEnergy;
    }
}
