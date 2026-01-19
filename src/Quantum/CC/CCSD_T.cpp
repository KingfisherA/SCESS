#include <Quantum/Q_PCH.h>
//
#include <Eigen/Core>
#include <unsupported/Eigen/CXX11/Tensor>
//
#include "Dahlia/Tools/Threading.h"
#include "Quantum/CC/CC.h"

//
//
namespace SCESS::CC {
    using namespace Eigen;
    //
    //
    CCSD_Wfn CCSD_T(const HFWfn &rhfResult,
                    const Integrals::MolInts &mints,
                    const DH::Thread *thread,
                    std::stringstream *console) {
        if (console)
            *console << "--- Coupled Cluster Singles & Doubles with perturbative Triples (CCSD(T)) ---\n";

        //Perform CCSD with DIIS acceleration and get t1 and t2 amplitudes
        CCSD_Wfn ccsdWfn = CCSD_DIIS(rhfResult, mints, thread, console);

        // if (ccsdWfn == nullptr) {
        //     if (console)
        //         *console << "Cannot calculate triples correction as the CCSD t amplitudes did not converge.\n";
        //
        //     return nullptr;
        // }

        if (console)
            *console << "Calculating triples correction to the energy ... \n";

        //Number of occupied spin orbitals
        Index nOcc = std::distance(rhfResult.Occ_Alpha.begin(),
                                   std::ranges::find(rhfResult.Occ_Alpha, 0));
        //For two electrons in each orbital
        nOcc *= 2;
        //Number of virtual spin orbitals
        Index nVirt = rhfResult.Occ_Alpha.size() * 2 - nOcc;

        //Get the Fock matrix in canonical MO basis set
        FockTensor F = TransformFockToMO(rhfResult);

        if (thread)
            if (thread->ShouldStop())
                return {};

        //Calculate the triples perturbation correction energy
        double E_T = CalculateTriplesCorrection(ccsdWfn, F, mints.eri, nOcc, nVirt, thread);

        if (console) {
            *console << "E(T) = " << E_T << "\n";
            *console << "E(Total) = E(Hartree Fock) + E(CCSD) + E((T)) = " <<
                rhfResult.E_HF + ccsdWfn.E_CCSD + E_T << '\n';
        }

        ccsdWfn.E_T = E_T;

        return ccsdWfn;
    }
    //
    //Computes the triple correction energy with the spin integrated out
    double CalculateTriplesCorrection(const CCSD_Wfn &ccsd,
                                      const FockTensor &F,
                                      const TEITensor &eri,
                                      Index nOcc,
                                      Index nVirt,
                                      const DH::Thread *thread) {
        using IdxP = Eigen::IndexPair<int>;


        //Construct W; the connected triples tensor
        auto MakeWijkabc = [&eri, &ccsd, nOcc, nVirt]() {
            //jkae,eibc -> jkaibc -> ijkabc
            Tensor<double, 6> Wijkabc = ccsd.t2.contract(eri.slice(array<Eigen::Index, 4>{nOcc, 0, nOcc, nOcc},
                                                                   array<Eigen::Index, 4>{
                                                                       nVirt, nOcc, nVirt, nVirt
                                                                   }),
                                                         array<IdxP, 1>{
                                                             IdxP{3, 0}
                                                         }).shuffle(array<int, 6>{3, 0, 1, 2, 4, 5});

            //imbc,majk-> ibcajk -> ijkabc
            Wijkabc -= ccsd.t2.contract(eri.slice(array<Eigen::Index, 4>{0, nOcc, 0, 0},
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
        auto MakeVijkabc = [&ccsd, &eri, nOcc, nVirt]() {
            const T1Tensor &t1 = ccsd.t1;

            //ia,jkbc -> iajkbc -> ijkabc
            Tensor<double, 6> Vijkabc = t1.contract(eri.slice(array<Eigen::Index, 4>{0, 0, nOcc, nOcc},
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

            if (thread)
                if (thread->ShouldStop())
                    return 0.0;

            Tensor<double, 6> sumTensor = (Wijkabc + MakeVijkabc()) / MakeDijkabc();

            if (thread)
                if (thread->ShouldStop())
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
    //
    //
}
