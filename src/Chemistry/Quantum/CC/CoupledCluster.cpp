//
// Created by charlie on 30/09/24.
//

#include "Chemistry/Quantum/CC/CoupledCluster.h"
#include "Chemistry/Quantum/HF/HartreeFock.h"

#include"Utilities/Utilities.h"
#include<iostream>
#include <boost/math/constants/constants.hpp>
#include <GL/glew.h>

#include "pch.h"


namespace SCESS {
    template<typename T>
    using MatrixType = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
    //
    template<typename Scalar, typename... Dims>
    auto MatrixToTensor(const MatrixType<Scalar> &matrix, Dims... dims) {
        constexpr int rank = sizeof...(Dims);

        Eigen::Tensor<double, rank> tensor(dims...);

        memcpy(tensor.data(), matrix.data(), matrix.size() * sizeof(Scalar));

        return tensor;
    }
    //
    std::optional<CCSDResult> CoupledCluster::CCSD(const RHFResult &rhfResult, std::stringstream *console) {
        constexpr double CC_CONVERGENCE = 1.0E-6;
        constexpr size_t MAX_ITERATIONS = 50;

        TEITensor moTEI = Integrals::Transforms::ToMolecularTEI(rhfResult.Integrals.tei, rhfResult.Coeffs);
        double mp2Energy = HartreeFock::MP2ClosedShell(rhfResult, moTEI, console);

        if (console) {
            *console << std::setprecision(12);
            *console << "--- Coupled Cluster Singles & Doubles (CCSD) ---\n";
        }

        //Antisymmetrise the two-electron integrals and integrate out the alpha and beta
        //spin functions for a closed-shell wavefunction
        moTEI = Integrals::Transforms::AntisymmetriseIntegrateOutSpinClosedShell(moTEI);

        //Number of occupied spin orbitals
        Eigen::Index nOcc = std::distance(rhfResult.Occupation.begin(),
                                          std::ranges::find(rhfResult.Occupation, 0));
        //For two electrons in each orbital
        nOcc *= 2;
        //Number of virtual spin orbitals
        Eigen::Index nVirt = rhfResult.Occupation.size() * 2 - nOcc;

        Eigen::Tensor<double, 2> F = [&rhfResult]() {
            //Transform the Fock matrix into the MO basis set and reshape so that the alpha and beta
            //spin electrons are in their own block diagonal
            Eigen::MatrixXd FockMat = rhfResult.Coeffs.transpose() * rhfResult.Fock * rhfResult.Coeffs;
            Eigen::Tensor<double, 2> FockTensor(FockMat.rows() * 2,
                                                FockMat.rows() * 2);
            FockTensor.setZero();
            for (Eigen::Index i = 0; i < FockTensor.dimension(0); i++)
                FockTensor(i, i) = FockMat(i / 2, i / 2);


            return FockTensor;
        }();

        //Slices of the Fock matricies for occupied and virtual orbitals
        Eigen::Tensor<double, 2> Focc = F.slice(Eigen::array<Eigen::Index, 2>{0, 0},
                                                Eigen::array<Eigen::Index, 2>{nOcc, nOcc});
        Eigen::Tensor<double, 2> Fvirt = F.slice(Eigen::array<Eigen::Index, 2>{nOcc, nOcc},
                                                 Eigen::array<Eigen::Index, 2>{nVirt, nVirt});

        //Denominator arrays
        Eigen::Tensor<double, 2> Dia(nOcc, nVirt);
        Eigen::Tensor<double, 4> Dijab(nOcc, nOcc, nVirt, nVirt);

        for (Eigen::Index i = 0; i < nOcc; i++)
            for (Eigen::Index a = 0; a < nVirt; a++)
                Dia(i, a) = Focc(i, i) - Fvirt(a, a);

        for (Eigen::Index i = 0; i < nOcc; i++)
            for (Eigen::Index j = 0; j < nOcc; j++)
                for (Eigen::Index a = 0; a < nVirt; a++)
                    for (Eigen::Index b = 0; b < nVirt; b++)
                        Dijab(i, j, a, b) = Focc(i, i) + Focc(j, j) - Fvirt(a, a) - Fvirt(b, b);


        //The first excitation amplitudes tensor
        //Dimensions of o,v
        Eigen::Tensor<double, 2> t1Amps(nOcc, nVirt);
        t1Amps.setZero();

        //The second excitation amplitudes tensor
        //Dimensions of o,o,v,v
        Eigen::Tensor<double, 4> t2Amps; {
            Eigen::Tensor<double, 4> MOijab = moTEI.slice(Eigen::array<Eigen::Index, 4>{0, 0, nOcc, nOcc},
                                                          Eigen::array<Eigen::Index, 4>{nOcc, nOcc, nVirt, nVirt});
            t2Amps = MOijab / Dijab;
        }


        auto MakeTildeTau = [&nOcc, &nVirt](const Eigen::Tensor<double, 2> &t1, const Eigen::Tensor<double, 4> &t2) {
            Eigen::Tensor<double, 4> tildeTau = t2;

            //ia,jb-> iajb -> ijab
            Eigen::Tensor<double, 4> tmp = 0.5 * t1.contract(t1, Eigen::array<Eigen::IndexPair<int>, 0>{}).shuffle(
                                               Eigen::array<int, 4>{0, 2, 1, 3});
            tildeTau += tmp;
            tildeTau -= tmp.shuffle(Eigen::array<int, 4>{0, 1, 3, 2});

            return tildeTau;
        };

        auto MakeTau = [&nOcc, &nVirt](const Eigen::Tensor<double, 2> &t1,
                                       const Eigen::Tensor<double, 4> &t2) {
            Eigen::Tensor<double, 4> tau = t2;

            //ia,jb-> iajb -> ijab
            Eigen::Tensor<double, 4> tmp = t1.contract(t1, Eigen::array<Eigen::IndexPair<int>, 0>{}).shuffle(
                Eigen::array<int, 4>{0, 2, 1, 3});
            tau += tmp;
            tau -= tmp.shuffle(Eigen::array<int, 4>{0, 1, 3, 2});


            return tau;
        };

        auto MakeFae = [ &moTEI, &F, &Fvirt, &nOcc, &nVirt](const Eigen::Tensor<double, 2> &t1,
                                                            const Eigen::Tensor<double, 4> &t2,
                                                            const Eigen::Tensor<double, 4> &tildeTau) {
            Eigen::Tensor<double, 2> Fae(nVirt, nVirt);
            Fae.setZero();

            //me,ma-> ea -> ae
            Fae -= 0.5 * F.slice(Eigen::array<Eigen::Index, 2>{0, nOcc},
                                 Eigen::array<Eigen::Index, 2>{nOcc, nVirt})
                    .contract(
                        t1, Eigen::array<Eigen::IndexPair<int>, 1>{
                            Eigen::IndexPair<int>{0, 0}
                        }).shuffle(Eigen::array<int, 2>{1, 0});

            //MO[o,v,v,v]
            //mf,mafe -> ae
            Fae += t1.contract(moTEI.slice(Eigen::array<Eigen::Index, 4>{0, nOcc, nOcc, nOcc},
                                           Eigen::array<Eigen::Index, 4>{nOcc, nVirt, nVirt, nVirt}),
                               Eigen::array<Eigen::IndexPair<int>, 2>{
                                   Eigen::IndexPair<int>{0, 0},
                                   Eigen::IndexPair<int>{1, 2}
                               });

            //MO[o, o, v, v]
            //mnaf,mnef->ae
            Fae -= 0.5 * tildeTau.contract(moTEI.slice(Eigen::array<Eigen::Index, 4>{0, 0, nOcc, nOcc},
                                                       Eigen::array<Eigen::Index, 4>{nOcc, nOcc, nVirt, nVirt}),
                                           Eigen::array<Eigen::IndexPair<int>, 3>{
                                               Eigen::IndexPair<int>{0, 0},
                                               Eigen::IndexPair<int>{1, 1},
                                               Eigen::IndexPair<int>{3, 3}
                                           });

            return Fae;
        };

        auto MakeFmi = [ &moTEI, &F, &Focc, &nOcc, &nVirt](const Eigen::Tensor<double, 2> &t1,
                                                           const Eigen::Tensor<double, 4> &t2,
                                                           const Eigen::Tensor<double, 4> &tildeTau) {
            Eigen::Tensor<double, 2> Fmi(nOcc, nOcc);
            Fmi.setZero();

            //F[o, v]
            //ie,me-> im -> mi
            Fmi += 0.5 * F.slice(Eigen::array<Eigen::Index, 2>{0, nOcc},
                                 Eigen::array<Eigen::Index, 2>{nOcc, nVirt})
                    .contract(
                        t1, Eigen::array<Eigen::IndexPair<int>, 1>{
                            Eigen::IndexPair<int>{1, 1}
                        }).shuffle(Eigen::array<int, 2>{1, 0});

            //MO[o, o, o, v]
            //ne,mnie->mi
            Fmi += t1.contract(moTEI.slice(Eigen::array<Eigen::Index, 4>{0, 0, 0, nOcc},
                                           Eigen::array<Eigen::Index, 4>{nOcc, nOcc, nOcc, nVirt}),
                               Eigen::array<Eigen::IndexPair<int>, 2>{
                                   Eigen::IndexPair<int>{0, 1},
                                   Eigen::IndexPair<int>{1, 3}
                               });

            //MO[o, o, v, v]
            //inef,mnef-> im -> mi
            Fmi += 0.5 * tildeTau.contract(moTEI.slice(Eigen::array<Eigen::Index, 4>{0, 0, nOcc, nOcc},
                                                       Eigen::array<Eigen::Index, 4>{nOcc, nOcc, nVirt, nVirt}),
                                           Eigen::array<Eigen::IndexPair<int>, 3>{
                                               Eigen::IndexPair<int>{1, 1},
                                               Eigen::IndexPair<int>{2, 2},
                                               Eigen::IndexPair<int>{3, 3}
                                           }).shuffle(Eigen::array<int, 2>{1, 0});

            return Fmi;
        };

        auto MakeFme = [&nOcc, &nVirt, &F, &moTEI](const Eigen::Tensor<double, 2> &t1,
                                                   const Eigen::Tensor<double, 4> &t2) {
            Eigen::Tensor<double, 2> Fme = F.slice(Eigen::array<Eigen::Index, 2>{0, nOcc},
                                                   Eigen::array<Eigen::Index, 2>{nOcc, nVirt});
            //nf,mnef->me
            Fme += t1.contract(moTEI.slice(Eigen::array<Eigen::Index, 4>{0, 0, nOcc, nOcc},
                                           Eigen::array<Eigen::Index, 4>{nOcc, nOcc, nVirt, nVirt}),
                               Eigen::array<Eigen::IndexPair<int>, 2>{
                                   Eigen::IndexPair<int>{0, 1},
                                   Eigen::IndexPair<int>{1, 3}
                               });

            return Fme;
        };

        auto MakeWmnij = [&moTEI, &nOcc, &nVirt, &MakeTau](const Eigen::Tensor<double, 2> &t1,
                                                           const Eigen::Tensor<double, 4> &t2,
                                                           const Eigen::Tensor<double, 4> &tau) {
            Eigen::Tensor<double, 4> Wmnij = moTEI.slice(Eigen::array<Eigen::Index, 4>{0, 0, 0, 0},
                                                         Eigen::array<Eigen::Index, 4>{nOcc, nOcc, nOcc, nOcc});

            //je,mnie-> jmni -> mnij
            Eigen::Tensor<double, 4> slicedMO = moTEI.slice(Eigen::array<Eigen::Index, 4>{0, 0, 0, nOcc},
                                                            Eigen::array<Eigen::Index, 4>{
                                                                nOcc, nOcc, nOcc, nVirt
                                                            });
            Eigen::Tensor<double, 4> Pij = t1.contract(slicedMO,
                                                       Eigen::array<Eigen::IndexPair<int>, 1>{
                                                           Eigen::IndexPair<int>{1, 3}
                                                       }).shuffle(Eigen::array<int, 4>{1, 2, 3, 0});
            Wmnij += Pij;
            Wmnij -= Pij.shuffle(Eigen::array<int, 4>{0, 1, 3, 2});

            //ijef,mnef-> ijmn -> mnij
            Wmnij += 0.25 * tau.contract(moTEI.slice(Eigen::array<Eigen::Index, 4>{0, 0, nOcc, nOcc},
                                                     Eigen::array<Eigen::Index, 4>{nOcc, nOcc, nVirt, nVirt}),
                                         Eigen::array<Eigen::IndexPair<int>, 2>{
                                             Eigen::IndexPair<int>{2, 2},
                                             Eigen::IndexPair<int>{3, 3}
                                         }).shuffle(Eigen::array<int, 4>{2, 3, 0, 1});


            return Wmnij;
        };

        auto MakeWabef = [&moTEI, &nVirt, &nOcc](const Eigen::Tensor<double, 2> &t1,
                                                 const Eigen::Tensor<double, 4> &t2,
                                                 const Eigen::Tensor<double, 4> &tau) {
            Eigen::Tensor<double, 4> Wabef = moTEI.slice(Eigen::array<Eigen::Index, 4>{nOcc, nOcc, nOcc, nOcc},
                                                         Eigen::array<Eigen::Index, 4>{nVirt, nVirt, nVirt, nVirt});

            //mb,amef-> baef -> abef
            Eigen::Tensor<double, 4> Pab = t1.contract(moTEI.slice(Eigen::array<Eigen::Index, 4>{nOcc, 0, nOcc, nOcc},
                                                                   Eigen::array<Eigen::Index, 4>{
                                                                       nVirt, nOcc, nVirt, nVirt
                                                                   }),
                                                       Eigen::array<Eigen::IndexPair<int>, 1>{
                                                           Eigen::IndexPair<int>{0, 1}
                                                       })
                    .shuffle(Eigen::array<int, 4>{1, 0, 2, 3});

            Wabef -= Pab;
            Wabef += Pab.shuffle(Eigen::array<int, 4>{1, 0, 2, 3});

            //MO[o, o, v, v]
            //mnab,mnef->abef
            Wabef += 0.25 * tau.contract(moTEI.slice(Eigen::array<Eigen::Index, 4>{0, 0, nOcc, nOcc},
                                                     Eigen::array<Eigen::Index, 4>{nOcc, nOcc, nVirt, nVirt}),
                                         Eigen::array<Eigen::IndexPair<int>, 2>{
                                             Eigen::IndexPair<int>{0, 0},
                                             Eigen::IndexPair<int>{1, 1}
                                         });


            return Wabef;
        };


        auto MakeWmbej = [&nOcc, &nVirt, &moTEI](const Eigen::Tensor<double, 2> &t1,
                                                 const Eigen::Tensor<double, 4> &t2) {
            Eigen::Tensor<double, 4> Wmbej = moTEI.slice(Eigen::array<Eigen::Index, 4>{0, nOcc, nOcc, 0},
                                                         Eigen::array<Eigen::Index, 4>{nOcc, nVirt, nVirt, nOcc});

            //jf,mbef-> jmbe -> mbej
            Wmbej += t1.contract(moTEI.slice(Eigen::array<Eigen::Index, 4>{0, nOcc, nOcc, nOcc},
                                             Eigen::array<Eigen::Index, 4>{nOcc, nVirt, nVirt, nVirt}),
                                 Eigen::array<Eigen::IndexPair<int>, 1>{
                                     Eigen::IndexPair<int>{1, 3}
                                 })
                    .shuffle(Eigen::array<int, 4>{1, 2, 3, 0});

            //nb,mnej-> bmej -> mbej
            Wmbej -= t1.contract(moTEI.slice(Eigen::array<Eigen::Index, 4>{0, 0, nOcc, 0},
                                             Eigen::array<Eigen::Index, 4>{nOcc, nOcc, nVirt, nOcc}),
                                 Eigen::array<Eigen::IndexPair<int>, 1>{
                                     Eigen::IndexPair<int>{0, 1}
                                 })
                    .shuffle(Eigen::array<int, 4>{1, 0, 2, 3});

            //jf,nb-> jfnb -> jnfb
            Eigen::Tensor<double, 4> tmp = (0.5 * t2) + t1.contract(t1, Eigen::array<Eigen::IndexPair<int>, 0>{})
                                           .shuffle(Eigen::array<int, 4>{0, 2, 1, 3});

            //MO[o, o, v, v]
            //jnfb,mnef-> jbme -> mbej
            Wmbej -= 0.25 * tmp.contract(moTEI.slice(Eigen::array<Eigen::Index, 4>{0, 0, nOcc, nOcc},
                                                     Eigen::array<Eigen::Index, 4>{nOcc, nOcc, nVirt, nVirt}),
                                         Eigen::array<Eigen::IndexPair<int>, 2>{
                                             Eigen::IndexPair<int>{1, 1},
                                             Eigen::IndexPair<int>{2, 3}
                                         })
                    .shuffle(Eigen::array<int, 4>{2, 1, 3, 0});


            return Wmbej;
        };

        double lastCorrelationEnergy = -DBL_MIN, correlationEnergy = 0.0;
        size_t ccIteration = 0;

        for (; ccIteration < MAX_ITERATIONS; ccIteration++) {
            Eigen::Tensor<double, 4>
                    tildeTau = MakeTildeTau(t1Amps, t2Amps),
                    tau = MakeTau(t1Amps, t2Amps);

            Eigen::Tensor<double, 2>
                    Fae = MakeFae(t1Amps, t2Amps, tildeTau),
                    Fmi = MakeFmi(t1Amps, t2Amps, tildeTau),
                    Fme = MakeFme(t1Amps, t2Amps);

            Eigen::Tensor<double, 4>
                    Wmnij = MakeWmnij(t1Amps, t2Amps, tau),
                    Wabef = MakeWabef(t1Amps, t2Amps, tau),
                    Wmbej = MakeWmbej(t1Amps, t2Amps);


            //Make the t1 RHS
            //Dimensions of nOcc and nVirt
            Eigen::Tensor<double, 2> t1RHS = F.slice(Eigen::array<Eigen::Index, 2>{0, nOcc},
                                                     Eigen::array<Eigen::Index, 2>{nOcc, nVirt});

            //ie,ae -> ia
            t1RHS += t1Amps.contract(Fae, Eigen::array<Eigen::IndexPair<int>, 1>{
                                         Eigen::IndexPair<int>{1, 1}
                                     });
            //ma,mi -> ai -> ia
            t1RHS -= t1Amps.contract(Fmi, Eigen::array<Eigen::IndexPair<int>, 1>{
                                         Eigen::IndexPair<int>{0, 0}
                                     }).shuffle(Eigen::array<int, 2>{1, 0});

            //imae,me -> ia
            t1RHS += t2Amps.contract(Fme, Eigen::array<Eigen::IndexPair<int>, 2>{
                                         Eigen::IndexPair<int>{1, 0},
                                         Eigen::IndexPair<int>{3, 1}
                                     });

            //nf,naif -> ai -> ia
            t1RHS -= t1Amps.contract(moTEI
                                     .slice(Eigen::array<Eigen::Index, 4>{0, nOcc, 0, nOcc},
                                            Eigen::array<Eigen::Index, 4>{nOcc, nVirt, nOcc, nVirt}),
                                     Eigen::array<Eigen::IndexPair<int>, 2>{
                                         Eigen::IndexPair<int>{0, 0},
                                         Eigen::IndexPair<int>{1, 3}
                                     })
                    .shuffle(Eigen::array<int, 2>{1, 0});


            //imef, maef -> ia
            t1RHS -= 0.5 * t2Amps.contract(moTEI
                                           .slice(Eigen::array<Eigen::Index, 4>{0, nOcc, nOcc, nOcc},
                                                  Eigen::array<Eigen::Index, 4>{nOcc, nVirt, nVirt, nVirt}),
                                           Eigen::array<Eigen::IndexPair<int>, 3>{
                                               Eigen::IndexPair<int>{1, 0},
                                               Eigen::IndexPair<int>{2, 2},
                                               Eigen::IndexPair<int>{3, 3}
                                           });

            //mnae, nmei -> ai -> ia
            t1RHS -= 0.5 * t2Amps.contract(moTEI
                                           .slice(Eigen::array<Eigen::Index, 4>{0, 0, nOcc, 0},
                                                  Eigen::array<Eigen::Index, 4>{nOcc, nOcc, nVirt, nOcc}),
                                           Eigen::array<Eigen::IndexPair<int>, 3>{
                                               Eigen::IndexPair<int>{0, 1},
                                               Eigen::IndexPair<int>{1, 0},
                                               Eigen::IndexPair<int>{3, 2}
                                           })
                    .shuffle(Eigen::array<int, 2>{1, 0});

            //Make the t2 RHS
            Eigen::Tensor<double, 4> t2RHS = moTEI.slice(Eigen::array<Eigen::Index, 4>{0, 0, nOcc, nOcc},
                                                         Eigen::array<Eigen::Index, 4>{nOcc, nOcc, nVirt, nVirt});
            //Fae - 1/2 Sum_m(tmb * Fme)
            //mb,me->be
            Eigen::Tensor<double, 2> tmp = Fae - 0.5 * Eigen::Tensor<double, 2>(t1Amps.contract(Fme,
                                               Eigen::array<Eigen::IndexPair<int>, 1>{
                                                   Eigen::IndexPair<int>{0, 0}
                                               }));

            //Do ijae,be->ijab = 0123,01->0123
            Eigen::Tensor<double, 4> Pab = t2Amps.contract(tmp, Eigen::array<Eigen::IndexPair<int>, 1>{
                                                               Eigen::IndexPair<int>{3, 1}
                                                           });
            t2RHS += Pab;
            t2RHS -= Pab.shuffle(Eigen::array<int, 4>{0, 1, 3, 2});

            //me,je->mj
            tmp = Fmi + 0.5 * Fme.contract(t1Amps,
                                           Eigen::array<Eigen::IndexPair<int>, 1>{
                                               Eigen::IndexPair<int>{1, 1}
                                           });
            //Do imab,mj->iabj->ijab
            Eigen::Tensor<double, 4> Pij = t2Amps.contract(tmp, Eigen::array<Eigen::IndexPair<int>, 1>{
                                                               Eigen::IndexPair<int>{1, 0}
                                                           }).shuffle(Eigen::array<int, 4>{0, 3, 1, 2});
            t2RHS -= Pij;
            t2RHS += Pij.shuffle(Eigen::array<int, 4>{1, 0, 2, 3});

            //mnij,mnab->ijab
            t2RHS += 0.5 * Wmnij.contract(tau, Eigen::array<Eigen::IndexPair<int>, 2>{
                                              Eigen::IndexPair<int>{0, 0},
                                              Eigen::IndexPair<int>{1, 1}
                                          });
            //ijef,abef->ijab
            t2RHS += 0.5 * tau.contract(Wabef, Eigen::array<Eigen::IndexPair<int>, 2>{
                                            Eigen::IndexPair<int>{2, 2},
                                            Eigen::IndexPair<int>{3, 3}
                                        });

            //mbej,ie -> mbji,ma->bjia -> ijab
            Eigen::Tensor<double, 4> tmp2 = Eigen::Tensor<double, 4>(moTEI.slice(
                Eigen::array<Eigen::Index, 4>{0, nOcc, nOcc, 0},
                Eigen::array<Eigen::Index, 4>{
                    nOcc, nVirt, nVirt, nOcc
                })).contract(t1Amps, Eigen::array<Eigen::IndexPair<int>, 1>{
                                 Eigen::IndexPair<int>{2, 1}
                             }).contract(
                t1Amps, Eigen::array<Eigen::IndexPair<int>, 1>{
                    Eigen::IndexPair<int>{0, 0}
                }).shuffle(Eigen::array<int, 4>{2, 1, 3, 0});


            //imae,mbej->iabj->ijab
            Eigen::Tensor<double, 4> Pijab = t2Amps.contract(Wmbej, Eigen::array<Eigen::IndexPair<int>, 2>{
                                                                 Eigen::IndexPair<int>{1, 0},
                                                                 Eigen::IndexPair<int>{3, 2}
                                                             }).shuffle(Eigen::array<int, 4>{0, 3, 1, 2});
            Pijab -= tmp2;
            //
            t2RHS += Pijab;
            t2RHS -= Pijab.shuffle(Eigen::array<int, 4>{0, 1, 3, 2});
            t2RHS -= Pijab.shuffle(Eigen::array<int, 4>{1, 0, 2, 3});
            t2RHS += Pijab.shuffle(Eigen::array<int, 4>{1, 0, 3, 2});

            //abej,ie -> abji -> ijab
            Pij = Eigen::Tensor<double, 4>(moTEI.slice(Eigen::array<Eigen::Index, 4>{nOcc, nOcc, nOcc, 0},
                                                       Eigen::array<Eigen::Index, 4>{
                                                           nVirt, nVirt, nVirt, nOcc
                                                       })).contract(t1Amps, Eigen::array<Eigen::IndexPair<int>, 1>{
                                                                        Eigen::IndexPair<int>{2, 1}
                                                                    }).shuffle(Eigen::array<int, 4>{3, 2, 0, 1});


            t2RHS += Pij;
            t2RHS -= Pij.shuffle(Eigen::array<int, 4>{1, 0, 2, 3});

            //mbij,ma -> bija -> ijab
            Pab = Eigen::Tensor<double, 4>(moTEI.slice(Eigen::array<Eigen::Index, 4>{0, nOcc, 0, 0},
                                                       Eigen::array<Eigen::Index, 4>{
                                                           nOcc, nVirt, nOcc, nOcc
                                                       })).contract(t1Amps, Eigen::array<Eigen::IndexPair<int>, 1>{
                                                                        Eigen::IndexPair<int>{0, 0}
                                                                    }).shuffle(Eigen::array<int, 4>{1, 2, 3, 0});

            t2RHS -= Pab;
            t2RHS += Pab.shuffle(Eigen::array<int, 4>{0, 1, 3, 2});

            t1Amps = t1RHS / Dia;
            t2Amps = t2RHS / Dijab;

            Eigen::Tensor<double, 0> correlationEnergyTensor =
                    F.slice(Eigen::array<Eigen::Index, 2>{0, nOcc},
                            Eigen::array<Eigen::Index, 2>{
                                nOcc, nVirt
                            }).contract(
                        t1Amps, Eigen::array<Eigen::IndexPair<int>, 2>{
                            Eigen::IndexPair<int>{0, 0},
                            Eigen::IndexPair<int>{1, 1}
                        });
            Eigen::Tensor<double, 4> moSlice = moTEI.slice(Eigen::array<Eigen::Index, 4>{0, 0, nOcc, nOcc},
                                                           Eigen::array<Eigen::Index, 4>{
                                                               nOcc, nOcc, nVirt, nVirt
                                                           });
            correlationEnergyTensor += 0.25 * moSlice.contract(t2Amps,
                                                               Eigen::array<Eigen::IndexPair<int>, 4>{
                                                                   Eigen::IndexPair<int>{0, 0},
                                                                   Eigen::IndexPair<int>{1, 1},
                                                                   Eigen::IndexPair<int>{2, 2},
                                                                   Eigen::IndexPair<int>{3, 3}
                                                               });
            correlationEnergyTensor += 0.5 * moSlice.contract(t1Amps,
                                                              Eigen::array<Eigen::IndexPair<int>, 2>{
                                                                  Eigen::IndexPair<int>{0, 0},
                                                                  Eigen::IndexPair<int>{2, 1},
                                                              }).contract(t1Amps,
                                                                          Eigen::array<Eigen::IndexPair<int>, 2>{
                                                                              Eigen::IndexPair<int>{0, 0},
                                                                              Eigen::IndexPair<int>{1, 1},
                                                                          });

            double dE = abs(correlationEnergyTensor(0) - lastCorrelationEnergy);
            if (console)
                *console << "   CCSD(" << ccIteration << "): E(corr) = " << correlationEnergyTensor(0) << "      dE = " << dE
                        << '\n';

            if (dE < CC_CONVERGENCE) {
                correlationEnergy = correlationEnergyTensor(0);
                if (console)
                    *console << "---- CCSD procedure converged! ----\n";
                break;
            }

            lastCorrelationEnergy = correlationEnergyTensor(0);
        }

        if (ccIteration == MAX_ITERATIONS - 1) {
            if (console)
                *console << "!!! CCSD procedure failed to converge !!!\n\n\n";
            return {};
        }

        if (console) {
            *console << "E(CCSD correlation) = " << correlationEnergy << " hartrees\n";
            *console << "E(Total) = E(Hartree Fock) + E(CCSD) = " << rhfResult.Energy + correlationEnergy <<
                    " hartrees\n\n\n";
        }
        return CCSDResult{mp2Energy, correlationEnergy, t1Amps, t2Amps};
    }
    //
    //     auto MakeTildeTau = [&nOcc, &nVirt](const Eigen::Tensor<double, 2> &t1, const Eigen::Tensor<double, 4> &t2) {
    //     Eigen::Tensor<double, 4> tildeTau(t2.dimensions());
    //
    //     for (Eigen::Index i = 0; i < nOcc; i++)
    //         for (Eigen::Index j = 0; j < nOcc; j++)
    //             for (Eigen::Index a = 0; a < nVirt; a++)
    //                 for (Eigen::Index b = 0; b < nVirt; b++) {
    //                     tildeTau(i, j, a, b) = t2(i, j, a, b) + 0.5 * (t1(i, a) * t1(j, b) - t1(i, b) * t1(j, a));
    //                 }
    //     return tildeTau;
    // };
    //
    // auto MakeTau = [&nOcc, &nVirt](const Eigen::Tensor<double, 2> &t1,
    //                                const Eigen::Tensor<double, 4> &t2) {
    //     Eigen::Tensor<double, 4> tau(t2.dimensions());
    //
    //     for (Eigen::Index i = 0; i < nOcc; i++)
    //         for (Eigen::Index j = 0; j < nOcc; j++)
    //             for (Eigen::Index a = 0; a < nVirt; a++)
    //                 for (Eigen::Index b = 0; b < nVirt; b++)
    //                     tau(i, j, a, b) = t2(i, j, a, b) + t1(i, a) * t1(j, b) - t1(i, b) * t1(j, a);
    //
    //     return tau;
    // };
    //
    // auto MakeFae = [ &moTEI, &F, &Fvirt, &nOcc, &nVirt](const Eigen::Tensor<double, 2> &t1,
    //                                                     const Eigen::Tensor<double, 4> &t2,
    //                                                     const Eigen::Tensor<double, 4> &tildeTau) {
    //     Eigen::Tensor<double, 2> Fae(nVirt, nVirt);
    //
    //     for (Eigen::Index a = 0; a < nVirt; a++)
    //         for (Eigen::Index e = 0; e < nVirt; e++) {
    //             double &element = Fae(a, e);
    //
    //             if (a != e)
    //                 element += Fvirt(a, e);
    //
    //             for (Eigen::Index m = 0; m < nOcc; m++) {
    //                 element -= 0.5 * F(m, e + nOcc - 1) * t1(m, a);
    //
    //                 for (Eigen::Index f = 0; f < nVirt; f++) {
    //                     element += t1(m, f) * moTEI(m, a + nOcc - 1, f + nOcc - 1, e + nOcc - 1);
    //                 }
    //
    //                 for (Eigen::Index n = 0; n < nOcc; n++)
    //                     for (Eigen::Index f = 0; f < nVirt; f++) {
    //                         element -= 0.5 * tildeTau(m, n, a, f) * moTEI(m, n, e + nOcc - 1, f + nOcc - 1);
    //                     }
    //             }
    //         }
    //
    //     return Fae;
    // };
    //
    // auto MakeFmi = [ &moTEI, &F, &Focc, &nOcc, &nVirt](const Eigen::Tensor<double, 2> &t1,
    //                                                    const Eigen::Tensor<double, 4> &t2,
    //                                                    const Eigen::Tensor<double, 4> &tildeTau) {
    //     Eigen::Tensor<double, 2> Fmi(nOcc, nOcc);
    //
    //     for (Eigen::Index m = 0; m < nOcc; m++)
    //         for (Eigen::Index i = 0; i < nOcc; i++) {
    //             double &element = Fmi(m, i);
    //
    //             if (m != i)
    //                 element += Focc(m, i);
    //
    //             for (Eigen::Index e = 0; e < nVirt; e++) {
    //                 element += 0.5 * F(m, e + nOcc - 1) * t1(e, i);
    //
    //                 for (Eigen::Index n = 0; n < nOcc; n++) {
    //                     element += t1(n, e) * moTEI(m, n, i, e + nOcc - 1);
    //                 }
    //
    //                 for (Eigen::Index n = 0; n < nOcc; n++)
    //                     for (Eigen::Index f = 0; f < nVirt; f++) {
    //                         element += 0.5 * tildeTau(i, n, e + nOcc - 1, f + nOcc - 1)
    //                                 * moTEI(m, n, e + nOcc - 1, f + nOcc - 1);
    //                     }
    //             }
    //         }
    //
    //     return Fmi;
    // };
    //
    // auto MakeFme = [&nOcc, &nVirt, &F, &moTEI](const Eigen::Tensor<double, 2> &t1,
    //                                            const Eigen::Tensor<double, 4> &t2) {
    //     Eigen::Tensor<double, 2> Fmi(nOcc, nVirt);
    //
    //     for (Eigen::Index m = 0; m < nOcc; m++)
    //         for (Eigen::Index e = 0; e < nVirt; e++) {
    //             double &element = Fmi(m, e);
    //
    //             element = F(m, e + nOcc - 1);
    //
    //             for (Eigen::Index n = 0; n < nOcc; n++)
    //                 for (Eigen::Index f = 0; f < nVirt; f++)
    //                     element += t1(n, f) * moTEI(m, n, e + nOcc - 1, f + nOcc - 1);
    //         }
    //     return Fmi;
    // };
    //
    // auto MakeWmnij = [&moTEI, &nOcc, &nVirt, &MakeTau](const Eigen::Tensor<double, 2> &t1,
    //                                                    const Eigen::Tensor<double, 4> &t2,
    //                                                    const Eigen::Tensor<double, 4> &tau) {
    //     Eigen::Tensor<double, 4> Wmnij(nOcc, nOcc, nOcc, nOcc);
    //     Eigen::Tensor<double, 4> Pij(Wmnij.dimensions());
    //
    //     for (Eigen::Index m = 0; m < nOcc; m++)
    //         for (Eigen::Index n = 0; n < nOcc; n++)
    //             for (Eigen::Index i = 0; i < nOcc; i++)
    //                 for (Eigen::Index j = 0; j < nOcc; j++)
    //                     for (Eigen::Index e = 0; e < nVirt; e++) {
    //                         Pij(m, n, i, j) += t1(j, e) * moTEI(m, n, i, e + nOcc - 1);
    //                     }
    //
    //     for (Eigen::Index m = 0; m < nOcc; m++)
    //         for (Eigen::Index n = 0; n < nOcc; n++)
    //             for (Eigen::Index i = 0; i < nOcc; i++)
    //                 for (Eigen::Index j = 0; j < nOcc; j++) {
    //                     double &element = Wmnij(m, n, i, j);
    //                     element = moTEI(m, n, i, j);
    //                     element += Pij(m, n, i, j);
    //                     element -= Pij(m, n, j, i);
    //
    //                     for (Eigen::Index e = 0; e < nVirt; e++)
    //                         for (Eigen::Index f = 0; f < nVirt; f++)
    //                             element += 0.25 * tau(i, j, e, f) * moTEI(m, n, e + nOcc - 1, f + nOcc - 1);
    //                 }
    //
    //     return Wmnij;
    // };
    //
    // auto MakeWabef = [&moTEI, &nVirt, &nOcc](const Eigen::Tensor<double, 2> &t1,
    //                                          const Eigen::Tensor<double, 4> &t2,
    //                                          const Eigen::Tensor<double, 4> &tau) {
    //     Eigen::Tensor<double, 4> Wabef(nVirt, nVirt, nVirt, nVirt);
    //     Eigen::Tensor<double, 4> Pab(Wabef.dimensions());
    //
    //     for (Eigen::Index a = 0; a < nVirt; a++)
    //         for (Eigen::Index b = 0; b < nVirt; b++)
    //             for (Eigen::Index e = 0; e < nVirt; e++)
    //                 for (Eigen::Index f = 0; f < nVirt; f++)
    //                     for (Eigen::Index m = 0; m < nOcc; m++) {
    //                         Pab(a, b, e, f) += t1(m, b) * moTEI(a + nOcc - 1, m, e + nOcc - 1, f + nOcc - 1);
    //                     }
    //
    //     for (Eigen::Index a = 0; a < nVirt; a++)
    //         for (Eigen::Index b = 0; b < nVirt; b++)
    //             for (Eigen::Index e = 0; e < nVirt; e++)
    //                 for (Eigen::Index f = 0; f < nVirt; f++) {
    //                     double &element = Wabef(a, b, e, f);
    //                     element = moTEI(a + nOcc - 1, b + nOcc - 1, e + nOcc - 1, f + nOcc - 1);
    //                     element -= Pab(a, b, e, f);
    //                     element += Pab(b, a, e, f);
    //
    //                     for (Eigen::Index m = 0; m < nOcc; m++)
    //                         for (Eigen::Index n = 0; n < nOcc; n++)
    //                             element += 0.25 * tau(m, n, a, b) * moTEI(m, n, e + nOcc - 1, f + nOcc - 1);
    //                 }
    //
    //     return Wabef;
    // };
    //
    //
    // auto MakeWmbej = [&nOcc, &nVirt, &moTEI](const Eigen::Tensor<double, 2> &t1,
    //                                          const Eigen::Tensor<double, 4> &t2) {
    //     Eigen::Tensor<double, 4> Wmbej(nOcc, nVirt, nVirt, nOcc);
    //
    //     for (Eigen::Index m = 0; m < nOcc; m++)
    //         for (Eigen::Index b = 0; b < nVirt; b++)
    //             for (Eigen::Index e = 0; e < nVirt; e++)
    //                 for (Eigen::Index j = 0; j < nOcc; j++) {
    //                     double &element = Wmbej(m, b, e, j);
    //                     element = moTEI(m, b + nOcc - 1, e + nOcc - 1, j);
    //
    //                     for (Eigen::Index f = 0; f < nVirt; f++)
    //                         element += t1(j, f) * moTEI(m, b + nOcc - 1, e + nOcc - 1, f + nOcc - 1);
    //
    //                     for (Eigen::Index n = 0; n < nOcc; n++)
    //                         element -= t1(n, b) * moTEI(m, n, e + nOcc - 1, j);
    //
    //                     for (Eigen::Index n = 0; n < nOcc; n++)
    //                         for (Eigen::Index f = 0; f < nVirt; f++)
    //                             element -= (0.5 * t2(j, n, f, b) + t1(j, f) * t1(n, b)) * moTEI(
    //                                 m, n, e + nOcc - 1, f + nOcc - 1);
    //                 }
    //
    //     return Wmbej;
    // };

    //
} // SCESS
