//
// Created by charlie on 30/09/24.
//

#include "Chemistry/Quantum/CC/CoupledCluster.h"
#include "Chemistry/Quantum/HF/HartreeFock.h"

#include"Utilities/Utilities.h"
#include<iostream>


namespace SCESS {
    //
    //A direct product decomposition implementation to solve the coupled cluster equations
    //to single and double exciations
    void CoupledCluster::CCSD(const RHFResult &rhfResult) {
        constexpr double CONVERGENCE = 1.0E-6;
        constexpr size_t MAX_ITERATIONS = 20;

        std::cout << "--- Coupled Cluster Singles & Doubles (CCSD) ---\n";

        const ConfigIntegrals &ints = rhfResult.Integrals;
        //Transform H into the MO basis set
        Eigen::MatrixXd H = rhfResult.Coeffs.transpose() * rhfResult.Hcore * rhfResult.Coeffs;

        //Resize H and make block diagonal
        // {
        //     Eigen::Index originalDim = H.rows();
        //     H.conservativeResize(H.rows() * 2, H.rows() * 2);
        //     H.block(originalDim, originalDim, originalDim, originalDim) = H.block(0, 0, originalDim, originalDim);
        // }
        //New Fock matrix in the MO basis set
        Eigen::MatrixXd F = rhfResult.Coeffs.transpose() * rhfResult.Fock * rhfResult.Coeffs;
        // {
        //     Eigen::Index originalDim = F.rows();
        //     F.conservativeResize(F.rows() * 2, F.rows() * 2);
        //     F.block(originalDim, originalDim, originalDim, originalDim) = F.block(0, 0, originalDim, originalDim);
        // }

        Eigen::Index nOcc = std::distance(rhfResult.Occupation.begin(),
                                          std::find(rhfResult.Occupation.begin(), rhfResult.Occupation.end(), 0));
        Eigen::Index nVirt = rhfResult.Occupation.size() - nOcc;

        //Fii
        Eigen::Block<Eigen::MatrixXd> Focc = F.block(0, 0, nOcc, nOcc);
        //Faa
        Eigen::Block<Eigen::MatrixXd> Fvirt = F.block(nOcc, nOcc, nVirt, nVirt);


        Eigen::MatrixXd Dia = Eigen::MatrixXd::Zero(nOcc, nVirt);

        for (Eigen::Index i = 0; i < nOcc; i++)
            for (Eigen::Index j = 0; j < nVirt; j++) {
                Dia(i, j) = Focc(i, i) - Fvirt(j, j);
            }


        Eigen::MatrixXd t1Amps = Eigen::MatrixXd::Zero(nOcc, nVirt);

        TEITensor moTEI = IntegralHelper::CalculateMolecularTEI(rhfResult.Integrals.tei, rhfResult.Coeffs);
        double mp2Energy = HF::MP2ClosedShell(rhfResult, moTEI);

        std::cout << "MP2 Correction: " << mp2Energy << '\n';
        std::cout << "E(Hartree Fock) + E(MP2) = " << rhfResult.Energy + mp2Energy << " hartrees\n";
        moTEI = IntegralHelper::AntisymmetriseTEI(moTEI);

        Eigen::Tensor<double, 4> Dijab(nOcc, nOcc, nVirt, nVirt);
        Eigen::Tensor<double, 4> MOijab(nOcc, nOcc, nVirt, nVirt);

        for (Eigen::Index i = 0; i < nOcc; i++)
            for (Eigen::Index j = 0; j < nOcc; j++)
                for (Eigen::Index k = 0; k < nVirt; k++)
                    for (Eigen::Index l = 0; l < nVirt; l++) {
                        Dijab(i, j, k, l) = Focc(i, i) + Focc(j, j) - Fvirt(k, k) - Fvirt(l, l);
                        MOijab(i, j, k, l) = moTEI(i, j, k, l);
                    }

        Eigen::Tensor<double, 4> t2Amps = MOijab / Dijab;

        auto MakeTildeTau = [&nOcc, &nVirt](const Eigen::MatrixXd &t1, const Eigen::Tensor<double, 4> &t2) {
            Eigen::Tensor<double, 4> tildeTau(t2.dimensions());

            for (Eigen::Index i = 0; i < nOcc; i++)
                for (Eigen::Index j = 0; j < nOcc; j++)
                    for (Eigen::Index a = 0; a < nVirt; a++)
                        for (Eigen::Index b = 0; b < nVirt; b++) {
                            tildeTau(i, j, a, b) = t2(i, j, a, b) + 0.5 * (t1(i, a) * t1(j, b) - t1(i, b) * t1(j, a));
                        }
            return tildeTau;
        };

        auto MakeFae = [ &moTEI, &F, &Fvirt, &nOcc, &nVirt](const Eigen::MatrixXd &t1,
                                                            const Eigen::Tensor<double, 4> &t2,
                                                            const Eigen::Tensor<double, 4> &tildeTau) {
            Eigen::MatrixXd Fae = Eigen::MatrixXd::Zero(nVirt, nVirt);

            for (Eigen::Index a = 0; a < nVirt; a++)
                for (Eigen::Index e = 0; e < nVirt; e++) {
                    double &element = Fae(a, e);

                    if (a != e)
                        element += Fvirt(a, e);

                    for (Eigen::Index m = 0; m < nOcc; m++) {
                        element -= 0.5 * F(m, e + nOcc - 1) * t1(m, a);

                        for (Eigen::Index f = 0; f < nVirt; f++) {
                            element += t1(m, f) * moTEI(m, a + nOcc - 1, f + nOcc - 1, e + nOcc - 1);
                        }

                        for (Eigen::Index n = 0; n < nOcc; n++)
                            for (Eigen::Index f = 0; f < nVirt; f++) {
                                element -= 0.5 * tildeTau(m, n, a, f) * moTEI(m, n, e + nOcc - 1, f + nOcc - 1);
                            }
                    }
                }

            return Fae;
        };

        auto MakeFmi = [ &moTEI, &F, &Focc, &nOcc, &nVirt](const Eigen::MatrixXd &t1,
                                                           const Eigen::Tensor<double, 4> &t2,
                                                           const Eigen::Tensor<double, 4> &tildeTau) {
            Eigen::MatrixXd Fmi = Eigen::MatrixXd::Zero(nOcc, nOcc);

            for (Eigen::Index m = 0; m < nOcc; m++)
                for (Eigen::Index i = 0; i < nOcc; i++) {
                    double &element = Fmi(m, i);

                    if (m != i)
                        element += Focc(m, i);

                    for (Eigen::Index e = 0; e < nVirt; e++) {
                        element += 0.5 * F(m, e + nOcc - 1) * t1(e, i);

                        for (Eigen::Index n = 0; n < nOcc; n++) {
                            element += t1(n, e) * moTEI(m, n, i, e + nOcc - 1);
                        }

                        for (Eigen::Index n = 0; n < nOcc; n++)
                            for (Eigen::Index f = 0; f < nVirt; f++) {
                                element += 0.5 * tildeTau(i, n, e + nOcc - 1, f + nOcc - 1)
                                        * moTEI(m, n, e + nOcc - 1, f + nOcc - 1);
                            }
                    }
                }

            return Fmi;
        };

        auto MakeFme = [&nOcc, &nVirt, &F, &moTEI](const Eigen::MatrixXd &t1, const Eigen::Tensor<double, 4> &t2) {
            Eigen::MatrixXd Fmi = Eigen::MatrixXd::Zero(nOcc, nVirt);

            for (Eigen::Index m = 0; m < nOcc; m++)
                for (Eigen::Index e = 0; e < nVirt; e++) {
                    double &element = Fmi(m, e);

                    element = F(m, e + nOcc - 1);

                    for (Eigen::Index n = 0; n < nOcc; n++)
                        for (Eigen::Index f = 0; f < nVirt; f++)
                            element += t1(n, f) * moTEI(m, n, e + nOcc - 1, f + nOcc - 1);
                }
            return Fmi;
        };


        for (size_t i = 0; i < MAX_ITERATIONS; i++) {
            Eigen::Tensor<double, 4> tildeTau = MakeTildeTau(t1Amps, t2Amps);

            Eigen::MatrixXd
                    Fae = MakeFae(t1Amps, t2Amps, tildeTau),
                    Fmi = MakeFmi(t1Amps, t2Amps, tildeTau),
                    Fme = MakeFme(t1Amps, t2Amps);

            Eigen::Tensor<double, 4>
                    Wmnij,
                    Wabef,
                    Wmbej;
        }
    }
    //
} // SCESS
