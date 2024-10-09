#include<optional>
#include <Eigen/Core>
#include <unsupported/Eigen/CXX11/Tensor>

#include "Chemistry/Quantum/CC/CoupledCluster.h"

namespace SCESS {
    std::optional<CCSDResult> CoupledCluster::CCSD_DIIS(const RHFResult &rhfResult, std::stringstream *console) {
        using namespace Eigen;

        constexpr double CC_CONVERGENCE = 1.0E-6;
        constexpr size_t MAX_ITERATIONS = 50;
        constexpr bool ENABLE_DIIS = true;
        constexpr size_t MAX_DIIS_VECS = 6;

        TEITensor moTEI = Integrals::Transforms::ToMolecularTEI(rhfResult.Integrals.tei, rhfResult.Coeffs);
        double mp2Energy = HartreeFock::MP2ClosedShell(rhfResult, moTEI, console);

        if (console) {
            *console << std::setprecision(12);
            *console << "--- Coupled Cluster Singles & Doubles (CCSD) with DIIS acceleration ---\n";
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

        FockTensor F = TransformFockToMO(rhfResult);

        //Denominator arrays
        Tensor<double, 2> Dia(nOcc, nVirt);
        Tensor<double, 4> Dijab(nOcc, nOcc, nVirt, nVirt);
        CreateDenominatorTensors(F, nOcc, nVirt, Dia, Dijab);

        //The first excitation amplitudes tensor
        //Dimensions of o,v
        T1Tensor t1Amps(nOcc, nVirt);
        t1Amps.setZero();

        //The second excitation amplitudes tensor
        //Dimensions of o,o,v,v
        T2Tensor t2Amps = CreateInitialT2Tensor(nOcc, nVirt, moTEI, Dijab);

        //Variables for holding operations for acceleration with DIIS
        std::vector<Tensor<double, 2> > t1AmpsArray;
        t1AmpsArray.reserve(6);
        std::vector<Tensor<double, 4> > t2AmpsArray;
        t2AmpsArray.reserve(6);
        std::vector<VectorXd> DIIS_ErrorVectors;
        DIIS_ErrorVectors.reserve(6);

        double lastCorrelationEnergy = -DBL_MIN, correlationEnergy = 0.0;
        size_t ccIteration = 0;

        for (; ccIteration < MAX_ITERATIONS; ccIteration++) {
            t1AmpsArray.push_back(t1Amps);
            t2AmpsArray.push_back(t2Amps);

            //Solve the CC equations
            CoupledCluster::IterateCCSD(moTEI, F, nOcc, nVirt, Dia, Dijab, t1Amps, t2Amps);

            //Do DIIS
            if (!t1AmpsArray.empty() && !t2AmpsArray.empty()) {
                Tensor<double, 2> t1Err = t1Amps - t1AmpsArray.back();
                Tensor<double, 4> t2Err = t2Amps - t2AmpsArray.back();
                VectorXd combinedErrVec(t1Err.size() + t2Err.size());

                Eigen::Index idx = 0;
                for (Eigen::Index i = 0; i < nOcc; i++)
                    for (Eigen::Index a = 0; a < nVirt; a++)
                        combinedErrVec(idx++) = t1Err(i, a);

                for (Eigen::Index i = 0; i < nOcc; i++)
                    for (Eigen::Index j = 0; j < nOcc; j++)
                        for (Eigen::Index a = 0; a < nVirt; a++)
                            for (Eigen::Index b = 0; b < nVirt; b++)
                                combinedErrVec(idx++) = t2Err(i, j, a, b);


                DIIS_ErrorVectors.push_back(combinedErrVec);
            }

            if (ENABLE_DIIS && ccIteration >= 3) {
                if (ccIteration == 3 && console)
                    *console << "Turning on DIIS\n";

                if (DIIS_ErrorVectors.size() > MAX_DIIS_VECS) {
                    DIIS_ErrorVectors.erase(DIIS_ErrorVectors.begin());
                    t1AmpsArray.erase(t1AmpsArray.begin());
                    t2AmpsArray.erase(t2AmpsArray.begin());
                }

                Eigen::Index DIIS_Size = DIIS_ErrorVectors.size();

                MatrixXd ErrorMat = MatrixXd::Constant(DIIS_Size + 1, DIIS_Size + 1, -1.0);
                ErrorMat(ErrorMat.rows() - 1, ErrorMat.cols() - 1) = 0.0;

                double maxVal = 0.0;
                for (Eigen::Index i = 0; i < DIIS_ErrorVectors.size(); i++)
                    for (Eigen::Index j = 0; j < DIIS_ErrorVectors.size(); j++) {
                        ErrorMat(i, j) = DIIS_ErrorVectors[i].dot(DIIS_ErrorVectors[j]);
                        maxVal = std::max(maxVal, abs(ErrorMat(i, j)));
                    }

                ErrorMat.block(0, 0, ErrorMat.rows() - 1, ErrorMat.cols() - 1) /= maxVal;

                VectorXd Residual = VectorXd::Zero(DIIS_Size + 1);
                Residual(Residual.rows() - 1) = -1.0;

                VectorXd DIISCoeffs = ErrorMat.lu().solve(Residual);

                t1Amps.setZero();
                t2Amps.setZero();
                for (Eigen::Index i = 0; i < DIIS_Size; i++) {
                    t1Amps += DIISCoeffs(i) * t1AmpsArray[i];
                    t2Amps += DIISCoeffs(i) * t2AmpsArray[i];
                }
            }

            correlationEnergy = CoupledCluster::CalculateEnergyCCSD(F, nOcc, nVirt, moTEI, t1Amps, t2Amps);

            double dE = abs(correlationEnergy - lastCorrelationEnergy);
            if (console)
                *console << "   CCSD(" << ccIteration << "): E(corr) = " << correlationEnergy << "      dE = "
                        << dE
                        << '\n';

            if (dE < CC_CONVERGENCE) {
                if (console)
                    *console << "---- CCSD procedure converged! ----\n";
                break;
            }

            lastCorrelationEnergy = correlationEnergy;
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
        return CCSDResult{mp2Energy, correlationEnergy, t1Amps, t2Amps, moTEI};
    }
}
