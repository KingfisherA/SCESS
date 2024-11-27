#include <pch.h>
#include <Eigen/Core>
#include <unsupported/Eigen/CXX11/Tensor>

#include "GraphicsEngine/Utilities/Threading.h"

#include "Quantum/CC/CoupledCluster.h"
#include "Quantum/Integrals/IntegralTransforms.h"

namespace SCESS {
    using namespace Eigen;

    uRef<CCSDWavefunction> CoupledCluster::CCSD_DIIS(const RHFWavefunction &rhfResult, const GEng::Thread *thread,
                                                     std::stringstream *console) {
        //The delta change in energy which must be achieved for convergence
        constexpr double CC_CONVERGENCE = 1.0E-6;
        //Maximum coupled cluster iterations
        constexpr size_t MAX_ITERATIONS = 50;
        //Maximum number of t amplitudes that will be stored to create the DIIS error
        //vector
        constexpr size_t MAX_DIIS_VECS = 6;

        //Transform the two electron repulsion integrals into the MO basis set
        //for MP2 calculation. Integrals are in chemist's notation
        TEITensor moTEI = Integrals::Transforms::TransformTEI(rhfResult.Integrals.eri, rhfResult.C);
        double mp2Energy = HartreeFock::MP2ClosedShell(rhfResult, moTEI, thread, console);

        if (console) {
            *console << std::setprecision(12);
            *console << "--- Coupled Cluster Singles & Doubles (CCSD) with DIIS acceleration ---\n";
        }

        //Antisymmetrise the two-electron integrals and integrate out the alpha and beta
        //spin functions for a closed-shell wavefunction. Also converts from chemist's
        //to physicits' notation
        moTEI = Integrals::Transforms::AntisymmetriseIntegrateOutSpinClosedShell(moTEI);

        //Number of occupied spin orbitals
        Eigen::Index nOcc = std::distance(rhfResult.Occupation.begin(),
                                          std::ranges::find(rhfResult.Occupation, 0));
        //For two electrons in each orbital
        nOcc *= 2;
        //Number of virtual spin orbitals
        Eigen::Index nVirt = rhfResult.Occupation.size() * 2 - nOcc;

        //Fock matrix in canonical MO basis set
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
            if(thread)
                if(thread->ShouldStop())
                    return nullptr;

            t1AmpsArray.push_back(t1Amps);
            t2AmpsArray.push_back(t2Amps);

            //Solve the CC equations
            CoupledCluster::IterateCCSD(moTEI, F, nOcc, nVirt, Dia, Dijab, t1Amps, t2Amps, thread);

            if(thread)
                if(thread->ShouldStop())
                    return nullptr;

            //Store the amplitudes to be used in DIIS calculation
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

            if (ccIteration >= 3) {
                if (ccIteration == 3 && console)
                    *console << "Turning on DIIS\n";

                //Remove vector from the end if we have stored more than
                //MAX_DIIS_VECS
                if (DIIS_ErrorVectors.size() > MAX_DIIS_VECS) {
                    DIIS_ErrorVectors.erase(DIIS_ErrorVectors.begin());
                    t1AmpsArray.erase(t1AmpsArray.begin());
                    t2AmpsArray.erase(t2AmpsArray.begin());
                }

                Eigen::Index DIIS_Size = DIIS_ErrorVectors.size();

                MatrixXd ErrorMat = MatrixXd::Constant(DIIS_Size + 1, DIIS_Size + 1, -1.0);
                ErrorMat(ErrorMat.rows() - 1, ErrorMat.cols() - 1) = 0.0;

                //Construct the DIIS error equations
                double maxVal = 0.0;
                for (Eigen::Index i = 0; i < DIIS_ErrorVectors.size(); i++)
                    for (Eigen::Index j = 0; j < DIIS_ErrorVectors.size(); j++) {
                        ErrorMat(i, j) = DIIS_ErrorVectors[i].dot(DIIS_ErrorVectors[j]);
                        maxVal = std::max(maxVal, abs(ErrorMat(i, j)));
                    }

                ErrorMat.block(0, 0, ErrorMat.rows() - 1, ErrorMat.cols() - 1) /= maxVal;

                VectorXd Residual = VectorXd::Zero(DIIS_Size + 1);
                Residual(Residual.rows() - 1) = -1.0;

                //Solve the set of linear equations to obtain coefficients to construct
                //next guesses for t1 and t2
                VectorXd DIISCoeffs = ErrorMat.lu().solve(Residual);

                //Construct next guess for t1 and t2
                t1Amps.setZero();
                t2Amps.setZero();
                for (Eigen::Index i = 0; i < DIIS_Size; i++) {
                    t1Amps += DIISCoeffs(i) * t1AmpsArray[i];
                    t2Amps += DIISCoeffs(i) * t2AmpsArray[i];
                }
            }

            //Calculate CC correlation energy
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
            return nullptr;
        }

        if (console) {
            *console << "E(CCSD correlation) = " << correlationEnergy << " hartrees\n";
            *console << "E(Total) = E(Hartree Fock) + E(CCSD) = " << rhfResult.Energy + correlationEnergy <<
                    " hartrees\n\n\n";
        }
        return std::make_unique<CCSDWavefunction>(mp2Energy, correlationEnergy, t1Amps, t2Amps, moTEI);
    }
}
