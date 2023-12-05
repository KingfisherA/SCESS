#include<pch.h>

#include "HartreeFock.h"
#include"Quantum/Basis/Gaussian.h"
#include"Quantum/Integrals/IntegralHelper.h"
#include"Quantum/Molecule.h"


//
//
namespace SCESS
{
    HFResult HF::ReducedHF(const Molecule& molecule,
                       const double convergeLimit,
                       const size_t scfMaxIterations,
                       const bool consoleOutput)
    {
        if (consoleOutput)
            std::cout << "Hartree Fock procedure begun \n Using STO-3G basis with " << molecule.GetOrbitals().size() << " basis set functions \n";

        Eigen::MatrixXd S = IntegralHelper::CalculateOverlapMatrix(molecule);
        Eigen::MatrixXd T = IntegralHelper::CalculateKineticMatrix(molecule);
        Eigen::MatrixXd V = IntegralHelper::CalculateNuclearMatrix(molecule);
        TETensor twoElectronTensor = IntegralHelper::CalculateTwoElectronTensor(molecule);
        
        const Eigen::Index dimensions = S.rows();
        double delta = 1.0;
        size_t scfCount = 0;
        double Vnn = CalculateNuclearRepulsion(molecule);
        double currentE = 0.0;
        size_t Nelectrons = 0;

        for (const Atom& a : molecule.GetAtoms())
        {
            //Nuclear charge=Number of electrons for reduced hartree fock
            //Assume a closed shell molecule
            Nelectrons += a.GetCharge();
        }

        Eigen::MatrixXd Hcore = T + V;

        Eigen::EigenSolver<Eigen::MatrixXd> eigensolver;
        eigensolver.compute(S);
        Eigen::VectorXd SEigenValues = eigensolver.eigenvalues().real();
        Eigen::MatrixXd SEigenVects = eigensolver.eigenvectors().real();

        Eigen::MatrixXd SValsInvSqrt = SEigenValues.array().rsqrt().matrix().asDiagonal();
        Eigen::MatrixXd TransformMat = SEigenVects * SValsInvSqrt;

        //eigensolver.compute(TransformMat.transpose()*Hcore*TransformMat);
        Eigen::MatrixXd Fock;
        Eigen::MatrixXd Coeffs;// = TransformMat * eigensolver.eigenvectors().real();
        Eigen::MatrixXd DensityMat = Eigen::MatrixXd::Zero(dimensions, dimensions);// = ConstructDensityRHF(Coeffs, Nelectrons);
        Eigen::VectorXd Energies;
        Eigen::MatrixXd OldDensity;


        if (consoleOutput)
        {
            std::cout << "--Overlap Matrix, S--\n";
            std::cout << S << "\n\n";
            std::cout << "--Kinetic Energy Matrix, T--\n";
            std::cout << T << "\n\n";
            std::cout << "--Nuclear Attraction Matrix, V--\n";
            std::cout << V << "\n\n";
            std::cout << "--Two electron Tensor--\n";

            std::stringstream ss;
            size_t trunctatedDim = std::min<size_t>(dimensions, 2);
            for (int i = 0; i < trunctatedDim; i++)
            {
                for (int j = 0; j < trunctatedDim; j++)
                {
                    for (int k = 0; k < trunctatedDim; k++)
                    {
                        for (int l = 0; l < trunctatedDim; l++)
                        {
                            char integralType = (i!=k && j!=l) ? 'K' : 'J';
                            ss << integralType <<"; <" << i + 1 << " " << j + 1 << " | " << k + 1 << " " << l + 1 << "> " << twoElectronTensor(i, j, k, l) << "\n";
                        }
                    }
                }
            }
            std::cout << ss.str();

            if (trunctatedDim < (size_t)dimensions)
                std::cout << "......\n";

            std::cout << "--Core Hamiltonian Matrix, Hcore--\n";
            std::cout << Hcore << "\n\n";

            std::cout << "--Nuclear Repulsion, Vnn--\n";
            std::cout << Vnn << "\n\n";
        }


        while (delta > convergeLimit && scfCount < scfMaxIterations)
        {
            scfCount++;
            //Generate the Fock matrix
            Fock = Hcore + ConstructG(DensityMat, twoElectronTensor);
            //Transform the Fock matrix in the orthonormal basis set
            Eigen::MatrixXd FockPrime = TransformMat.transpose() * Fock * TransformMat;

            //Diagonalise the transformed Fock matrix
            eigensolver.compute(FockPrime);
            Energies = eigensolver.eigenvalues().real();

            //Transform coefficients back into normal basis
            Coeffs = TransformMat * eigensolver.eigenvectors().real();

            //Recompute the density matrix
            OldDensity = DensityMat;
            DensityMat = ConstructDensityRHF(Coeffs, Nelectrons);

            //Calculate delta change to the density matrix
            delta = CalculateDelta(OldDensity, DensityMat);

            currentE = CalculateEnergy(DensityMat, Hcore, Fock) + Vnn;

            //Calculate the current energy
            if(consoleOutput) std::cout << "SCF(" << scfCount << "); Delta=" << delta << ", E=" << currentE << " hartrees\n";
        }

        currentE = CalculateEnergy(DensityMat, Hcore, Fock) + Vnn;
        if (consoleOutput)
        {
            if (scfCount < scfMaxIterations)
                std::cout << "SCF procedure converged after " << scfCount << " iterations.\n HF Energy=" << currentE << " hartrees\n";
            else
                std::cout << "SCF procedure failed to achieve convergence criteria in "<<scfMaxIterations<<" iterations \n";
        }

        HFResult result(Fock, Hcore, DensityMat, Coeffs, Energies, Vnn, currentE);
        return result;
    }

    Eigen::MatrixXd HF::ConstructG(const Eigen::MatrixXd& Density, const TETensor& teiTensor)
    {
        Eigen::MatrixXd G = Eigen::MatrixXd::Zero(Density.rows(), Density.rows());

        for (Eigen::DenseIndex i = 0; i < G.rows(); i++)
        {
            for (Eigen::DenseIndex j = 0; j < G.rows(); j++)
            {
                for (Eigen::DenseIndex k = 0; k < G.rows(); k++)
                {
                    for (Eigen::DenseIndex l = 0; l < G.rows(); l++)
                    {
                        //F += P*(J-0.5K)
                        G(i, j) += Density(k, l) * (teiTensor(i, j, l, k) - 0.5 * teiTensor(i, k, l, j));
                    }
                }
            }
        }

        return G;
    }

    Eigen::MatrixXd HF::ConstructDensityRHF(const Eigen::MatrixXd& Coeffs,
        const size_t Nelecs)
    {
        Eigen::MatrixXd Density = Eigen::MatrixXd::Zero(Coeffs.rows(), Coeffs.rows());
        for (int i = 0; i < Density.rows(); i++)
        {
            for (int j = 0; j < Density.rows(); j++)
            {
                for (int k = 0; k < (Nelecs / 2); k++)
                {
                    Density(i, j) += 2 * Coeffs(i, k) * Coeffs(j, k);
                }
            }
        }
        return Density;
    }

    double HF::CalculateDelta(const Eigen::MatrixXd& oldDensity,
        const Eigen::MatrixXd& density)
    {
        double delta = 0.0;
        for (int i = 0; i < density.rows(); i++)
        {
            for (int j = 0; j < density.rows(); j++)
            {
                delta += std::pow(density(i, j) - oldDensity(i,j), 2);
            }
        }

        return std::sqrt(delta);
    }

    double HF::CalculateNuclearRepulsion(const Molecule& molecule)
    {
        double Vnn = 0.0;
        const std::vector<Atom>& atoms = molecule.GetAtoms();
        for (size_t i = 0; i < atoms.size(); i++) {
            for (size_t j = i + 1; j < atoms.size(); j++)
            {
                double r = (atoms[i].GetPos() - atoms[j].GetPos()).norm();
                if(r>0.0) Vnn += atoms[i].GetCharge() * atoms[j].GetCharge() / r;
            }
        }
        return Vnn;
    }

    double HF::CalculateEnergy(const Eigen::MatrixXd& Density,
        const Eigen::MatrixXd& Hcore, const Eigen::MatrixXd& Fock)
    {
        double result = 0.0;

        for (int i = 0; i < Fock.rows(); i++)
        {
            for (int j = 0; j < Fock.rows(); j++)
            {
                result += Density(i, j) * (Hcore(i, j) + Fock(i, j));
            }
        }

        return 0.5 * result;
    }
}