#include<pch.h>

#include "HartreeFock.h"
#include"Quantum/Basis/Gaussian.h"
#include"Quantum/Integrals/IntegralHelper.h"
#include"Quantum/Molecule.h"


//
//
namespace SCESS
{
	RHFResult HF::RestrictedHF(const Molecule& molecule,
		const double convergeLimit,
		const size_t scfMaxIterations,
		std::stringstream* consoleStream)
	{
		if (consoleStream)
			*consoleStream << "---Restricted Hartree Fock---\n Using STO-3G basis with " << molecule.GetOrbitals().size() << " basis set functions \nBeginning integral evaluation...\n";

		ConfigIntegrals integrals = IntegralHelper::CalculateConfiguration(molecule);
		const Eigen::MatrixXd& S = integrals.S;
		const Eigen::MatrixXd& T = integrals.T;
		const Eigen::MatrixXd& V = integrals.V;
		const TETensor& tei = integrals.tei; //Two-electron tensor

		if (consoleStream)
			*consoleStream << "Integral evaluation completed\n";

		const Eigen::Index dimensions = S.rows();
		double delta = 1.0;
		size_t scfCount = 0;
		double Vnn = CalculateNuclearRepulsion(molecule);
		double currentE = 0.0;

		Eigen::MatrixXd Hcore = T + V;

		Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver;
		eigensolver.compute(S);
		const Eigen::VectorXd& SEigenValues = eigensolver.eigenvalues().real();
		const Eigen::MatrixXd& SEigenVects = eigensolver.eigenvectors().real();
		Eigen::MatrixXd SValsInvSqrt = SEigenValues.cwiseInverse().cwiseSqrt().asDiagonal();
		Eigen::MatrixXd TransformMat = SEigenVects * SValsInvSqrt * SEigenVects.transpose();
		TransformMat = (TransformMat + TransformMat.transpose()) * 0.5;

		Eigen::MatrixXd Fock, Coeffs, Energies, OldDensity;
		Eigen::MatrixXd DensityMat = Eigen::MatrixXd::Zero(dimensions, dimensions);

		if (consoleStream)
			*consoleStream << "Beginning SCF iterations...\n";


		while (delta > convergeLimit && scfCount < scfMaxIterations)
		{
			scfCount++;
			//Generate the Fock matrix
			Fock = Hcore + ConstructGRHF(DensityMat, tei);
			Fock = (Fock + Fock.transpose()) * 0.5;
			//Transform the Fock matrix in the orthonormal basis set
			Eigen::MatrixXd FockPrime = TransformMat.transpose() * Fock * TransformMat;

			//Diagonalise the transformed Fock matrix
			Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver2(FockPrime);
			//eigensolver2.compute(FockPrime);
			Energies = eigensolver2.eigenvalues().real();

			//Transform coefficients back into normal basis
			const Eigen::MatrixXd& CPrime = eigensolver2.eigenvectors().real();
			Coeffs = TransformMat * CPrime;

			//Recompute the density matrix
			OldDensity = DensityMat;
			DensityMat = ConstructDensityRHF(Coeffs, molecule.GetNoElectrons());

			//Calculate delta change to the density matrix
			const Eigen::MatrixXd DensityMatriciesDelta = DensityMat - OldDensity;
			delta = std::sqrt(DensityMatriciesDelta.cwiseProduct(DensityMatriciesDelta).sum());

			currentE = CalculateEnergyRHF(DensityMat, Hcore, Fock) + Vnn;

			//Calculate the current energy
			if (consoleStream) *consoleStream << "SCF(" << scfCount << "); Delta=" << delta << ", E=" << currentE << " hartrees\n";
		}

		currentE = CalculateEnergyRHF(DensityMat, Hcore, Fock) + Vnn;
		if (consoleStream)
		{
			if (scfCount < scfMaxIterations)
				*consoleStream << "SCF procedure achieved required convergence criteria!\n SCF converged after " << scfCount << " iterations.\n Restricted Hartree Fock Energy=" << currentE << " hartrees\n";
			else
				*consoleStream << "SCF procedure failed to achieve convergence criteria in " << scfMaxIterations << " iterations \n";
		}

		Eigen::VectorXi Occupation = Eigen::VectorXi::Zero(Hcore.rows());

		for (long i = (molecule.GetNoElectrons() / 2) -1; i >= 0; --i)
		{
			Occupation(i) = 2;
		}


		return RHFResult(Fock, Hcore, DensityMat, Coeffs, Energies, Occupation, Vnn, currentE);
	}

	UHFResult HF::UnrestrictedHF(const Molecule& molecule, const double convergeLimit, const size_t scfMaxIterations,
		std::stringstream* consoleStream)
	{
		if (consoleStream != nullptr)
			*consoleStream << "---Unrestricted Hartree Fock---\n Using STO-3G basis with " << molecule.GetOrbitals().size() << " basis set functions \nBeginning integral evaluation...\n";

		//Spin-up electrons
		size_t alphaElectrons = molecule.GetNoElectrons() / 2;
		//Spin-down electrons
		size_t betaElectrons = molecule.GetNoElectrons() - alphaElectrons;

		ConfigIntegrals integrals = IntegralHelper::CalculateConfiguration(molecule);
		const Eigen::MatrixXd& S = integrals.S;
		const Eigen::MatrixXd& T = integrals.T;
		const Eigen::MatrixXd& V = integrals.V;
		const TETensor& tei = integrals.tei;

		const size_t dimensions = S.rows();

		if (consoleStream)
			*consoleStream << "Integral evaluation completed\n";

		double delta = 1.0;
		size_t scfCount = 0;
		double Vnn = CalculateNuclearRepulsion(molecule);
		double currentE = 0.0;

		const Eigen::MatrixXd& Hcore = T + V;

		Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver;
		eigensolver.compute(S);
		const Eigen::VectorXd& SEigenValues = eigensolver.eigenvalues().real();
		const Eigen::MatrixXd& SEigenVects = eigensolver.eigenvectors().real();
		Eigen::MatrixXd SValsInvSqrt = SEigenValues.cwiseInverse().cwiseSqrt().asDiagonal();
		Eigen::MatrixXd TransformMat = SEigenVects * SValsInvSqrt * SEigenVects.transpose();
		TransformMat = (TransformMat + TransformMat.transpose()) * 0.5;

		Eigen::MatrixXd DensityAlpha = Eigen::MatrixXd::Zero(dimensions, dimensions);
		Eigen::MatrixXd DensityBeta = Eigen::MatrixXd::Zero(dimensions, dimensions);
		Eigen::MatrixXd OldDensityAlpha, OldDensityBeta;
		Eigen::MatrixXd FockAlpha, FockBeta;
		Eigen::MatrixXd CoeffsAlpha, CoeffsBeta;
		Eigen::VectorXd EnergiesAlpha, EnergiesBeta;

		while (delta > convergeLimit && scfCount < scfMaxIterations)
		{
			const auto [GAlpha, GBeta] = ConstructGUHF(DensityAlpha, DensityBeta, tei);
			FockAlpha = Hcore + GAlpha;
			FockAlpha = (FockAlpha + FockAlpha.transpose()) * 0.5;
			FockBeta = Hcore + GBeta;
			FockBeta = (FockBeta + FockBeta.transpose()) * 0.5;

			Eigen::MatrixXd FockAlphaPrime = TransformMat.transpose() * FockAlpha * TransformMat;
			Eigen::MatrixXd FockBetaPrime = TransformMat.transpose() * FockBeta * TransformMat;

			Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> alphaEigen(FockAlphaPrime);
			Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> betaEigen(FockBetaPrime);

			EnergiesAlpha = alphaEigen.eigenvalues().real();
			EnergiesBeta = betaEigen.eigenvalues().real();

			CoeffsAlpha = TransformMat * alphaEigen.eigenvectors().real();
			CoeffsBeta = TransformMat * betaEigen.eigenvectors().real();

			OldDensityAlpha = DensityAlpha;
			OldDensityBeta = DensityBeta;

			DensityAlpha = Eigen::MatrixXd::Zero(Hcore.rows(), Hcore.rows());
			DensityBeta = Eigen::MatrixXd::Zero(Hcore.rows(), Hcore.rows());
			for (size_t i = 0; i < Hcore.rows(); ++i)
			{
				for (size_t j = 0; j < Hcore.cols(); ++j)
				{
					size_t k = 0;
					for (k = 0; k < alphaElectrons; ++k)
					{
						DensityAlpha(i, j) += CoeffsAlpha(i, k) * CoeffsAlpha(j, k);
					}

					for (k = 0; k < betaElectrons; ++k)
					{
						DensityBeta(i, j) += CoeffsBeta(i, k) * CoeffsBeta(j, k);
					}
				}
			}

			Eigen::MatrixXd DensityMatriciesDelta = DensityAlpha - OldDensityAlpha;
			delta = DensityMatriciesDelta.cwiseProduct(DensityMatriciesDelta).sum();
			DensityMatriciesDelta = DensityBeta - OldDensityBeta;
			delta += DensityMatriciesDelta.cwiseProduct(DensityMatriciesDelta).sum();
			delta = sqrt(delta);

			currentE = CalculateEnergyUHF(DensityAlpha,
				DensityBeta,
				Hcore,
				FockAlpha,
				FockBeta) + Vnn;

			if (consoleStream) *consoleStream << "SCF(" << scfCount << "); Delta=" << delta << ", E=" << currentE << " hartrees\n";

			scfCount++;
		}

		currentE = CalculateEnergyUHF(DensityAlpha,
			DensityBeta,
			Hcore,
			FockAlpha,
			FockBeta) + Vnn;

		if (consoleStream)
		{
			if (scfCount < scfMaxIterations)
				*consoleStream << "SCF procedure achieved required convergence criteria!\n SCF converged after " << scfCount << " iterations.\n Unrestricted Hartree Fock Energy=" << currentE << " hartrees\n";
			else
				*consoleStream << "SCF procedure failed to achieve convergence criteria in " << scfMaxIterations << " iterations \n";
		}

		Eigen::VectorXi OccupationAlpha = Eigen::VectorXi::Zero(Hcore.rows());
		Eigen::VectorXi OccupationBeta = Eigen::VectorXi::Zero(Hcore.rows());
		for (long i = alphaElectrons - 1; i >= 0; i--)
		{
			OccupationAlpha(i) = 1;
		}
		for (long i = betaElectrons - 1; i >= 0; i--)
		{
			OccupationBeta(i) = 1;
		}

		return UHFResult(FockAlpha, FockBeta, Hcore, DensityAlpha, DensityBeta, 
			CoeffsAlpha, CoeffsBeta, EnergiesAlpha, EnergiesBeta, OccupationAlpha, OccupationBeta, Vnn, currentE);
	}

	Eigen::MatrixXd HF::ConstructGRHF(const Eigen::MatrixXd& Density,
		const TETensor& teiTensor)
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
						//Coulomb integral
						double J = teiTensor(i, j, k, l);
						//Exchange integral
						double K = teiTensor(i, l, k, j);
						G(i, j) += Density(k, l) * (J - 0.5 * K);
					}
				}
			}
		}

		return G;
	}

	std::pair<Eigen::MatrixXd, Eigen::MatrixXd> HF::ConstructGUHF(const Eigen::MatrixXd& DensityAlpha,
		const Eigen::MatrixXd& DensityBeta, const TETensor& teiTensor)
	{
		Eigen::MatrixXd GAlpha = Eigen::MatrixXd::Zero(DensityAlpha.rows(), DensityAlpha.rows());
		Eigen::MatrixXd GBeta = Eigen::MatrixXd::Zero(DensityAlpha.rows(), DensityAlpha.rows());

		for (size_t i = 0; i < DensityAlpha.rows(); i++)
		{
			for (size_t j = 0; j < DensityAlpha.rows(); j++)
			{
				for (size_t k = 0; k < DensityAlpha.rows(); k++)
				{
					for (size_t l = 0; l < DensityAlpha.rows(); l++)
					{
						const double J = teiTensor(i, j, k, l);
						const double K = teiTensor(i, l, k, j);
						GAlpha(i, j) += DensityAlpha(k, l) * (J - K) + DensityBeta(k, l) * J;
						GBeta(i, j) += DensityBeta(k, l) * (J - K) + DensityAlpha(k, l) * J;
					}
				}
			}
		}
		return { GAlpha, GBeta };
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
					Density(i, j) += 2.0 * Coeffs(i, k) * Coeffs(j, k);
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
				delta += std::pow(density(i, j) - oldDensity(i, j), 2);
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
				if (r > 0.0) Vnn += atoms[i].GetCharge() * atoms[j].GetCharge() / r;
			}
		}
		return Vnn;
	}

	double HF::CalculateEnergyRHF(const Eigen::MatrixXd& Density,
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

	double HF::CalculateEnergyUHF(
		const Eigen::MatrixXd& DensityAlpha,
		const Eigen::MatrixXd& DensityBeta,
		const Eigen::MatrixXd& Hcore,
		const Eigen::MatrixXd& FockAlpha,
		const Eigen::MatrixXd& FockBeta)
	{
		double energy = 0.0;
		for (size_t i = 0; i < Hcore.rows(); i++)
		{
			for (size_t j = 0; j < Hcore.cols(); j++)
			{
				energy += (DensityAlpha(i, j) + DensityBeta(i, j)) * Hcore(i, j);
				energy += DensityAlpha(i, j) * FockAlpha(i, j) + DensityBeta(i, j) * FockBeta(i, j);
			}
		}
		energy *= 0.5;
		return energy;
	}

}