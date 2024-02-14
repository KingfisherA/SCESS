#include <pch.h>

#include "DFT.h"
#include<Eigen/Core>
#include "Quantum/Integrals/IntegralHelper.h"
#include "Quantum/Molecule.h"
#include "XCFunctional.h"
#include "DFTGrid.h"

//
//
namespace SCESS
{
	double DFT::RunDFT(const Molecule& mol, std::stringstream* consoleStream)
	{
		if (consoleStream)
			*consoleStream << "---Density Functional Theory---\n Using STO-3G basis with " << mol.GetOrbitals().size() << " basis set functions \nBeginning integral evaluation...\n";

		ConfigIntegrals integrals = IntegralHelper::CalculateConfiguration(mol);
		const Eigen::MatrixXd& S = integrals.S;
		const Eigen::MatrixXd& T = integrals.T;
		const Eigen::MatrixXd& V = integrals.V;
		const TETensor& tei = integrals.tei; //Two-electron tensor


		DFTMolecularGrid dftGrid(mol);
		if (consoleStream)
			*consoleStream << "Constructing DFT Grid with 15 radial points and a lebbedev order of 7\n";
		dftGrid.Construct(20, 7);
		if (consoleStream)
			*consoleStream << "DFT grid constructed!\n";

		if (consoleStream)
			*consoleStream << "Integral evaluation completed\n";

		Eigen::MatrixXd Hcore = T + V;

		Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver;
		eigensolver.compute(S);
		const Eigen::VectorXd& SEigenValues = eigensolver.eigenvalues().real();
		const Eigen::MatrixXd& SEigenVects = eigensolver.eigenvectors().real();
		Eigen::MatrixXd SValsInvSqrt = SEigenValues.cwiseInverse().cwiseSqrt().asDiagonal();
		Eigen::MatrixXd TransformMat = SEigenVects * SValsInvSqrt * SEigenVects.transpose();
		TransformMat = (TransformMat + TransformMat.transpose()) * 0.5;


		Eigen::MatrixXd DensityMat = Eigen::MatrixXd::Zero(S.rows(), S.rows());
		Eigen::MatrixXd XCMat = Eigen::MatrixXd::Zero(S.rows(), S.rows());
		Eigen::MatrixXd JMat = Eigen::MatrixXd::Zero(S.rows(), S.rows());
		Eigen::MatrixXd KSMat, Coeffs, OldDensity;

		double delta = 1.0,
			Vnn = mol.CalculateNuclearRepulsion(),
			totalE = 0.0,
			exchangeCorreE = 0.0;
		size_t scfIterations = 0, scfMaxIterations = 200;

		if (consoleStream)
			*consoleStream << "Beginning SCF iterations...\n";

		while (delta > 1E-5 && scfIterations < scfMaxIterations)
		{
			KSMat = Hcore + 2.0 * JMat + XCMat;
			KSMat = TransformMat.transpose() * KSMat * TransformMat;

			Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver2(KSMat);
			Eigen::MatrixXd CoeffsPrime = eigensolver2.eigenvectors().real();

			Coeffs = TransformMat * CoeffsPrime;

			OldDensity = DensityMat;
			DensityMat = ConstructDensity(mol.GetNoElectrons(), Coeffs);

			dftGrid.SetDensity(DensityMat);
			dftGrid.CorrectDensities();

			JMat = ConstructCoulombMatrix(DensityMat, tei);
			EvaluateExchangeCorrelation(dftGrid, mol, XCMat, exchangeCorreE);

			//Calculate delta change to the density matrix
			const Eigen::MatrixXd DensityMatriciesDelta = DensityMat - OldDensity;
			delta = std::sqrt(DensityMatriciesDelta.cwiseProduct(DensityMatriciesDelta).sum());

			totalE = CalculateEnergy(DensityMat, Hcore, JMat, Vnn, exchangeCorreE);

			// Calculate the current energy
			if (consoleStream) *consoleStream << "SCF(" << scfIterations << "); Delta=" << delta << ", E=" << totalE << " hartrees\n";
			scfIterations++;
		}

		if (consoleStream)
		{
			if (scfIterations < scfMaxIterations)
				*consoleStream << "SCF procedure achieved required convergence criteria!\n SCF converged after " << scfIterations << " iterations.\n DFT=" << totalE << " hartrees\n";
			else
				*consoleStream << "SCF procedure failed to achieve convergence criteria in " << scfMaxIterations << " iterations \n";
		}

		return totalE;
	}

	Eigen::MatrixXd DFT::ConstructDensity(const size_t Nelecs, const Eigen::MatrixXd& coeffs)
	{
		Eigen::MatrixXd DensityMat = Eigen::MatrixXd::Zero(coeffs.rows(), coeffs.rows());

		for (size_t i = 0; i < DensityMat.rows(); i++)
		{
			for (size_t j = 0; j < DensityMat.cols(); j++)
			{
				for (size_t k = 0; k < Nelecs / 2; k++)
				{
					DensityMat(i, j) += coeffs(i, k) * coeffs(j, k);
				}
			}
		}

		//Mix the old and new density matrices??
		return DensityMat;
	}
	Eigen::MatrixXd DFT::ConstructCoulombMatrix(const Eigen::MatrixXd& DensityMat, const TETensor& tei)
	{
		Eigen::MatrixXd J = Eigen::MatrixXd::Zero(DensityMat.rows(), DensityMat.rows());
		for (size_t i = 0; i < J.rows(); i++)
		{
			for (size_t j = 0; j < J.rows(); j++)
			{
				for (size_t k = 0; k < J.rows(); k++)
				{
					for (size_t l = 0; l < J.rows(); l++)
					{
						J(i, j) += DensityMat(k, l) * tei(i, j, k, l);
					}
				}
			}
		}
		return J;
	}

	void DFT::EvaluateExchangeCorrelation(const DFTMolecularGrid& dftGrid,
		const Molecule& mol,
		Eigen::MatrixXd& XCMatrix, double& energyOut)
	{
		Eigen::VectorXd densities = dftGrid.GetDensities();
		Eigen::VectorXd weights = dftGrid.GetWeights();
		Eigen::MatrixXd amplitudes = dftGrid.GetAmplitudes();

		Eigen::VectorXd exchangeE, exchangePotentialAlpha, exchangePotentialBeta;
		Eigen::VectorXd correlationE, correlationPotentialAlpha, correlationPotentialBeta;

		Eigen::VectorXd densityAlpha = densities * 0.5,
			densityBeta = densities * 0.5;

		XCFunctional::Exchange_AlphaX(densityAlpha, densityBeta, exchangeE, exchangePotentialAlpha, exchangePotentialBeta);
		XCFunctional::Correlation_VMW(densityAlpha, densityBeta, correlationE, correlationPotentialAlpha, correlationPotentialBeta);
		
		energyOut = weights.dot(exchangeE + correlationE);

		Eigen::VectorXd wva = weights.cwiseProduct((exchangePotentialAlpha + exchangePotentialBeta
			+ correlationPotentialAlpha + correlationPotentialBeta) * 0.5);

		for (size_t i = 0; i < mol.GetOrbitals().size(); i++) {
			Eigen::VectorXd row = amplitudes.row(i);
			Eigen::VectorXd wva_i = wva.cwiseProduct(row);
			for (size_t j = 0; j < mol.GetOrbitals().size(); j++) {
				XCMatrix(i, j) = wva_i.dot(amplitudes.row(j));
			}
		}
	}

	double DFT::CalculateEnergy(const Eigen::MatrixXd& DensityMat,
		const Eigen::MatrixXd& Hcore,
		const Eigen::MatrixXd& JMat,
		const double nuclearNuclearE,
		const double exchangeCorEnergy)
	{
		double result = 0.0;
		for (int i = 0; i < DensityMat.rows(); i++)
		{
			for (int j = 0; j < DensityMat.rows(); j++)
			{
				result += DensityMat(i, j) * (Hcore(i, j) + JMat(i, j));
			}
		}

		return result * 2.0 + nuclearNuclearE + exchangeCorEnergy;
	}
}