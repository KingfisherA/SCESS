#pragma once
#include<strstream>
#include<Eigen/Core>

#include "Chemistry/Quantum/Integrals/Integrals.h"

namespace SCESS
{
	class Molecule;
	class DFTMolecularGrid;

	class DFT
	{
	public:
		static double RunDFT(const Molecule& mol, std::stringstream* consoleStream = nullptr, const std::string_view basisSet = "");

		static Eigen::MatrixXd ConstructDensity(const Eigen::Index Nelecs, const Eigen::MatrixXd& coeffs);

		static Eigen::MatrixXd ConstructCoulombMatrix(const Eigen::MatrixXd& DensityMat, const TEITensor& tei);

		static void EvaluateExchangeCorrelation(const DFTMolecularGrid& dftGrid, const Molecule& mol,
			Eigen::MatrixXd& XCMatrix, double& energyOut);

		static double CalculateEnergy(const Eigen::MatrixXd& DensityMat,
			const Eigen::MatrixXd& Hcore,
			const Eigen::MatrixXd& JMat,
			const double nuclearNuclearE,
			const double exchangeCorEnergy);
	};
}