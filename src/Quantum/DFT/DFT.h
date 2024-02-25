#pragma once
#include<strstream>
#include<Eigen/Core>

namespace SCESS
{
	class Molecule;
	struct TETensor;
	class DFTMolecularGrid;

	class DFT
	{
	public:
		static double RunDFT(const Molecule& mol, std::stringstream* consoleStream = nullptr, const std::string_view basisSet = "");

		static Eigen::MatrixXd ConstructDensity(const size_t Nelecs, const Eigen::MatrixXd& coeffs);

		static Eigen::MatrixXd ConstructCoulombMatrix(const Eigen::MatrixXd& DensityMat, const TETensor& tei);

		static void EvaluateExchangeCorrelation(const DFTMolecularGrid& dftGrid, const Molecule& mol,
			Eigen::MatrixXd& XCMatrix, double& energyOut);

		static double CalculateEnergy(const Eigen::MatrixXd& DensityMat,
			const Eigen::MatrixXd& Hcore,
			const Eigen::MatrixXd& JMat,
			const double nuclearNuclearE,
			const double exchangeCorEnergy);
	};
}