#pragma once
#include<strstream>
#include<Eigen/Core>

namespace SCESS
{
	class Molecule;
	class TETensor;
	class DFTMolecularGrid;

	class DFT
	{
	public:
		static double RunDFT(const Molecule& mol, std::stringstream* consoleStream = nullptr);

		static Eigen::MatrixXd ConstructDensity(const size_t Nelecs, const Eigen::MatrixXd& coeffs);

		static Eigen::MatrixXd ConstructCoulombMatrix(const Eigen::MatrixXd& DensityMat, const TETensor& tei);

		static void EvaluateExchangeCorrelation(const DFTMolecularGrid& dftGrid, const Molecule& mol,
			Eigen::MatrixXd& XCMatrix, double& energyOut);
	};
}