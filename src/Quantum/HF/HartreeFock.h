#pragma once

#include "Eigen/Dense"

namespace SCESS
{
	class GaussianCollection;
	struct HFResult;
	class Molecule;
	struct TETensor;

	class HF abstract
	{
	public:

		static HFResult ReducedHF(const Molecule& molecule, const double convergeLimit=1E-5, 
			const size_t scfMaxIterations=100, const bool consoleOutput = true);

		static double CalculateEnergy(const Eigen::MatrixXd& Density,
			const Eigen::MatrixXd& Hcore, const Eigen::MatrixXd& Fock);

	private:
		static Eigen::MatrixXd ConstructG(const Eigen::MatrixXd& Density, const TETensor& teiTensor);

		static Eigen::MatrixXd ConstructDensityRHF(const Eigen::MatrixXd& Coeffs,
			const size_t Nelecs);

		static double CalculateDelta(const Eigen::MatrixXd& oldDensity,
			const Eigen::MatrixXd& density);

		static double CalculateNuclearRepulsion(const Molecule& molecule);
	};

	struct HFResult
	{
		HFResult():Vnn(0), Energy(0){}
		HFResult(Eigen::MatrixXd Fock,
		Eigen::MatrixXd Hcore,
		Eigen::MatrixXd Density,
		Eigen::MatrixXd Coeffs,
		Eigen::VectorXd Energies,
		const double Vnn,
		const double Energy) :Fock(Fock), Hcore(Hcore),
			Density(Density), Coeffs(Coeffs), Energies(Energies),
			Vnn(Vnn), Energy(Energy) {}
		HFResult(const HFResult& res) :Fock(res.Fock), Hcore(res.Hcore), 
			Density(res.Density), Coeffs(res.Coeffs), Energies(res.Energies),
			Vnn(res.Vnn), Energy(res.Energy)
		{}
		Eigen::MatrixXd Fock;
		Eigen::MatrixXd Hcore;
		Eigen::MatrixXd Density;
		Eigen::MatrixXd Coeffs;
		Eigen::VectorXd Energies;
		const double Vnn;
		const double Energy;
	};
}