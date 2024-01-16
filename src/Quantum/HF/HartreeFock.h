#pragma once

#include "Eigen/Dense"

namespace SCESS
{
	class GaussianCollection;
	struct RHFResult;
	struct UHFResult;
	class Molecule;
	struct TETensor;

	class HF abstract
	{
	public:

		static RHFResult RestrictedHF(const Molecule& molecule, const double convergeLimit = 1E-5,
			const size_t scfMaxIterations = 100, std::stringstream* consoleStream = nullptr);

		static UHFResult UnrestrictedHF(const Molecule& molecule, const double convergeLimit = 1E-5,
			const size_t scfMaxIterations = 100, std::stringstream* consoleStream = nullptr);

		static double CalculateEnergyRHF(const Eigen::MatrixXd& Density,
			const Eigen::MatrixXd& Hcore, const Eigen::MatrixXd& Fock);

		static double CalculateEnergyUHF(
			const Eigen::MatrixXd& DensityAlpha,
			const Eigen::MatrixXd& DensityBeta,
			const Eigen::MatrixXd& Hcore,
			const Eigen::MatrixXd& FockAlpha,
			const Eigen::MatrixXd& FockBeta);

	private:
		static Eigen::MatrixXd ConstructGRHF(const Eigen::MatrixXd& Density, const TETensor& teiTensor);
		static std::pair< Eigen::MatrixXd, Eigen::MatrixXd> ConstructGUHF(const Eigen::MatrixXd& DensityAlpha,
			const Eigen::MatrixXd& DensityBeta, const TETensor& teiTensor);


		static Eigen::MatrixXd ConstructDensityRHF(const Eigen::MatrixXd& Coeffs,
			const size_t Nelecs);

		static double CalculateDelta(const Eigen::MatrixXd& oldDensity,
			const Eigen::MatrixXd& density);

		static double CalculateNuclearRepulsion(const Molecule& molecule);
	};

	struct RHFResult
	{
		RHFResult() :Vnn(0), Energy(0) {}
		RHFResult(
			const Eigen::MatrixXd& Fock,
			const Eigen::MatrixXd& Hcore,
			const Eigen::MatrixXd& Density,
			const Eigen::MatrixXd& Coeffs,
			const Eigen::VectorXd& Energies,
			const Eigen::VectorXi& Occupation,
			const double Vnn,
			const double Energy) :Fock(Fock), Hcore(Hcore),
			Density(Density), Coeffs(Coeffs), Energies(Energies),Occupation(Occupation),
			Vnn(Vnn), Energy(Energy) {}
		RHFResult(const RHFResult& res) :Fock(res.Fock), Hcore(res.Hcore),
			Density(res.Density), Coeffs(res.Coeffs), Energies(res.Energies), Occupation(res.Occupation),
			Vnn(res.Vnn), Energy(res.Energy)
		{}
		Eigen::MatrixXd Fock;
		Eigen::MatrixXd Hcore;
		Eigen::MatrixXd Density;
		Eigen::MatrixXd Coeffs;
		Eigen::VectorXd Energies;
		Eigen::VectorXi Occupation;
		const double Vnn;
		const double Energy;
	};

	struct UHFResult
	{
		UHFResult() :Vnn(0), Energy(0) {}
		UHFResult(
			const Eigen::MatrixXd& FockAlpha,
			const Eigen::MatrixXd& FockBeta,
			const Eigen::MatrixXd& Hcore,
			const Eigen::MatrixXd& DensityAlpha,
			const Eigen::MatrixXd& DensityBeta,
			const Eigen::MatrixXd& CoeffsAlpha,
			const Eigen::MatrixXd& CoeffsBeta,
			const Eigen::VectorXd& EnergiesAlpha,
			const Eigen::VectorXd& EnergiesBeta,
			const Eigen::VectorXi& OccupationAlpha,
			const Eigen::VectorXi& OccupationBeta,
			const double Vnn,
			const double Energy) :FockAlpha(FockAlpha), FockBeta(FockBeta), Hcore(Hcore),
			DensityAlpha(DensityAlpha), DensityBeta(DensityBeta), CoeffsAlpha(CoeffsAlpha),
			CoeffsBeta(CoeffsBeta), EnergiesAlpha(EnergiesAlpha), EnergiesBeta(EnergiesBeta),
			OccupationAlpha(OccupationAlpha), OccupationBeta(OccupationBeta),
			Vnn(Vnn), Energy(Energy) {}
		UHFResult(const UHFResult& res) :FockAlpha(res.FockAlpha), FockBeta(res.FockBeta), Hcore(res.Hcore),
			DensityAlpha(res.DensityAlpha), DensityBeta(res.DensityBeta), CoeffsAlpha(res.CoeffsAlpha),
			CoeffsBeta(res.CoeffsBeta), EnergiesAlpha(res.EnergiesAlpha), EnergiesBeta(res.EnergiesBeta),
			OccupationAlpha(res.OccupationAlpha), OccupationBeta(res.OccupationBeta),
			Vnn(res.Vnn), Energy(res.Energy)
		{}
		Eigen::MatrixXd FockAlpha, FockBeta;
		Eigen::MatrixXd Hcore;
		Eigen::MatrixXd DensityAlpha, DensityBeta;
		Eigen::MatrixXd CoeffsAlpha, CoeffsBeta;
		Eigen::VectorXd EnergiesAlpha, EnergiesBeta;
		Eigen::VectorXi OccupationAlpha, OccupationBeta;
		const double Vnn;
		const double Energy;
	};
}