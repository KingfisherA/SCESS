#pragma once
#include<Eigen/Core>

namespace SCESS
{
	class Atom;
	class Molecule;

	struct DFTGridPoint
	{
	public:
		DFTGridPoint(const Eigen::Vector3d& rIn, const Eigen::Vector3d& rBeckeIn);

		void SetBasisFuncAmplitude(const Molecule& mol);

		void MultiplyWeight(double x)
		{
			mWeight *= x;
		}


		void SetWeight(double w)
		{
			mWeight = w;
		}

		void SetDensity(const Eigen::MatrixXd& DensityMat);
		void MultiplyDensity(double x)
		{
			mDensity *= x;
		}

		double GetWeight() const { return mWeight; }

		const Eigen::Vector3d& GetPosBecke() const { return mPosBecke; }
		const Eigen::Vector3d& GetPos() const { return mPos; }
		const double GetDensity() const { return mDensity; }
		const Eigen::VectorXd& GetBasisAmplitudes() const
		{
			return mBasisAmplitudes;
		}

	private:
		Eigen::Vector3d mPos, mPosBecke;
		Eigen::VectorXd mBasisAmplitudes;
		double mWeight = 0.0, mDensity = 0.0;
	};

	class DFTAtomicGrid
	{
	public:
		DFTAtomicGrid(const Molecule& molIn, const Atom& atomIn) :mMolecule(molIn), mAtom(atomIn) {}

		void Construct(size_t radialPoints, size_t angularPoints, size_t lebedevOffset);

		void SetDensity(const Eigen::MatrixXd& DensityMat);
		void AdjustDensity(double correctionFactor);

		double GetDensity();

		Eigen::VectorXd GetWeights() const;
		Eigen::VectorXd GetDensities() const;
		Eigen::MatrixXd GetPositions() const;
		Eigen::MatrixXd GetAmplitudes() const;

		inline size_t GetGridSize() const
		{
			return mPoints.size();
		}
	private:
		double mDensity = 0.0;
		bool mDensityInvalid = false;
		const Molecule& mMolecule;
		const Atom& mAtom;
		Eigen::VectorXd mRn;
		std::vector<DFTGridPoint> mPoints;
	};

	class DFTMolecularGrid
	{
	public:
		DFTMolecularGrid(const Molecule& molIn);

		void Construct(size_t radialPoints, size_t lebedevOrder);

		void SetDensity(const Eigen::MatrixXd& DensityMat);

		void CorrectDensities();

		Eigen::VectorXd GetWeights() const;
		Eigen::VectorXd GetDensities() const;
		Eigen::MatrixXd GetAmplitudes() const;

	private:
		size_t mGridSize = 0;
		const Molecule& mMolecule;
		std::vector<DFTAtomicGrid> mGrids;
		std::vector<double> mAtomicDensities;
	};
}