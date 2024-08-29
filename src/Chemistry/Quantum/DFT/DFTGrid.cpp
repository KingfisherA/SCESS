#include <pch.h>

#include <numbers>
#include "DFTGrid.h"
#include "Chemistry/Quantum/Molecule.h"
#include "Maths/Quadrature.h"

namespace SCESS
{
	DFTGridPoint::DFTGridPoint(const Eigen::Vector3d &rIn, const Eigen::Vector3d &rBeckeIn)
		: mPos(rIn), mPosBecke(rBeckeIn)
	{
	}

	void DFTGridPoint::SetBasisFuncAmplitude(const Molecule &mol)
	{
		mBasisAmplitudes = Eigen::VectorXd::Zero(mol.GetNoOrbitals());

		for (size_t i = 0; i < mol.GetNoOrbitals(); i++)
		{
			mBasisAmplitudes[i] = EvaluateGaussianFunction<double>(mol.GetOrbital(i), mPos.x(), mPos.y(), mPos.z());
		}
	}

	void DFTGridPoint::SetDensity(const Eigen::MatrixXd &DensityMat)
	{
		mDensity = 2.0 * mBasisAmplitudes.dot(DensityMat * mBasisAmplitudes);
	}

	void DFTAtomicGrid::Construct(size_t radialPoints, size_t angularPoints, size_t lebedevOffset)
	{
		double f = std::numbers::pi / double(radialPoints + 1);
		mRn = Eigen::VectorXd::Zero(radialPoints);

		for (size_t i = 1; i <= radialPoints; i++)
		{
			// Gauss-Chebyshev
			double w = f * pow(sin(f * double(i)), 2.0);

			double x = cos(f * double(i));

			double r = (1.0 + x) / (1.0 - x);
			mRn(i - 1) = r;

			w = w / sqrt(1.0 - x * x) * 2.0 / ((1.0 - x) * (1.0 - x));

			for (size_t j = lebedevOffset; j < angularPoints + lebedevOffset; j++)
			{
				Eigen::Vector3d p2(
					Quadrature::lebedev_coefficients[j][0],
					Quadrature::lebedev_coefficients[j][1],
					Quadrature::lebedev_coefficients[j][2]);
				p2 *= r;

				DFTGridPoint &currentP = mPoints.emplace_back(mAtom.GetPos() + p2, p2);
				currentP.SetWeight(w * Quadrature::lebedev_coefficients[j][3]);
				currentP.SetBasisFuncAmplitude(mMolecule);
			}
		}

		for (size_t i = 0; i < mPoints.size(); i++)
		{
			mPoints[i].MultiplyWeight((mPoints[i].GetPosBecke()).squaredNorm() * 4.0 * std::numbers::pi);
		}
	}

	void DFTAtomicGrid::SetDensity(const Eigen::MatrixXd &DensityMat)
	{
		for (DFTGridPoint &point : mPoints)
		{
			point.SetDensity(DensityMat);
		}

		mDensityInvalid = true;
	}

	void DFTAtomicGrid::AdjustDensity(double correctionFactor)
	{
		for (DFTGridPoint &pnt : mPoints)
		{
			pnt.MultiplyDensity(correctionFactor);
		}

		mDensityInvalid = true;
	}

	double DFTAtomicGrid::GetDensity()
	{
		if (mDensityInvalid)
		{
			mDensity = 0.0;

			for (DFTGridPoint &pnt : mPoints)
			{
				mDensity += pnt.GetWeight() * pnt.GetDensity();
			}

			mDensityInvalid = false;
		}

		return mDensity;
	}

	Eigen::VectorXd DFTAtomicGrid::GetWeights() const
	{
		Eigen::VectorXd weights = Eigen::VectorXd::Zero(mPoints.size());

		for (size_t i = 0; i < mPoints.size(); i++)
		{
			const DFTGridPoint &pnt = mPoints[i];
			weights(i) = pnt.GetWeight();
		}

		return weights;
	}

	Eigen::VectorXd DFTAtomicGrid::GetDensities() const
	{
		Eigen::VectorXd densities = Eigen::VectorXd::Zero(mPoints.size());

		for (size_t i = 0; i < mPoints.size(); i++)
		{
			const DFTGridPoint &pnt = mPoints[i];
			densities(i) = pnt.GetDensity();
		}

		return densities;
	}

	Eigen::MatrixXd DFTAtomicGrid::GetPositions() const
	{
		Eigen::MatrixXd positions = Eigen::MatrixXd::Zero(mPoints.size(), 3);
		for (size_t i = 0; i < mPoints.size(); i++)
		{
			const DFTGridPoint &pnt = mPoints[i];
			for (size_t j = 0; j < 3; j++)
			{
				positions(i, j) = pnt.GetPos()[j];
			}
		}

		return positions;
	}

	Eigen::MatrixXd DFTAtomicGrid::GetAmplitudes() const
	{
		Eigen::MatrixXd amplitudes = Eigen::MatrixXd::Zero(mMolecule.GetOrbitals().size(), mPoints.size());

		for (unsigned int i = 0; i < mPoints.size(); i++)
		{
			Eigen::VectorXd basisFuncAmplitudes = mPoints[i].GetBasisAmplitudes();
			for (unsigned int j = 0; j < mMolecule.GetOrbitals().size(); j++)
			{
				amplitudes(j, i) = basisFuncAmplitudes(j);
			}
		}

		return amplitudes;
	}

	DFTMolecularGrid::DFTMolecularGrid(const Molecule &molIn) : mMolecule(molIn)
	{
	}

	void DFTMolecularGrid::Construct(size_t radialPoints, size_t lebedevOrder)
	{
		size_t angularPoints = Quadrature::num_lebedev_points[lebedevOrder];

		size_t lebedevOffset = 0;
		for (size_t i = 0; i < lebedevOrder; i++)
		{
			lebedevOffset += Quadrature::num_lebedev_points[i];
		}

		for (size_t i = 0; i < mMolecule.GetAtoms().size(); i++)
		{
			mAtomicDensities.push_back(0.0);

			mGrids.emplace_back(mMolecule, mMolecule.GetAtom(i));
			DFTAtomicGrid &grid = mGrids.back();
			grid.Construct(radialPoints, angularPoints, lebedevOffset);
			mGridSize += grid.GetGridSize();
		}
	}

	void DFTMolecularGrid::SetDensity(const Eigen::MatrixXd &DensityMat)
	{
		for (DFTAtomicGrid &grid : mGrids)
		{
			grid.SetDensity(DensityMat);
		}
	}

	void DFTMolecularGrid::CorrectDensities()
	{
		double sum = 0.0;
		int charge = 0;

		for (size_t i = 0; i < mGrids.size(); i++)
		{
			sum += mGrids[i].GetDensity();
			charge += mMolecule.GetAtom(i).GetCharge();
		}

		double correctionFactor = double(charge) / sum;

		for (size_t i = 0; i < mGrids.size(); i++)
		{
			mGrids[i].AdjustDensity(correctionFactor);
		}
	}

	Eigen::VectorXd DFTMolecularGrid::GetWeights() const
	{
		Eigen::VectorXd weights = Eigen::VectorXd::Zero(mGridSize);
		size_t i = 0;
		for (const DFTAtomicGrid &grid : mGrids)
		{
			Eigen::VectorXd atomicW = grid.GetWeights();

			for (size_t j = 0; j < atomicW.size(); j++)
			{
				weights(i + j) = atomicW(j);
			}

			i += atomicW.size();
		}
		return weights;
	}

	Eigen::VectorXd DFTMolecularGrid::GetDensities() const
	{
		Eigen::VectorXd densities = Eigen::VectorXd::Zero(mGridSize);
		size_t i = 0;
		for (const DFTAtomicGrid &grid : mGrids)
		{
			Eigen::VectorXd atomicD = grid.GetDensities();

			for (size_t j = 0; j < atomicD.size(); j++)
			{
				densities(i + j) = atomicD(j);
			}

			i += atomicD.size();
		}
		return densities;
	}

	Eigen::MatrixXd DFTMolecularGrid::GetAmplitudes() const
	{
		Eigen::MatrixXd amplitudes = Eigen::MatrixXd::Zero(mMolecule.GetOrbitals().size(), mGridSize);
		size_t i = 0;

		for (const DFTAtomicGrid &grid : mGrids)
		{
			const Eigen::MatrixXd &atomicA = grid.GetAmplitudes();

			for (size_t j = 0; j < grid.GetGridSize(); j++)
			{
				for (size_t k = 0; k < mMolecule.GetNoOrbitals(); k++)
				{
					amplitudes(k, i + j) = atomicA(k, j);
				}
			}

			i += grid.GetGridSize();
		}

		return amplitudes;
	}

}