#pragma once
#include <deque>
#include <vector>
#include <Eigen/Dense>

//
//
namespace SCESS::HF {
    //
    //
    //Pulay DIIS extrapolation of Fock matrices. Each entry holds one Fock matrix per spin (one for RHF and DFT, two
    //for UHF) with its error FDS - SDF in the orthonormal basis; the error of an entry is the sum over its spins.
    class FockDIIS {
    public:
        explicit FockDIIS(size_t maxVectors = 6) : mMaxVectors(maxVectors) {
        }

        //Stores the Fock matrices and errors of an iteration, dropping the oldest entry beyond the maximum
        void Add(std::vector<Eigen::MatrixXd> focks, std::vector<Eigen::MatrixXd> errors) {
            mFocks.push_back(std::move(focks));
            mErrors.push_back(std::move(errors));
            while (mFocks.size() > mMaxVectors)
                DiscardOldest(1);
        }

        void DiscardOldest(size_t n) {
            for (size_t i = 0; i < n && !mFocks.empty(); i++) {
                mFocks.pop_front();
                mErrors.pop_front();
            }
        }

        [[nodiscard]] size_t Size() const { return mFocks.size(); }

        //The linear combination of the stored Fock matrices which minimises the norm of the combined error
        [[nodiscard]] std::vector<Eigen::MatrixXd> Extrapolate() const {
            using namespace Eigen;
            const Index n = (Index) mFocks.size();
            MatrixXd B = MatrixXd::Constant(n + 1, n + 1, -1.0);
            B(n, n) = 0.0;

            double maxVal = 0.0;
            for (Index i = 0; i < n; i++)
                for (Index j = i; j < n; j++) {
                    double dot = 0.0;
                    for (size_t s = 0; s < mErrors[i].size(); s++)
                        dot += mErrors[i][s].cwiseProduct(mErrors[j][s]).sum();
                    B(i, j) = B(j, i) = dot;
                    maxVal = std::max(maxVal, std::abs(dot));
                }
            if (maxVal > 0.0)
                B.topLeftCorner(n, n) /= maxVal;

            VectorXd rhs = VectorXd::Zero(n + 1);
            rhs(n) = -1.0;
            const VectorXd coeffs = B.fullPivLu().solve(rhs);

            std::vector<MatrixXd> result;
            for (size_t s = 0; s < mFocks.back().size(); s++) {
                MatrixXd F = MatrixXd::Zero(mFocks.back()[s].rows(), mFocks.back()[s].cols());
                for (Index i = 0; i < n; i++)
                    F += coeffs(i) * mFocks[i][s];
                result.push_back(std::move(F));
            }
            return result;
        }

    private:
        size_t mMaxVectors;
        std::deque<std::vector<Eigen::MatrixXd> > mFocks, mErrors;
    };
    //
    //
}
