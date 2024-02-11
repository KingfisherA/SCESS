#pragma once
#include<Eigen/Core>

//Taken from https://github.com/ifilot/dftcxx/blob/master/src/functionals.cpp

namespace SCESS::XCFunctional
{
    void Exchange_AlphaX(const Eigen::VectorXd& densitiesa,
        const Eigen::VectorXd& densitiesb,
        Eigen::VectorXd& ex,
        Eigen::VectorXd& vxa,
        Eigen::VectorXd& vxb);

    void Correlation_VMW(const Eigen::VectorXd& densitiesa,
        const Eigen::VectorXd& densitiesb,
        Eigen::VectorXd& ec,
        Eigen::VectorXd& vca,
        Eigen::VectorXd& vcb);
}