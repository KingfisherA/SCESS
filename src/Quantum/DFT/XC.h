#pragma once
//
#include<Eigen/Core>

//
//
namespace SCESS::DFT::XC {
    //
    //
    void xalpha_x_functional(const Eigen::VectorXd &densitiesa,
                             const Eigen::VectorXd &densitiesb,
                             Eigen::VectorXd &ex,
                             Eigen::VectorXd &vxa,
                             Eigen::VectorXd &vxb);
    //
    void vwm_c_functional(const Eigen::VectorXd &densitiesa,
                          const Eigen::VectorXd &densitiesb,
                          Eigen::VectorXd &ec,
                          Eigen::VectorXd &vca,
                          Eigen::VectorXd &vcb);
    //
    double vwn_xx(double x, double b, double c);
    //
    double vwn_epsp(double x);
    //
    double vwn_epsf(double x);
    //
    double vwn_eps(double x, double a, double x0, double b, double c);
    //
    double vwn_depsp(double x);
    double vwn_depsf(double x);
    //
    double vwn_deps(double x, double a, double x0, double b, double c);
    //
    double vwn_g(double z);
    //
    double vwn_dg(double z);
    //
    //
} // SCESS
