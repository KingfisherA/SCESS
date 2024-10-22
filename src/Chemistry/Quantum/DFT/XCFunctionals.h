#pragma once

#include<Eigen/Core>

namespace SCESS {
    class XCFunctionals {
    public:
        static void xalpha_x_functional(const Eigen::VectorXd &densitiesa,
                                        const Eigen::VectorXd &densitiesb,
                                        Eigen::VectorXd &ex,
                                        Eigen::VectorXd &vxa,
                                        Eigen::VectorXd &vxb);

        static void vwm_c_functional(const Eigen::VectorXd &densitiesa,
                                     const Eigen::VectorXd &densitiesb,
                                     Eigen::VectorXd &ec,
                                     Eigen::VectorXd &vca,
                                     Eigen::VectorXd &vcb);

        static double vwn_xx(double x, double b, double c);

        static double vwn_epsp(double x);

        static double vwn_epsf(double x);

        static double vwn_eps(double x, double a, double x0, double b, double c);

        static double vwn_depsp(double x);
        static double vwn_depsf(double x);

        static double vwn_deps(double x, double a, double x0, double b, double c);

        static double vwn_g(double z);

        static double vwn_dg(double z);
    };
} // SCESS
