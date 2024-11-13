#include<pch.h>
#include"IntegralTransforms.h"

namespace SCESS::Integrals {
    TEITensor Transforms::ToMolecularTEI(const TEITensor &tei, const Eigen::MatrixXd &C) {
        //Note that the tei tensors are in chemists’ notation

        TEITensor transformed1(tei.dimensions());
        transformed1.setZero();
        TEITensor transformed2(tei.dimensions());
        transformed2.setZero();
        TEITensor transformed3(tei.dimensions());
        transformed3.setZero();
        TEITensor transformed4(tei.dimensions());
        transformed4.setZero();

        Eigen::Index dim = tei.dimension(0);

        for (Eigen::Index s = 0; s < dim; ++s) {
            for (Eigen::Index mu = 0; mu < dim; ++mu)
                for (Eigen::Index nu = 0; nu < dim; ++nu)
                    for (Eigen::Index kap = 0; kap < dim; ++kap)
                        for (Eigen::Index lam = 0; lam < dim; ++lam)
                            //From <mu nu | kap lam> to <mu nu | kap s>
                                //Einsum("mnkl,ls->mnks")
                                    transformed1(mu, nu, kap, s) += C(lam, s) * tei(mu, nu, kap, lam);

            for (Eigen::Index r = 0; r < dim; ++r) {
                for (Eigen::Index mu = 0; mu < dim; ++mu)
                    for (Eigen::Index nu = 0; nu < dim; ++nu)
                        for (Eigen::Index kap = 0; kap < dim; ++kap)
                            //From <mu nu | kap s> to <mu nu | r s>
                                //Einsum("mnks,kr->mnrs")
                                    transformed2(mu, nu, r, s) += C(kap, r) * transformed1(mu, nu, kap, s);

                for (Eigen::Index q = 0; q < dim; ++q) {
                    for (Eigen::Index mu = 0; mu < dim; ++mu)
                        for (Eigen::Index nu = 0; nu < dim; ++nu)
                            //From <mu nu | r s> to <mu q | r s>
                                //Einsum("mnrs,nq->mqrs")
                                    transformed3(mu, q, r, s) += C(nu, q) * transformed2(mu, nu, r, s);

                    for (Eigen::Index p = 0; p < dim; ++p) {
                        for (Eigen::Index mu = 0; mu < dim; ++mu)
                            //From <mu q | r s> to <p q | r s>
                                //Einsum("mqrs,mp->pqrs")
                                    transformed4(p, q, r, s) += C(mu, p) * transformed3(mu, q, r, s);
                    }
                }
            }
        }

        return transformed4;
    }

    TEITensor Transforms::AntisymmetriseTEI(const TEITensor &tei) {
        TEITensor transformed(tei.dimensions());

        for (Eigen::Index p = 0; p < transformed.dimension(0); ++p)
            for (Eigen::Index q = 0; q < transformed.dimension(1); ++q)
                for (Eigen::Index r = 0; r < transformed.dimension(2); ++r)
                    for (Eigen::Index s = 0; s < transformed.dimension(3); ++s)
                        transformed(p, q, r, s) = tei(p, r, q, s) - tei(p, s, q, r);

        return transformed;
    }

    TEITensor Transforms::AntisymmetriseIntegrateOutSpinClosedShell(const TEITensor &moTEI) {
        Eigen::DSizes<Eigen::Index, 4> dimensions = moTEI.dimensions();

        //Make the tensor twice as big to hold the spin-integrated integrals
        //Each integral in the moTEI actually represents two electrons
        for (Eigen::Index &dim: dimensions)
            dim *= 2;

        TEITensor transformed(dimensions);

        for (Eigen::Index p = 0; p < transformed.dimension(0); ++p)
            for (Eigen::Index q = 0; q < transformed.dimension(1); ++q)
                for (Eigen::Index r = 0; r < transformed.dimension(2); ++r)
                    for (Eigen::Index s = 0; s < transformed.dimension(3); ++s) {
                        //Only allow the integral to exist if the alpha is integrated with an alpha
                        //or beta with a beta (spin function orthogonality)
                        Eigen::Index coloumbMask = (p % 2 == r % 2) * (q % 2 == s % 2);
                        //One if the spin functions commute. Zero otherwise
                        Eigen::Index exchangeMask = (p % 2 == s % 2) * (q % 2 == r % 2);

                        //Convert from chemist's to physicists' notation and apply spin integration
                        transformed(p, q, r, s) =
                                moTEI(p / 2, r / 2, q / 2, s / 2) * static_cast<double>(coloumbMask)
                                - moTEI(p / 2, s / 2, q / 2, r / 2) * static_cast<double>(exchangeMask);
                    }

        return transformed;
    }

    TEITensor Transforms::ConvertChemistsToPhysicists(const TEITensor &tei) {
        TEITensor transformed(tei.dimensions());
        Eigen::Index dim = tei.dimension(0);

        for (Eigen::Index mu = 0; mu < dim; ++mu)
            for (Eigen::Index nu = 0; nu < dim; ++nu)
                for (Eigen::Index kap = 0; kap < dim; ++kap)
                    for (Eigen::Index lam = 0; lam < dim; ++lam)
                        transformed(mu, nu, kap, lam) = tei(mu, kap, nu, lam);

        return transformed;
    }
}