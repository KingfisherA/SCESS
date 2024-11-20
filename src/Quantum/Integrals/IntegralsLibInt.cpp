#include <pch.h>
#include "IntegralsLibInt.h"
#include "Integrals.h"
#include "IntegralsF12.h"

#include<any>

namespace SCESS::Integrals {
    using namespace libint2;
    using namespace Utils;

    using LIAtom = libint2::Atom;

    void LibInt::InitLibInt() {
        libint2::initialize(true);
    }

    void LibInt::CloseLibInt() {
        libint2::finalize();
    }

    static TEITensor EvaluateTwoElectronTensor(const BasisSet &basisSet, Operator op,
                                               const std::any &params = {}) {
        TEITensor teiTensor(
            (Eigen::Index) basisSet.nbf(),
            (Eigen::Index) basisSet.nbf(),
            (Eigen::Index) basisSet.nbf(),
            (Eigen::Index) basisSet.nbf());


        Engine integralEngine;
        if (params.has_value())
            integralEngine = Engine(op,
                                    basisSet.max_nprim(),
                                    static_cast<int>(basisSet.max_l()),
                                    0,
                                    std::numeric_limits<double>::epsilon(),
                                    params);
        else
            integralEngine = Engine(op,
                                    basisSet.max_nprim(),
                                    static_cast<int>(basisSet.max_l()));

        const auto &shell2bf = basisSet.shell2bf();
        const auto &results = integralEngine.results();

        for (Eigen::Index s1 = 0; s1 < basisSet.size(); s1++) {
            uint32_t n1 = basisSet[s1].size();
            uint32_t bf1 = shell2bf[s1];
            for (Eigen::Index s2 = 0; s2 < basisSet.size(); s2++) {
                uint32_t n2 = basisSet[s2].size();
                uint32_t bf2 = shell2bf[s2];
                for (Eigen::Index s3 = 0; s3 < basisSet.size(); s3++) {
                    uint32_t n3 = basisSet[s3].size();
                    uint32_t bf3 = shell2bf[s3];
                    for (Eigen::Index s4 = 0; s4 < basisSet.size(); s4++) {
                        uint32_t n4 = basisSet[s4].size();
                        uint32_t bf4 = shell2bf[s4];

                        integralEngine.compute(
                            basisSet[s1], basisSet[s2], basisSet[s3], basisSet[s4]);
                        const double *integrals = results[0];
                        if (integrals == nullptr)
                            continue;

                        for (uint32_t f1 = 0, f1234 = 0; f1 < n1; f1++) {
                            uint32_t i = f1 + bf1;
                            for (uint32_t f2 = 0; f2 < n2; f2++) {
                                uint32_t j = f2 + bf2;
                                for (uint32_t f3 = 0; f3 < n3; f3++) {
                                    uint32_t k = f3 + bf3;
                                    for (uint32_t f4 = 0; f4 < n4; f4++, ++f1234) {
                                        uint32_t l = f4 + bf4;
                                        teiTensor(i, j, k, l) =
                                                integrals[f1234];
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        return teiTensor;
    }

    static Eigen::MatrixXd
    EvaluateOneElectronMatrix(const BasisSet &basisSet, Operator op, const std::any &params = {}) {
        Eigen::MatrixXd resultMatrix(basisSet.nbf(), basisSet.nbf());

        Engine intEngine(op, basisSet.max_nprim(), static_cast<int>(basisSet.max_l()));

        const auto &shell2bf = basisSet.shell2bf();
        const auto &results = intEngine.results();


        for (Eigen::Index s1 = 0; s1 < basisSet.size(); s1++) {
            uint32_t n1 = basisSet[s1].size();
            uint32_t bf1 = shell2bf[s1];
            for (Eigen::Index s2 = s1; s2 < basisSet.size(); s2++) {
                uint32_t n2 = basisSet[s2].size();
                uint32_t bf2 = shell2bf[s2];

                intEngine.compute(basisSet[s1], basisSet[s2]);

                const double *integrals = results[0];

                for (uint32_t f1 = 0, f12 = 0; f1 < n1; f1++) {
                    uint32_t i = f1 + bf1;
                    for (uint32_t f2 = 0; f2 < n2; f2++, f12++) {
                        uint32_t j = f2 + bf2;
                        resultMatrix(i, j) = resultMatrix(j, i) = integrals[f12];
                    }
                }
            }
        }

        return resultMatrix;
    }

    Eigen::MatrixXd LibInt::CalculateOverlapMatrixLibInt(const BasisSet &basisSet) {
        return EvaluateOneElectronMatrix(basisSet, Operator::overlap);
    }

    Eigen::MatrixXd LibInt::CalculateNuclearMatrixLibInt(const std::vector<LIAtom> &atoms, const BasisSet &basisSet) {
        return EvaluateOneElectronMatrix(basisSet, Operator::nuclear, make_point_charges(atoms));
    }

    Eigen::MatrixXd LibInt::CalculateKineticMatrix(const libint2::BasisSet &basisSet) {
        return EvaluateOneElectronMatrix(basisSet, Operator::kinetic);
    }

    TEITensor LibInt::TwoElectronRepulsionTensor(const libint2::BasisSet &basisSet) {
        return EvaluateTwoElectronTensor(basisSet, Operator::coulomb);
    }

    F12::LibIntCGF F12::CalculateSlaterContraction(double gamma) {
        if (gamma != 1.0)
            throw std::runtime_error(
                "Ability to calculate contraction of exp(-gamma*r12) for any gamma which is not one is not yet supported.");

        //Alpha exponent followed by the Gaussian coefficient
        LibIntCGF cgf({
            {0.2209, 0.3144},
            {1.004, 0.3037},
            {3.622, 0.1681},
            {12.16, 0.09811},
            {45.87, 0.06024},
            {254.4, 0.03726}
        });

        for (auto &ele: cgf)
            ele.second *= -1.0 / gamma;

        return cgf;
    }

    F12::LibIntCGF F12::CalculateSlaterSquaredContraction(const LibIntCGF &slaterContraction) {
        size_t cgfLength = slaterContraction.size();
        LibIntCGF squaredCGF(cgfLength * cgfLength, {0.0, 0.0});

        for (uint32_t i = 0; i < slaterContraction.size(); i++)
            for (uint32_t j = 0; j < slaterContraction.size(); j++) {
                squaredCGF[i * cgfLength + j].first = slaterContraction[i].first * slaterContraction[j].first;
                squaredCGF[i * cgfLength + j].second = slaterContraction[i].second + slaterContraction[j].second;
            }

        return squaredCGF;
    }

    constexpr double gamma = 1.0;

    //F12 integrals
    TEITensor F12::TwoElectronCGFF12Tensor(const BasisSet &basisSet, const LibIntCGF &cgf) {
        return EvaluateTwoElectronTensor(basisSet, Operator::cgtg, cgf);
    }

    TEITensor F12::TwoElectronCGFF12RepulsionTensor(const BasisSet &basisSet, const LibIntCGF &cgf) {
        return EvaluateTwoElectronTensor(basisSet, Operator::cgtg_x_coulomb, cgf);
    }

    TEITensor F12::TwoElectronCGFDelF12SquaredTensor(const libint2::BasisSet &basisSet, const LibIntCGF &cgf) {
        return EvaluateTwoElectronTensor(basisSet, Operator::delcgtg2, cgf);
    }
}
