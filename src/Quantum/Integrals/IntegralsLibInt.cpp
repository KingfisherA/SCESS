#include <Quantum/QuantumPCH.h>
#include "IntegralsLibInt.h"
#include "Integrals.h"
#include "IntegralsF12.h"
#include "Quantum/Basis/BasisSet.h"
#include "Quantum/TensorOperations.h"

#include<any>

namespace SCESS::Integrals {
    using namespace libint2;
    using namespace Utils;

    using LIAtom = libint2::Atom;

#pragma region LibInt Init

    void LibInt::InitLibInt() {
        libint2::initialize(true);
    }

    void LibInt::CloseLibInt() {
        libint2::finalize();
    }

#pragma endregion

#pragma region Evaluation Functions

    static TEITensor EvaluateTwoElectronTensor(const libint2::BasisSet &basisSet, Operator op,
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

    static TEITensor EvaluateTwoElectronTensor(const F12BasisSet &basisSet, Operator op,
                                               const std::any &params = {}) {
        TEITensor teiTensor(
            (Eigen::Index) basisSet.GetNumOfFuncs(),
            (Eigen::Index) basisSet.GetNumOfFuncs(),
            (Eigen::Index) basisSet.GetNumOfFuncs(),
            (Eigen::Index) basisSet.GetNumOfFuncs());


        Engine integralEngine;
        if (params.has_value())
            integralEngine = Engine(op,
                                    basisSet.GetNumOfPrimitives(),
                                    basisSet.GetMaxL(),
                                    0,
                                    std::numeric_limits<double>::epsilon(),
                                    params);
        else
            integralEngine = Engine(op,
                                    basisSet.GetNumOfPrimitives(),
                                    basisSet.GetMaxL(),
                                    0,
                                    std::numeric_limits<double>::epsilon());

        const auto &shells = basisSet.GetShells();

        const auto &shell2bf = basisSet.ShellToBasisIdx();
        const auto &results = integralEngine.results();

        for (Eigen::Index s1 = 0; s1 < basisSet.GetNumShells(); s1++) {
            uint32_t n1 = shells[s1].size();
            uint32_t bf1 = shell2bf[s1];
            for (Eigen::Index s2 = 0; s2 < basisSet.GetNumShells(); s2++) {
                uint32_t n2 = shells[s2].size();
                uint32_t bf2 = shell2bf[s2];
                for (Eigen::Index s3 = 0; s3 < basisSet.GetNumShells(); s3++) {
                    uint32_t n3 = shells[s3].size();
                    uint32_t bf3 = shell2bf[s3];
                    for (Eigen::Index s4 = 0; s4 < basisSet.GetNumShells(); s4++) {
                        uint32_t n4 = shells[s4].size();
                        uint32_t bf4 = shell2bf[s4];

                        integralEngine.compute(
                            shells[s1], shells[s2], shells[s3], shells[s4]);
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

    static Eigen::MatrixXd EvaluateOneElectronMatrix(const libint2::BasisSet &basisSet,
                                                     Operator op, const std::any &params = {}) {
        Eigen::MatrixXd resultMatrix(basisSet.nbf(), basisSet.nbf());

        Engine integralEngine;
        if (params.has_value())
            integralEngine = Engine(op,
                                    basisSet.nbf(),
                                    basisSet.max_l(),
                                    0,
                                    std::numeric_limits<double>::epsilon(),
                                    params);
        else
            integralEngine = Engine(op,
                                    basisSet.nbf(),
                                    basisSet.max_l());

        const auto &shell2bf = basisSet.shell2bf();
        const auto &results = integralEngine.results();


        for (Eigen::Index s1 = 0; s1 < basisSet.size(); s1++) {
            uint32_t n1 = basisSet[s1].size();
            uint32_t bf1 = shell2bf[s1];
            for (Eigen::Index s2 = s1; s2 < basisSet.size(); s2++) {
                uint32_t n2 = basisSet[s2].size();
                uint32_t bf2 = shell2bf[s2];

                integralEngine.compute(basisSet[s1], basisSet[s2]);

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

    static Eigen::MatrixXd EvaluateOneElectronMatrix(const F12BasisSet &basisSet,
                                                     Operator op, const std::any &params = {}) {
        Eigen::MatrixXd resultMatrix(basisSet.GetNumOfFuncs(), basisSet.GetNumOfFuncs());

        Engine integralEngine;
        if (params.has_value())
            integralEngine = Engine(op,
                                    basisSet.GetNumOfPrimitives(),
                                    basisSet.GetMaxL(),
                                    0,
                                    std::numeric_limits<double>::epsilon(),
                                    params);
        else
            integralEngine = Engine(op,
                                    basisSet.GetNumOfPrimitives(),
                                    basisSet.GetMaxL());

        const auto &shell2bf = basisSet.ShellToBasisIdx();
        const auto &results = integralEngine.results();

        for (Eigen::Index s1 = 0; s1 < basisSet.GetNumShells(); s1++) {
            uint32_t n1 = basisSet[s1].size();
            uint32_t bf1 = shell2bf[s1];
            for (Eigen::Index s2 = 0; s2 < basisSet.GetNumShells(); s2++) {
                uint32_t n2 = basisSet[s2].size();
                uint32_t bf2 = shell2bf[s2];

                integralEngine.compute(basisSet[s1], basisSet[s2]);

                const double *integrals = results[0];

                for (uint32_t f1 = 0, f12 = 0; f1 < n1; f1++) {
                    uint32_t i = f1 + bf1;
                    for (uint32_t f2 = 0; f2 < n2; f2++, f12++) {
                        uint32_t j = f2 + bf2;
                        resultMatrix(i, j) = integrals[f12];
                    }
                }
            }
        }

        return resultMatrix;
    }

#pragma endregion

#pragma region IntegralsLibInt

    MolInts LibInt::CalculateConfiguration(const libint2::BasisSet &basisSet,
                                           const AtomListLI &atoms) {
        MolInts mints;

        mints.S = CalculateOverlapMatrix(basisSet);
        mints.T = CalculateKineticMatrix(basisSet);
        mints.V = CalculateNuclearMatrix(basisSet, atoms);
        mints.eri = TwoElectronRepulsionTensor(basisSet);

        return mints;
    }

    Eigen::MatrixXd LibInt::CalculateOverlapMatrix(const libint2::BasisSet &basisSet) {
        return EvaluateOneElectronMatrix(basisSet, Operator::overlap);
    }

    Eigen::MatrixXd LibInt::CalculateNuclearMatrix(const libint2::BasisSet &basisSet, const AtomListLI &atoms) {
        return EvaluateOneElectronMatrix(basisSet, Operator::nuclear, make_point_charges(atoms));
    }

    Eigen::MatrixXd LibInt::CalculateKineticMatrix(const libint2::BasisSet &basisSet) {
        return EvaluateOneElectronMatrix(basisSet, Operator::kinetic);
    }

    TEITensor LibInt::TwoElectronRepulsionTensor(const libint2::BasisSet &basisSet) {
        return EvaluateTwoElectronTensor(basisSet, Operator::coulomb);
    }

#pragma endregion

#pragma region F12 Integrals

    MolInts F12::EvaluateMolecularIntegrals(const F12BasisSet &basisSet, const AtomListLI &atoms) {
        MolInts mints;

        mints.S = CalculateOverlapMatrix(basisSet);
        mints.T = CalculateKineticMatrix(basisSet);
        mints.V = CalculateNuclearMatrix(basisSet, atoms);
        mints.eri = TwoElectronRepulsionTensor(basisSet);

        return mints;
    }

    MolInts F12::ExtractHFIntegrals(const MolInts &f12MolInts, long startIdx, long hfBasisSetSize) {
        MolInts hfInts;
        hfInts.S = f12MolInts.S.block(startIdx, startIdx, hfBasisSetSize, hfBasisSetSize);
        hfInts.T = f12MolInts.T.block(startIdx, startIdx, hfBasisSetSize, hfBasisSetSize);
        hfInts.V = f12MolInts.V.block(startIdx, startIdx, hfBasisSetSize, hfBasisSetSize);
        hfInts.eri = f12MolInts.eri.Slice4({startIdx, startIdx, startIdx, startIdx},
                                           {hfBasisSetSize, hfBasisSetSize, hfBasisSetSize, hfBasisSetSize});
        return hfInts;
    }

    Eigen::MatrixXd F12::CalculateOverlapMatrix(const F12BasisSet &basisSet) {
        return EvaluateOneElectronMatrix(basisSet, Operator::overlap);
    }

    Eigen::MatrixXd F12::CalculateNuclearMatrix(const F12BasisSet &basisSet,
                                                const std::vector<libint2::Atom> &atoms) {
        return EvaluateOneElectronMatrix(basisSet, Operator::nuclear, make_point_charges(atoms));
    }

    Eigen::MatrixXd F12::CalculateKineticMatrix(const F12BasisSet &basisSet) {
        return EvaluateOneElectronMatrix(basisSet, Operator::kinetic);
    }

    TEITensor F12::TwoElectronRepulsionTensor(const F12BasisSet &basisSet) {
        return EvaluateTwoElectronTensor(basisSet, Operator::coulomb);
    }

    F12::F12Integrals F12::EvaluateF12Integrals(const F12BasisSet &basisSet, double gamma) {
        LibIntCGF slaterGeminal = CalculateSlaterContraction(gamma);
        LibIntCGF slaterGeminalSquared = CalculateSlaterSquaredContraction(slaterGeminal);

        F12Integrals mints;
        mints.f12 = TwoElectronCGFF12Tensor(basisSet, slaterGeminal);
        mints.f12Sq = TwoElectronCGFF12Tensor(basisSet, slaterGeminalSquared);
        mints.f12DelSq = TwoElectronCGFDelF12SquaredTensor(basisSet, slaterGeminal);
        mints.f12Repulsion = TwoElectronCGFF12RepulsionTensor(basisSet, slaterGeminal);

        return mints;
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
                squaredCGF[i * cgfLength + j].first = slaterContraction[i].first + slaterContraction[j].first;
                squaredCGF[i * cgfLength + j].second = slaterContraction[i].second * slaterContraction[j].second;
            }

        return squaredCGF;
    }

    constexpr double gamma = 1.0;

    //F12 integrals
    TEITensor F12::TwoElectronCGFF12Tensor(const F12BasisSet &basisSet, const LibIntCGF &cgf) {
        return EvaluateTwoElectronTensor(basisSet, Operator::cgtg, cgf);
    }

    TEITensor F12::TwoElectronCGFF12RepulsionTensor(const F12BasisSet &basisSet, const LibIntCGF &cgf) {
        return EvaluateTwoElectronTensor(basisSet, Operator::cgtg_x_coulomb, cgf);
    }

    TEITensor F12::TwoElectronCGFDelF12SquaredTensor(const F12BasisSet &basisSet, const LibIntCGF &cgf) {
        return EvaluateTwoElectronTensor(basisSet, Operator::delcgtg2, cgf);
    }

#pragma endregion
}
