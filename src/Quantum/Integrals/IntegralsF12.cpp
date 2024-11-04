#include <pch.h>
#include "Integrals.h"
#include "IntegralsF12.h"
#include <libint2.hpp>

namespace SCESS::Integrals::F12 {
    using namespace libint2;
    using LIAtom = libint2::Atom;

    static TEITensor EvaluateTwoElectronF12Tensor(const Molecule &mol, Operator op, double gamma) {
        std::vector<LIAtom> atoms(mol.GetNoAtoms());

        for (uint32_t i = 0; i < mol.GetNoAtoms(); i++) {
            const SCESS::Atom &atm = mol.GetAtom(i);
            atoms[i] = {
                atm.GetCharge(),
                atm.GetPos().x(), atm.GetPos().y(), atm.GetPos().z()
            };
        }

        //Use sto-3g basis set in this example
        BasisSet basisSet("sto-3g", atoms);

        TEITensor f12Tensor(
            (Eigen::Index) basisSet.size(),
            (Eigen::Index) basisSet.size(),
            (Eigen::Index) basisSet.size(),
            (Eigen::Index) basisSet.size());

        //Slater type geminal operator
        Engine integralEngine(op, basisSet.max_nprim(), static_cast<int>(basisSet.max_l()));
        integralEngine.set_params(gamma);

        const auto &results = integralEngine.results();

        std::unordered_map<size_t, double> integralMap;
        integralMap.reserve(f12Tensor.dimension(0) * f12Tensor.dimension(1));

        for (Eigen::Index i = 0; i < f12Tensor.dimension(0); i++)
            for (Eigen::Index j = 0; j < f12Tensor.dimension(1); j++)
                for (Eigen::Index k = 0; k < f12Tensor.dimension(2); k++)
                    for (Eigen::Index l = 0; l < f12Tensor.dimension(3); l++) {
                        size_t idx = TEIYoshimineIdx(i, j, k, l);

                        auto it = integralMap.find(idx);

                        if (it == integralMap.end()) {
                            integralEngine.compute(
                                basisSet[i],
                                basisSet[j],
                                basisSet[k],
                                basisSet[l]);
                            const double *integralVal = results.front();
                            if (integralVal == nullptr)
                                continue;
                            integralMap.emplace(idx, *integralVal);
                            f12Tensor(i, j, k, l) = *integralVal;
                        } else
                            f12Tensor(i, j, k, l) = it->second;
                    }

        // for (Eigen::Index i = 0; i < f12Tensor.dimension(0); i++)
        //     for (Eigen::Index j = 0; j < f12Tensor.dimension(1); j++)
        //         for (Eigen::Index k = 0; k < f12Tensor.dimension(2); k++)
        //             for (Eigen::Index l = 0; l < f12Tensor.dimension(3); l++) {
        //                 integralEngine.compute(
        //                     basisSet[i],
        //                     basisSet[j],
        //                     basisSet[k],
        //                     basisSet[l]);
        //
        //                 const double *integrals = results[0];
        //                 if (integrals == nullptr)
        //                     continue;
        //
        //                 uint32_t n1 = basisSet[i].size();
        //                 uint32_t n2 = basisSet[j].size();
        //                 uint32_t n3 = basisSet[k].size();
        //                 uint32_t n4 = basisSet[l].size();
        //
        //                 for (uint32_t f1 = 0, f1234 = 0; f1 < n1; f1++)
        //                     for (uint32_t f2 = 0; f2 < n2; f2++)
        //                         for (uint32_t f3 = 0; f3 < n3; f3++)
        //                             for (uint32_t f4 = 0; f4 < n4; f4++, f1234++) {
        //                                 double integral = integrals[f1234];
        //                                 f12Tensor(i, j, k, l) = integral;
        //                             }
        //             }


        return f12Tensor;
    }

    TEITensor TwoElectronSlaterF12Tensor(const Molecule &mol, double gamma) {
        return EvaluateTwoElectronF12Tensor(mol, Operator::stg, gamma);
    }

    TEITensor TwoElectronSlaterF12RepulsionTensor(const Molecule &mol, double gamma) {
        return EvaluateTwoElectronF12Tensor(mol, Operator::yukawa, gamma);
    }
}
