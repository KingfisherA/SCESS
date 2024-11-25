#pragma once

#include "Resources/Resources.h"
#define DATADIR SCESS::Res::AssetsFolder::Get().GetResourcePath("Basis").string()
#include <libint2.hpp>

#include "Quantum/Basis/BasisUtils.h"

namespace SCESS::BasisUtils {
    class F12BasisSet;
}

namespace SCESS::Integrals::LibInt {
    using TEITensor = Eigen::Tensor<double, 4>;



    void InitLibInt();
    void CloseLibInt();

    Eigen::MatrixXd CalculateOverlapMatrix(const libint2::BasisSet &basisSet);
    Eigen::MatrixXd CalculateTotalOverlapMatrix(const BasisUtils::F12BasisSet& basisSet);

    Eigen::MatrixXd CalculateNuclearMatrix(const std::vector<libint2::Atom> &atoms,
                                                 const libint2::BasisSet &basisSet);
    Eigen::MatrixXd CalculateKineticMatrix(const libint2::BasisSet &basisSet);
    TEITensor TwoElectronRepulsionTensor(const libint2::BasisSet &basisSet);
    TEITensor TwoElectronRepulsionTensor(const BasisUtils::F12BasisSet &basisSet);

}
