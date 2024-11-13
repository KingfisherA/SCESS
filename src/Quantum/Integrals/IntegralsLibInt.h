#pragma once

#include "Resources/Resources.h"
#define DATADIR SCESS::Res::AssetsFolder::Get().GetResourcePath("Basis").string()
#include <libint2.hpp>
#undef DATADIR

namespace SCESS::Integrals::LibInt {
    using TEITensor = Eigen::Tensor<double, 4>;

    void InitLibInt();
    void CloseLibInt();

    Eigen::MatrixXd CalculateOverlapMatrixLibInt(const libint2::BasisSet &basisSet);
    Eigen::MatrixXd CalculateNuclearMatrixLibInt(const std::vector<libint2::Atom> &atoms,
                                                 const libint2::BasisSet &basisSet);
    Eigen::MatrixXd CalculateKineticMatrix(const libint2::BasisSet &basisSet);
    TEITensor TwoElectronRepulsionTensor(const libint2::BasisSet &basisSet);
}
