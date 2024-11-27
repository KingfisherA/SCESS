#pragma once

#include "Resources/Resources.h"
#define DATADIR SCESS::Res::AssetsFolder::Get().GetResourcePath("Basis").string()
#include <libint2.hpp>

#include "Quantum/Basis/BasisUtils.h"

namespace SCESS {
    struct ConfigIntegrals;
}

namespace SCESS::BasisUtils {
    class F12BasisSet;
}

namespace SCESS::Integrals::LibInt {
    using TEITensor = Eigen::Tensor<double, 4>;
    using AtomListLI = std::vector<libint2::Atom>;

    void InitLibInt();
    void CloseLibInt();

    ConfigIntegrals CalculateConfiguration(const libint2::BasisSet &basisSet, const AtomListLI &atoms);

    Eigen::MatrixXd CalculateOverlapMatrix(const libint2::BasisSet &basisSet);
    Eigen::MatrixXd CalculateNuclearMatrix(const libint2::BasisSet &basisSet, const AtomListLI &atoms);
    Eigen::MatrixXd CalculateKineticMatrix(const libint2::BasisSet &basisSet);
    TEITensor TwoElectronRepulsionTensor(const libint2::BasisSet &basisSet);
}
