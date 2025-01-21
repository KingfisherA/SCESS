#pragma once

#include "Resources/AssetFolder.h"
#define DATADIR SCESS::Res::AssetsFolder::Get().GetResourcePath("Basis").string()
#include <libint2.hpp>
#include <unsupported/Eigen/CXX11/Tensor>
#include "Quantum/Basis/BasisSet.h"


namespace SCESS {
    namespace BasisUtils {
        class F12BasisSet;
    }

    namespace Integrals {
        struct MolInts;
    }

    namespace Integrals::LibInt {
        using TEITensor = Eigen::Tensor<double, 4>;
        using AtomListLI = std::vector<libint2::Atom>;

        void InitLibInt();
        void CloseLibInt();

        MolInts CalculateConfiguration(const libint2::BasisSet &basisSet, const AtomListLI &atoms);

        Eigen::MatrixXd CalculateOverlapMatrix(const libint2::BasisSet &basisSet);
        Eigen::MatrixXd CalculateNuclearMatrix(const libint2::BasisSet &basisSet, const AtomListLI &atoms);
        Eigen::MatrixXd CalculateKineticMatrix(const libint2::BasisSet &basisSet);
        TEITensor TwoElectronRepulsionTensor(const libint2::BasisSet &basisSet);
    }
}