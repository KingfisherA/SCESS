#pragma once
//
#include "Quantum/Molecule.h"
#include "Quantum/Integrals/Grid/MolGrid.h"

namespace DH {
    class Thread;
}

namespace SCESS {
    class Molecule;
    class BasisSet;
}

namespace SCESS::Integrals {
    struct MolInts;
}

//
//
namespace SCESS::DFT {
    using TEITensor = Eigen::Tensor<double, 4>;
    //
    //
    //A structure representing a Kohn-Sham wavefunction
    struct KS_Wfn {
        //The total electronic energy of the Kohn-Sham calculation including the
        //nuclear-nuclear repulsion term Vnn
        double E_KS = 0.0;

        //The Kohn-Sham Fock
        Eigen::MatrixXd F;

        //The Kohn-Sham density matrix
        Eigen::MatrixXd D;

        //The Kohn-Sham molecular orbital coefficients
        Eigen::MatrixXd C;

        //The Kohn-Sham Fock matrix orbital eigenvalues
        Eigen::VectorXd Ep;

        //The alpha and beta electron occupation numbers
        Eigen::VectorXi Occ;

        //
        int GetTotalNoElectrons() const {
            int elec = 0;

            for (const auto &ele: Occ)
                elec += ele;

            return elec;
        }
    };
    //
    //
    //Perform Kohn-Sham DFT
    KS_Wfn KS_DFT(const Molecule &mol,
                  const BasisSet &basisSet,
                  Quadrature::MolGrid::Resolution gridRes,
                  const DH::Thread *thread = nullptr,
                  std::stringstream *console = nullptr);
    //
    //
}
