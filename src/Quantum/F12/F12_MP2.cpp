#include<pch.h>

#include<unsupported/Eigen/CXX11/Tensor>
#include<Eigen/Core>

#include "F12.h"
#include "Quantum/Integrals/Integrals.h"
#include "Quantum/Integrals/IntegralsF12.h"
#include "Quantum/HF/HartreeFock.h"

namespace SCESS {
    using namespace Eigen;

    double F12::F12_MP2(const Molecule &mol, const RHFWavefunction &rhfWfn) {
        // const ConfigIntegrals hfBasisInts = rhfWfn.Integrals;
        //
        // TEITensor hfERIInts = Integrals::Transforms::AntisymmetriseIntegrateOutSpinClosedShell(
        //     Integrals::Transforms::ToMolecularTEI(hfBasisInts.tei, rhfWfn.Coeffs));
        //
        // //CABS should contain angular momentum functions at least 3*l{occ}
        //
        // double gamma = 0.5;
        //
        // TEITensor f12Ints = Integrals::Transforms::AntisymmetriseIntegrateOutSpinClosedShell(
        //     Integrals::Transforms::ToMolecularTEI(Integrals::F12::TwoElectronSlaterF12Tensor(mol, gamma, "sto-3g"),
        //                                           rhfWfn.Coeffs));
        // TEITensor f12SquaredInts = Integrals::Transforms::AntisymmetriseIntegrateOutSpinClosedShell(
        //     Integrals::Transforms::ToMolecularTEI(
        //         Integrals::F12::TwoElectronSlaterF12Tensor(mol, gamma * 2.0, "sto-3g"),
        //         rhfWfn.Coeffs));
        // TEITensor f12RepulsionInts = Integrals::Transforms::AntisymmetriseIntegrateOutSpinClosedShell(
        //     Integrals::Transforms::ToMolecularTEI(
        //         Integrals::F12::TwoElectronSlaterF12RepulsionTensor(mol, gamma, "sto-3g"),
        //         rhfWfn.Coeffs));

        return 0.0;
    }

    F12::Tensor4 F12::MakeBxyvwTensor() {
    }

    F12::Tensor4 F12::MakeXxyvwTensor() {
    }

    F12::Tensor4 F12::MakeCabxyTensor() {
    }

    F12::Tensor4 F12::MakeVijxyTensor() {
    }
}
