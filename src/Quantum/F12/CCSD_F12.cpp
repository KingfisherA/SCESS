#include<pch.h>
#include<unsupported/Eigen/CXX11/Tensor>
#include<Eigen/Core>

#include "F12.h"
#include "Quantum/Integrals/Integrals.h"
#include "Quantum/Integrals/IntegralsF12.h"
#include "Quantum/Integrals/IntegralTransforms.h"
#include "Quantum/HF/HartreeFock.h"

namespace SCESS::F12 {
    using namespace Eigen;

    double CCSD::CCSD_F12(const UHFWavefunction &refWfn,
                         Integrals::F12::F12Integrals &f12Ints) {
        f12Ints.eri = Integrals::Transforms::AntisymmetriseIntegrateOutSpinClosedShell(f12Ints.eri);
        f12Ints.f12 = Integrals::Transforms::AntisymmetriseIntegrateOutSpinClosedShell(f12Ints.f12);
        f12Ints.f12Sq = Integrals::Transforms::AntisymmetriseIntegrateOutSpinClosedShell(f12Ints.f12Sq);
        f12Ints.f12DelSq = Integrals::Transforms::AntisymmetriseIntegrateOutSpinClosedShell(f12Ints.f12DelSq);
        f12Ints.f12Repulsion = Integrals::Transforms::AntisymmetriseIntegrateOutSpinClosedShell(f12Ints.f12Repulsion);


        return 0.0;
    }
}
