#pragma once
#include<Eigen/Core>
#include<unsupported/Eigen/CXX11/Tensor>

#include "Chemistry/Quantum/Integrals/IntegralHelper.h"

namespace SCESS {
    struct RHFResult;
}

namespace SCESS {

class CoupledCluster {

public:

    static void CCSD(const RHFResult& rhfResult);

private:

};

} // SCESS

