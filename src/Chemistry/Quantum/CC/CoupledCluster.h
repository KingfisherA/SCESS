#pragma once
#include<Eigen/Core>
#include<unsupported/Eigen/CXX11/Tensor>

#include "Chemistry/Quantum/Integrals/Integrals.h"
#include<optional>

namespace SCESS {
    struct RHFResult;
}

namespace SCESS {
    struct CCSDResult {
    public:
        double mp2CorrelationEnergy =0.0;
        double ccsdCorrelationEnergy = 0.0;
        Eigen::Tensor<double, 2> t1Amps;
        Eigen::Tensor<double, 4> t2Amps;
    };

    class CoupledCluster {
    public:
        /**
        * A direct product decomposition implementation to solve the coupled cluster equations
        * to single and double excitations for a closed shell molecule starting from a RHF determinant
         * @param rhfResult The RHF determinant
         */
        static std::optional<CCSDResult> CCSD(const RHFResult &rhfResult, std::stringstream* console = nullptr);
    };
} // SCESS
