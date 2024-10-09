#pragma once
#include<Eigen/Core>
#include<unsupported/Eigen/CXX11/Tensor>

#include "Chemistry/Quantum/Integrals/Integrals.h"
#include<optional>

#include "Chemistry/Quantum/HF/HartreeFock.h"

namespace SCESS {
    struct RHFResult;
}

namespace SCESS {
    struct CCSDResult {
    public:
        double mp2CorrelationEnergy = 0.0;
        double ccsdCorrelationEnergy = 0.0;
        Eigen::Tensor<double, 2> t1Amps;
        Eigen::Tensor<double, 4> t2Amps;
        TEITensor moTEI;
    };

    class CoupledCluster {
    public:
        using T1Tensor = Eigen::Tensor<double, 2>;
        using T2Tensor = Eigen::Tensor<double, 4>;
        using FTensor = Eigen::Tensor<double, 2>;
        using WTensor = Eigen::Tensor<double, 4>;
        using TauTensor = Eigen::Tensor<double, 4>;
        using FockTensor = Eigen::Tensor<double, 2>;

        /**
        * An implementation to solve the coupled cluster equations
        * to single and double excitations for a closed shell molecule starting from a RHF determinant
         * @param rhfResult The RHF determinant
         */
        static std::optional<CCSDResult> CCSD(const RHFResult &rhfResult, std::stringstream *console = nullptr);

        /**
        * An implementation to solve the coupled cluster equations
        * to single and double excitations for a closed shell molecule starting from a RHF determinant with DIIS acceleration.
         * @param rhfResult The RHF determinant
         */
        static std::optional<CCSDResult> CCSD_DIIS(const RHFResult &rhfResult, std::stringstream *console = nullptr);

        static std::optional<CCSDResult> CCSD_T(const RHFResult& rhfResult, std::stringstream* console = nullptr);

        static double CalculateTriplesCorrection(const CCSDResult& ccsd, const FockTensor& F, const Eigen::Index nOcc, const Eigen::Index nVirt);

        static FockTensor TransformFockToMO(const RHFResult &rhfResult);
        static void CreateDenominatorTensors(const FockTensor &F, const Eigen::Index nOcc, const Eigen::Index nVirt,
                                             Eigen::Tensor<double, 2> &Dia, Eigen::Tensor<double, 4> &Dijab);
        static T2Tensor CreateInitialT2Tensor(const Eigen::Index nOcc, const Eigen::Index nVirt, const TEITensor &moTEI,
                                              const Eigen::Tensor<double, 4> &Dijab);
        static TauTensor MakeTildeTau(const T1Tensor &t1, const T2Tensor &t2);
        static TauTensor MakeTau(const T1Tensor &t1, const T2Tensor &t2);
        static FTensor MakeFae(const TEITensor &moTEI, const FockTensor &F, const Eigen::Index nOcc,
                               const Eigen::Index nVirt, const T1Tensor &t1, const T2Tensor &t2,
                               const TauTensor &tildeTau);
        static FTensor MakeFmi(const TEITensor &moTEI, const FockTensor &F, const Eigen::Index nOcc,
                               const Eigen::Index nVirt, const T1Tensor &t1, const T2Tensor &t2,
                               const TauTensor &tildeTau);
        static FTensor MakeFme(const TEITensor &moTEI, const FockTensor &F, const Eigen::Index nOcc,
                               const Eigen::Index nVirt, const T1Tensor &t1, const T2Tensor &t2);
        static WTensor MakeWmnij(const TEITensor &moTEI, const Eigen::Index nOcc,
                                 const Eigen::Index nVirt, const T1Tensor &t1, const T2Tensor &t2,
                                 const TauTensor &tau);
        static WTensor MakeWabef(const TEITensor &moTEI, const Eigen::Index nOcc,
                                 const Eigen::Index nVirt, const T1Tensor &t1, const T2Tensor &t2,
                                 const TauTensor &tau);
        static WTensor MakeWmbej(const TEITensor &moTEI, const Eigen::Index nOcc,
                                 const Eigen::Index nVirt, const T1Tensor &t1, const T2Tensor &t2);
        static void IterateCCSD(const TEITensor &moTEI, const FockTensor &F, const Eigen::Index nOcc,
                                const Eigen::Index nVirt, const Eigen::Tensor<double, 2> &Dia,
                                const Eigen::Tensor<double, 4> &Dijab, T1Tensor &t1, T2Tensor &t2);

        static double CalculateEnergyCCSD(const FockTensor &F, const Eigen::Index nOcc, const Eigen::Index nVirt,
                                      const TEITensor &moTEI,
                                      const T1Tensor &t1, const T2Tensor &t2);
    };
} // SCESS
