#pragma once

#include<Eigen/Core>
#include <GraphicsEngine/Utilities/Memory.h>
#include<unsupported/Eigen/CXX11/Tensor>

namespace SCESS {
    class F12BasisSet;
    struct UHFWavefunction;
}

namespace SCESS::Integrals::F12 {
    struct F12Integrals;
}

namespace SCESS {
    struct ConfigIntegrals;
}

namespace SCESS::BasisUtils {
    class F12BasisSet;
}

namespace SCESS::Integrals {
    using TEITensor = Eigen::Tensor<double, 4>;

    namespace Transforms {
        /**
         * Transform the two-electron integrals (in chemist's notation) from the atomic orbital to the molecular orbital basis set
         * @param tei Two-electron integral tensor in atomic orbital basis (in chemist's notation)
         * @param C Wavefunction coefficients
         * @return Transformed TEI tensor (in chemist's notation)
         */
        TEITensor TransformTEI(const TEITensor &tei, const Eigen::MatrixXd &C);

        TEITensor TransformTEI(
            const TEITensor &tei,
            const Eigen::MatrixXd &A,
            const Eigen::MatrixXd &B,
            const Eigen::MatrixXd &C,
            const Eigen::MatrixXd &D);

        /**
         * Transforms all integrals (f12 and eri) in the OBS into MO integrals and:
         * - Project out linear dependence of the HF and CABS basis sets
         * - Makes CABS canonicalise the Fock operator
         * These are the requirements for the CABS to be used in orbital-sensitive CCSD-F12
         */
        void TransformF12Basis(const F12BasisSet &basisSet,
                               ConfigIntegrals &mints,
                               F12::F12Integrals &f12Ints,
                               uRef<UHFWavefunction> *hfWfnOut = nullptr);

        /**
         * Produce a tensor with < i j || k l>=<i j | k l> - <i j | l k>
         * Automatically converts chemist's to physicists notation
         * @param tei Two-electron integral tensor (in chemist's notation)
         * @return Two-electron integral tensor with antisymmetric integrals (in physicists' notation)
         */
        TEITensor AntisymmetriseTEI(const TEITensor &tei);
        /**
         * Transforms a molecular two electron integral tensor (in chemist's notation) into a tensor of antisymmetric
         * two-electron integrals with the spin functions integrated out.
         * A closed-shell wavefunction is assumed such that alpha and beta functions alternate.
         * @param moTEI The TEI tensor transformed into the MO basis (in chemist's notation)
         * @return A tensor with elements < i j || k l>= <i j | k l><alpha/beta ...> - <i j | l k><alpha/beta ...>
         */
        TEITensor AntisymmetriseIntegrateOutSpinClosedShell(const TEITensor &moTEI);
        //
        /**
         *Converts tei from chemists' to physicists' notation i.e. <i j | k l> = ( i k | j l)
         * @param tei The TEI tensor in chemist's notation
         * @return The TEI tensor in physicists' notation
         */
        TEITensor ConvertChemistsToPhysicists(const TEITensor &tei);
    }
};
