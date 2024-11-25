#pragma once

#include<Eigen/Core>
#include<unsupported/Eigen/CXX11/Tensor>

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
        TEITensor ToMolecularTEI(const TEITensor &tei, const Eigen::MatrixXd &C);

        /**
         * Transform the two-electron integrals (in chemist's notation) from the atomic orbital to the molecular orbital basis set
         * @param tei Two-electron integral tensor in atomic orbital basis (in chemist's notation)
         * @param C Wavefunction coefficients
         * @return Transformed TEI tensor (in chemist's notation)
         */
        TEITensor ToMolecularTEI(const TEITensor &tei, const BasisUtils::F12BasisSet &f12BasisSet,
                                 const Eigen::MatrixXd &C);

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
