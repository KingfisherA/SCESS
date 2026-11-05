#pragma once
#include<Eigen/Core>
#include<unsupported/Eigen/CXX11/Tensor>


namespace SCESS {
    class F12BasisSet;
}

namespace SCESS::Integrals {
    struct MolInts;
}

namespace SCESS::Integrals::F12 {
    struct F12Integrals;
}

//
//
namespace SCESS::Integrals::Transforms {
    //
    //
    using TEITensor = Eigen::Tensor<double, 4>;
    /**
     * Transform the two-electron integrals (in chemist's notation) from the atomic orbital to the molecular orbital basis set
     * @param tei Two-electron integral tensor in atomic orbital basis (in chemist's notation)
     * @param C Wavefunction coefficients
     * @return Transformed TEI tensor (in chemist's notation)
     */
    TEITensor TransformTEI(const TEITensor &tei, const Eigen::MatrixXd &C, bool convertChemToPhys = false);
    //
    TEITensor TransformTEI(
        const TEITensor &tei,
        const Eigen::MatrixXd &A,
        const Eigen::MatrixXd &B,
        const Eigen::MatrixXd &C,
        const Eigen::MatrixXd &D,
        bool convertChemToPhys = false);
    //
    /**
     * Coefficients of the canonical OBS + CABS molecular orbitals in the AO basis of the OBS + CABS functions:
     * - The first nOBS columns are the Hartree-Fock orbitals (which only have OBS components)
     * - The OBS is projected out of the OBS + CABS functions and the remainder orthonormalised, giving the CABS
     * - The CABS is then canonicalised with the Fock operator
     * These are the requirements for the CABS to be used in orbital-sensitive F12 methods. Should the CABS functions
     * be linearly dependent on the OBS the redundant CABS orbitals are left as zero columns.
     * @param S AO overlap over the OBS + CABS functions (the OBS functions first)
     * @param F_AO AO Fock matrix over the OBS + CABS functions built from the OBS reference density
     * @param C_HF The Hartree Fock Wavefunction coefficients (nOBS x nOBS)
     */
    Eigen::MatrixXd BuildCABSCoefficients(const Eigen::MatrixXd &S,
                                          const Eigen::MatrixXd &F_AO,
                                          const Eigen::MatrixXd &C_HF,
                                          Eigen::Index nOBS,
                                          Eigen::Index nCABS);
    //
    /**
     * Produce a tensor with < i j || k l>=<i j | k l> - <i j | l k>
     * Automatically converts chemist's to physicists notation
     * @param tei Two-electron integral tensor (in chemist's notation)
     * @return Two-electron integral tensor with antisymmetric integrals (in physicists' notation)
     */
    TEITensor AntisymmetriseTEI(const TEITensor &tei);
    //
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
    //
    TEITensor ConvertPhysicistsToChemists(const TEITensor& tei);
    //
    //
}
