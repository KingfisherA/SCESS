#pragma once
#include<libint2/atom.h>
#include<libint2/shell.h>
#include "Integrals.h"

namespace SCESS {
    class F12BasisSet;
}

namespace SCESS::Integrals {
    struct MolInts;
}

//
//Evaluate integrals with f(12) correlation kernels
namespace SCESS::Integrals::F12 {
    //
    //
    using TEITensor = Eigen::Tensor<double, 4>;
    using LibIntCGF = std::vector<std::pair<double, double> >;
    using AtomListLI = std::vector<libint2::Atom>;
    //
    //
    struct F12Integrals {
    public:
        //f12 = ( ij | f(12) | kl)
        TEITensor f12;
        //f12Sq = ( ij | f(12)^2 | kl)
        TEITensor f12Sq;
        //f12DelSq = ( ij | |Del f(12)|^2 | kl)
        TEITensor f12DelSq;
        //f12Repulsion = ( ij | f(12)/r12 | kl)
        TEITensor f12Repulsion;
    };
    //
    //The integrals needed by the F12 methods in the canonical OBS + CABS MO basis (orbitals ordered occupied, OBS
    //virtual, CABS) in physicists' notation. Only blocks with an occupied pair and the OBS-only ERIs are kept: every
    //MP2-F12/3C and CCSD-F12b intermediate can be written in terms of these, and they are O(nOcc^2 N^2) rather
    //than O(N^4) for N = nOBS + nCABS.
    struct F12MOIntegrals {
        //Core Hamiltonian, Fock and exchange matrices over all N orbitals. K_pq = sum_i <ii|pq>.
        Eigen::MatrixXd h, F, K;
        //<pq|rs> over the OBS orbitals
        TEITensor eriOBS;
        //<ij|pq> for occupied i, j and all orbitals p, q, indexed (i, j, p, q)
        TEITensor eri, f12, f12Sq, f12DelSq, f12Repulsion;
        //<mx|pq> for occupied m, CABS x and OBS p, q, indexed (m, x - nOBS, p, q). Only built when requested
        //(it is needed by CCSD-QMC[MP2-F12] alone)
        TEITensor eri_oxpq;
    };
    //
    //Two-electron operators of the F12 methods. The geminal is the Gaussian-contracted Slater function
    //f12 = -exp(-gamma r12) / gamma
    enum class Operator2e { Coulomb, F12, F12Squared, F12Coulomb, DelF12Squared };
    //
    //Functions an AO index runs over: the OBS functions (which come first) or every OBS + CABS function
    enum class AOSpace { OBS, All };
    //
    /**
     * AO integrals (mu nu | lam sig) in chemists' notation over the first nShells shells, evaluated on every
     * hardware thread with the 8-fold permutational symmetry of the (real, symmetric) operators
     */
    TEITensor AOTwoElectronTensor(const std::vector<libint2::Shell> &shells, size_t nShells, Operator2e op,
                                  double gamma = 1.0);
    //
    //AO ERIs (mu nu | lam sig) over the OBS functions only
    TEITensor OBSTwoElectronRepulsionTensor(const F12BasisSet &basisSet);
    //
    /**
     * Half-transformed integrals H(p, nu, r, sig) = sum_{mu lam} C1(mu, p) C3(lam, r) (mu nu | lam sig), chemists'
     * notation, evaluated directly from the shell quartets on every hardware thread. The N^4 AO tensor is never
     * formed.
     * @param C1,C3 Coefficients with one row per OBS + CABS function. Shells whose rows are all zero are skipped,
     *              so occupied orbitals (which have no CABS component) only need the OBS shells.
     * @param space2,space4 The functions nu and sig run over, which sets dimensions 1 and 3 of the result
     */
    TEITensor HalfTransformTwoElectron(const F12BasisSet &basisSet, Operator2e op,
                                       const Eigen::MatrixXd &C1, AOSpace space2,
                                       const Eigen::MatrixXd &C3, AOSpace space4,
                                       double gamma = 1.0);
    //
    //Completes HalfTransformTwoElectron: (p q | r s) = sum_{nu sig} H(p, nu, r, sig) C2(nu, q) C4(sig, s)
    TEITensor FinishTransform(const TEITensor &halfTransformed, const Eigen::MatrixXd &C2, const Eigen::MatrixXd &C4);
    //
    //Coulomb matrix J_{mu nu} = sum_{lam sig} (mu nu | lam sig) D_{lam sig} over all OBS + CABS functions mu, nu
    //for a density D over the OBS functions
    Eigen::MatrixXd CoulombMatrix(const F12BasisSet &basisSet, const Eigen::MatrixXd &D_OBS);
    //
    //
    MolInts EvaluateMolecularIntegrals(const F12BasisSet &basisSet, const AtomListLI &atoms);
    MolInts ExtractHFIntegrals(const MolInts &f12MolInts, long startIdx, long hfBasisSetSize);
    //
    Eigen::MatrixXd CalculateOverlapMatrix(const F12BasisSet &basisSet);
    Eigen::MatrixXd CalculateNuclearMatrix(const F12BasisSet &basisSet, const AtomListLI &atoms);
    Eigen::MatrixXd CalculateKineticMatrix(const F12BasisSet &basisSet);
    TEITensor TwoElectronRepulsionTensor(const F12BasisSet &basisSet);
    //
    //Chemist's notation
    F12Integrals EvaluateF12Integrals(const F12BasisSet &basisSet, double gamma = 1.0);
    /**
     * Import an integral file from MRCC
     * @param fileName The file path of integral dump to read from
     * @param nbf Number of basis set functions in the OBS and CABS
     * @return The imported F12 integrals
     */
    void ImportF12IntegralsMRCC(const std::string &fileName, uint32_t nbf, F12Integrals &f12Ints, TEITensor &eri);
    //
    LibIntCGF CalculateSlaterContraction(double gamma = 1.0);
    LibIntCGF CalculateSlaterSquaredContraction(const LibIntCGF &slaterContraction);
    //
    /**
     * Evaluate the integral (ij | exp(-gamma * r12) | kl) for a Gaussian-contracted Slater type f(12)
     * @param basisSet Basis set to calculate matrix elements for
     * @param gamma The exponent of the f(12) correlation factor
     * @return A tensor of the required integrals
     */
    TEITensor TwoElectronCGFF12Tensor(const F12BasisSet &basisSet, const LibIntCGF &cgf);
    //
    /**
     * Evaluate the integral (ij | exp(-gamma * r12) * 1/r12 | kl) for a Gaussian-contracted Slater type f(12)
     * @param basisSet Basis set to calculate matrix elements for
     * @param gamma The exponent of the f(12) correlation factor
     * @return A tensor of the required integrals
     */
    TEITensor TwoElectronCGFF12RepulsionTensor(const F12BasisSet &basisSet, const LibIntCGF &cgf);
    //
    /**
     * Evaluate the integral (ij | |Del exp(-gamma * r12)|^2 | kl) for a Gaussian-contracted Slater type f(12)
     * @param basisSet
     * @param cgf
     * @return
     */
    TEITensor TwoElectronCGFDelF12SquaredTensor(const F12BasisSet &basisSet, const LibIntCGF &cgf);
    //
    //
}
