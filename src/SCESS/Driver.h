#pragma once

namespace SCESS {
    class Molecule;
}

namespace SCESS::Driver {

    // -- Self-consistent field theories --

    //Kohn-Sham Density Functional Theory driver
    void DFT(const char *basisSetName, const Molecule& mol, const char *moleculeName = nullptr);
    //Restricted Hartree-Fock driver
    void RHF(const char *basisSetName, const Molecule& mol, const char *moleculeName = nullptr);
    //Unrestricted Hartree-Fock driver
    void UHF(const char *basisSetName, const Molecule& mol, const char *moleculeName = nullptr);


    // -- Moller-Plesset perturbation theory --

    //Closed-Shell Moller-Plesset second order driver
    void MP2(const char *basisSetName, const Molecule& mol, const char *moleculeName = nullptr);


    //-- Configuration interaction theories --

    //Full Configuration Interaction Driver
    void FCI(const char *basisSetName, const Molecule& mol, const char *moleculeName = nullptr);
    //Full Configuration Interaction Quantum Monte Carlo driver
    void FCIQMC(const char *basisSetName, const Molecule& mol, const char *moleculeName = nullptr);


    // -- Coupled Cluster Theories --

    //Coupled-cluster singles and doubles driver
    void CCSD(const char *basisSetName, const Molecule& mol, const char *moleculeName = nullptr);
    //Coupled cluster singles and doubles with perturbative triples driver
    void CCSD_T(const char *basisSetName, const Molecule& mol, const char *moleculeName = nullptr);
    //Coupled cluster quantum monte carlo driver
    void CCQMC(const char *basisSetName, const Molecule& mol, const char *moleculeName = nullptr);

    // -- F12 theories --

    //Moller-plesset second order F12 fixed amplitudes driver
    void MP2_F12_Fixed(const char *obsName,
                       const char *cabsName,
                       const Molecule &mol,
                       const char *moleculeName = nullptr);
    //Moller-plesset second order F12 optimise amplitudes driver
    void MP2_F12_Optimise(const char *obsName,
                          const char *cabsName,
                          const Molecule &mol,
                          const char *moleculeName = nullptr);
    //CCSD-F12b (3C(FIX) ansatz, CABS approximation) driver
    void CCSD_F12b(const char *obsName,
                   const char *cabsName,
                   const Molecule &mol,
                   const char *moleculeName = nullptr);


    // -- F12 QMC theories --

    //CCSD-QMC[MP2-F12] driver
    void CCSD_QMC_MP2_F12(const char *obsName,
                       const char *cabsName,
                       const Molecule &mol,
                       const char *moleculeName = nullptr);


}
