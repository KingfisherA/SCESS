#pragma once

namespace SCESS {
    class Molecule;
}

namespace SCESS::Driver {
    //
    void MP2_F12_Fixed_H2O();
    void MP2_F12_Fixed_LiH();
    //
    void MP2_F12_Fixed(const char *obsName,
                       const char *cabsName,
                       const Molecule &mol,
                       const char *moleculeName = "");
    //
    void CCQMC_MP2_F12_H20();
    void CCQMC_MP2_F12_LiH();
    //
    void CCQMC_MP2_F12(const char *obsName, const char *cabsName,
                       const Molecule &mol, const char *moleculeName);
}
