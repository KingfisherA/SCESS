#pragma once

//
//
namespace SCESS::Tests {
//
//
void ProfileMolInts();
void TestMolecularInts();
void TestDFT();
void TestCoupledCluster();
void TestCCSD_T();
void TestFullCI();
void TestFCIQMC();
void TestCCQMC();
void TestHelium();
void TestF12Integrals();
void TestMP2_F12_Optimise();
void TestMP2_F12_Fixed();
void TestCCQMC_MP2_F12();
// void TestUHF();
//
void TestTensorEinsum();
//
// Run the deterministic methods against reference energies without user
// interaction. Returns the number of failures.
int RunRegressionTests();
//
//
} // namespace SCESS::Tests
