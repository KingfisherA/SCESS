# QuantumSCESS: shared library with the quantum chemistry code, built into build/QuantumSCESS

add_library(QuantumSCESS SHARED
  "${SCESS_SRC_DIR}/Resources/AssetFolder.cpp"
  "${SCESS_SRC_DIR}/Resources/AssetFolder.h"
  "${SCESS_SRC_DIR}/Quantum/Basis/BasisSet.cpp"
  "${SCESS_SRC_DIR}/Quantum/Basis/BasisSet.h"
  "${SCESS_SRC_DIR}/Quantum/Basis/CGF.cpp"
  "${SCESS_SRC_DIR}/Quantum/Basis/CGF.h"
  "${SCESS_SRC_DIR}/Quantum/CC/CCSD-DIIS.cpp"
  "${SCESS_SRC_DIR}/Quantum/CC/CCSD.cpp"
  "${SCESS_SRC_DIR}/Quantum/CC/CCSD_T.cpp"
  "${SCESS_SRC_DIR}/Quantum/CC/RCCSD.cpp"
  "${SCESS_SRC_DIR}/Quantum/CC/RCCSD_T.cpp"
  "${SCESS_SRC_DIR}/Quantum/CC/CC.h"
  "${SCESS_SRC_DIR}/Quantum/CI/CI.h"
  "${SCESS_SRC_DIR}/Quantum/DFT/DFT.cpp"
  "${SCESS_SRC_DIR}/Quantum/DFT/DFT.h"
  "${SCESS_SRC_DIR}/Quantum/DFT/XC.cpp"
  "${SCESS_SRC_DIR}/Quantum/DFT/XC.h"
  "${SCESS_SRC_DIR}/Quantum/F12/F12.h"
  "${SCESS_SRC_DIR}/Quantum/F12/MP2_F12_Optimise.cpp"
  "${SCESS_SRC_DIR}/Quantum/F12/CCSD_F12b.cpp"
  "${SCESS_SRC_DIR}/Quantum/F12/F12Setup.cpp"
  "${SCESS_SRC_DIR}/Quantum/HF/HF.h"
  "${SCESS_SRC_DIR}/Quantum/HF/DIIS.h"
  "${SCESS_SRC_DIR}/Quantum/HF/RHF.cpp"
  "${SCESS_SRC_DIR}/Quantum/HF/UHF.cpp"
  "${SCESS_SRC_DIR}/Quantum/Integrals/Grid/AtomicGrid.cpp"
  "${SCESS_SRC_DIR}/Quantum/Integrals/Grid/AtomicGrid.h"
  "${SCESS_SRC_DIR}/Quantum/Integrals/Grid/GridPoint.cpp"
  "${SCESS_SRC_DIR}/Quantum/Integrals/Grid/GridPoint.h"
  "${SCESS_SRC_DIR}/Quantum/Integrals/Grid/MolGrid.cpp"
  "${SCESS_SRC_DIR}/Quantum/Integrals/Grid/MolGrid.h"
  "${SCESS_SRC_DIR}/Quantum/Integrals/Grid/Quadrature.cpp"
  "${SCESS_SRC_DIR}/Quantum/Integrals/Grid/Quadrature.h"
  "${SCESS_SRC_DIR}/Quantum/Integrals/IntegralBlocks.cpp"
  "${SCESS_SRC_DIR}/Quantum/Integrals/IntegralImporter.cpp"
  "${SCESS_SRC_DIR}/Quantum/Integrals/Integrals.cpp"
  "${SCESS_SRC_DIR}/Quantum/Integrals/Integrals.h"
  "${SCESS_SRC_DIR}/Quantum/Integrals/IntegralsF12.h"
  "${SCESS_SRC_DIR}/Quantum/Integrals/IntegralsLibInt.cpp"
  "${SCESS_SRC_DIR}/Quantum/Integrals/IntegralsLibInt.h"
  "${SCESS_SRC_DIR}/Quantum/Integrals/IntegralTransforms.cpp"
  "${SCESS_SRC_DIR}/Quantum/Integrals/IntegralTransforms.h"
  "${SCESS_SRC_DIR}/Quantum/Molecule.cpp"
  "${SCESS_SRC_DIR}/Quantum/Molecule.h"
  "${SCESS_SRC_DIR}/Quantum/QMC/CCQMC.cpp"
  "${SCESS_SRC_DIR}/Quantum/QMC/FCIQMC.cpp"
  "${SCESS_SRC_DIR}/Quantum/QMC/QMC.h"
  "${SCESS_SRC_DIR}/Quantum/F12/F12SinglesCorrection.cpp"
  "${SCESS_SRC_DIR}/Tensor/TensorOperations.h"
  "${SCESS_SRC_DIR}/Tensor/TensorEinsum.h"
  "${SCESS_SRC_DIR}/Tensor/TensorParallel.h"
  "${SCESS_SRC_DIR}/Tensor/TensorSlice.h"
  "${SCESS_SRC_DIR}/Quantum/F12/MP2_F12_Fixed.cpp"
  "${SCESS_SRC_DIR}/Quantum/F12_QMC/F12_QMC.h"
  "${SCESS_SRC_DIR}/Quantum/F12_QMC/CCSD_QMC_MP2_F12.cpp"
  "${SCESS_SRC_DIR}/Quantum/CI/CI.cpp"
)
scess_target_defaults(QuantumSCESS)

target_include_directories(QuantumSCESS PUBLIC
  ${SCESS_SRC_DIR}
  ${DAHLIA_INCLUDE_DIR}
  ${SCESS_LIB_DIR}/cppitertools
  ${SCESS_LIB_DIR}/JSON/include
  ${SCESS_LIB_DIR}/Eigen
)

target_compile_definitions(QuantumSCESS
  PUBLIC
    # Fallback assets location when no settings.json or bundled Assets folder is found
    SCESS_SOURCE_ASSETS_DIR="${SCESS_ASSETS_DIR}"
    # Enables Eigen's ThreadPoolDevice for the multithreaded tensor contractions (Tensor/TensorParallel.h). Must be
    # defined identically in every translation unit which includes the Eigen Tensor module
    EIGEN_USE_THREADS
  PRIVATE
    EIGEN_USE_BLAS
    EIGEN_USE_LAPACKE
)

target_link_libraries(QuantumSCESS
  PUBLIC
    Boost::filesystem
    Libint2::libint2_cxx
    Libxc::xc
    fmt::fmt
  PRIVATE
    DahliaTools
    ${LIB_LAPACKE}
    ${LAPACK_LIBRARIES}
    $<$<PLATFORM_ID:Darwin>:$<LINK_LIBRARY:FRAMEWORK,Accelerate>>
)
