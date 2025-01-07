# Project defining the shared library for the quantum calculation code for SCESS 

add_library("QuantumSCESS" SHARED
  "${SCESS_SRC_DIR}/Quantum/Basis/BasisSet.cpp"
  "${SCESS_SRC_DIR}/Quantum/Basis/BasisSet.h"
  "${SCESS_SRC_DIR}/Quantum/Basis/Gaussian.cpp"
  "${SCESS_SRC_DIR}/Quantum/Basis/Gaussian.h"
  "${SCESS_SRC_DIR}/Quantum/CC/CCSD-DIIS.cpp"
  "${SCESS_SRC_DIR}/Quantum/CC/CCSD.cpp"
  "${SCESS_SRC_DIR}/Quantum/CC/CCSD_T.cpp"
  "${SCESS_SRC_DIR}/Quantum/CC/CoupledCluster.h"
  "${SCESS_SRC_DIR}/Quantum/CI/ConfigurationInteraction.h"
  "${SCESS_SRC_DIR}/Quantum/CI/FullCI.cpp"
  "${SCESS_SRC_DIR}/Quantum/DFT/DensityFunctional.cpp"
  "${SCESS_SRC_DIR}/Quantum/DFT/DensityFunctional.h"
  "${SCESS_SRC_DIR}/Quantum/DFT/XCFunctionals.cpp"
  "${SCESS_SRC_DIR}/Quantum/DFT/XCFunctionals.h"
  "${SCESS_SRC_DIR}/Quantum/F12/CCSD_F12.cpp"
  "${SCESS_SRC_DIR}/Quantum/F12/CCSD_F12.h"
  "${SCESS_SRC_DIR}/Quantum/HF/HartreeFock.h"
  "${SCESS_SRC_DIR}/Quantum/HF/RestrictedHF.cpp"
  "${SCESS_SRC_DIR}/Quantum/HF/UnrestrictedHF.cpp"
  "${SCESS_SRC_DIR}/Quantum/Integrals/Grid/AtomicGrid.cpp"
  "${SCESS_SRC_DIR}/Quantum/Integrals/Grid/AtomicGrid.h"
  "${SCESS_SRC_DIR}/Quantum/Integrals/Grid/GridPoint.cpp"
  "${SCESS_SRC_DIR}/Quantum/Integrals/Grid/GridPoint.h"
  "${SCESS_SRC_DIR}/Quantum/Integrals/Grid/MolecularGrid.cpp"
  "${SCESS_SRC_DIR}/Quantum/Integrals/Grid/MolecularGrid.h"
  "${SCESS_SRC_DIR}/Quantum/Integrals/Grid/Quadrature.cpp"
  "${SCESS_SRC_DIR}/Quantum/Integrals/Grid/Quadrature.h"
  "${SCESS_SRC_DIR}/Quantum/Integrals/IntegralTransforms.cpp"
  "${SCESS_SRC_DIR}/Quantum/Integrals/IntegralTransforms.h"
  "${SCESS_SRC_DIR}/Quantum/Integrals/Integrals.cpp"
  "${SCESS_SRC_DIR}/Quantum/Integrals/Integrals.h"
  "${SCESS_SRC_DIR}/Quantum/Integrals/IntegralsF12.h"
  "${SCESS_SRC_DIR}/Quantum/Integrals/IntegralsLibInt.cpp"
  "${SCESS_SRC_DIR}/Quantum/Integrals/IntegralsLibInt.h"
  "${SCESS_SRC_DIR}/Quantum/Molecule.cpp"
  "${SCESS_SRC_DIR}/Quantum/Molecule.h"
  "${SCESS_SRC_DIR}/Quantum/PeriodicTable.cpp"
  "${SCESS_SRC_DIR}/Quantum/PeriodicTable.h"
  "${SCESS_SRC_DIR}/Quantum/QMC/CCQMC.cpp"
  "${SCESS_SRC_DIR}/Quantum/QMC/FCIQMC.cpp"
  "${SCESS_SRC_DIR}/Quantum/QMC/QMC.h"
  "${SCESS_SRC_DIR}/Resources/AssetFolder.cpp"
  "${SCESS_SRC_DIR}/Resources/Resources.h"
  "${SCESS_SRC_DIR}/Resources/ShaderSources.cpp"
  "${SCESS_SRC_DIR}/Utilities/Assets.h"
  "${SCESS_SRC_DIR}/Utilities/Intrinsics.h"
  "${SCESS_SRC_DIR}/Utilities/Utilities.cpp"
  "${SCESS_SRC_DIR}/Utilities/Utilities.h"
)

find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
find_library(LIB_LAPACKE lapacke REQUIRED)

target_include_directories("QuantumSCESS" PRIVATE
  ${SCESS_SRC_DIR}
  ${SCESS_LIB_DIR}/GraphicsEngine-1.0.0/include
  ${SCESS_LIB_DIR}/cppitertools
  ${SCESS_LIB_DIR}/JSON/include
  ${SCESS_LIB_DIR}/Eigen
)

target_link_libraries("QuantumSCESS" PRIVATE
  LAPACK::LAPACK
  ${LIB_LAPACKE}
  BLAS::BLAS
  Boost::filesystem
  ${LIB_GraphicsEngine}
  Libint2::libint2
  Libxc::xc
)
target_compile_definitions("QuantumSCESS" PRIVATE
  OS_LINUX
  EIGEN_USE_BLAS
  EIGEN_USE_LAPACKE
  $<$<CONFIG:DEBUG>:CONFIG_DEBUG>
  $<$<CONFIG:DEBUG>:_DEBUG>
  $<$<CONFIG:RELEASE>:CONFIG_RELEASE>
)
target_precompile_headers("QuantumSCESS" PUBLIC ./src/pch.h)
target_compile_options("QuantumSCESS" PRIVATE
    -m64
    -fPIC
    -std=c++20
)

if(CMAKE_BUILD_TYPE STREQUAL Debug)
  set_target_properties("QuantumSCESS" PROPERTIES
    CXX_STANDARD 20
    CXX_STANDARD_REQUIRED YES
    CXX_EXTENSIONS NO
    POSITION_INDEPENDENT_CODE True
    INTERPROCEDURAL_OPTIMIZATION False
    OUTPUT_NAME "QuantumSCESS"
    ARCHIVE_OUTPUT_DIRECTORY "${SCESS_BIN_DIR}/Debug-x64"
    LIBRARY_OUTPUT_DIRECTORY "${SCESS_BIN_DIR}/Debug-x64"
    RUNTIME_OUTPUT_DIRECTORY "${SCESS_BIN_DIR}/Debug-x64"
  )

  target_compile_options("QuantumSCESS" PRIVATE
    -g
    -O0
  )
endif()


if(CMAKE_BUILD_TYPE STREQUAL Release)
  set_target_properties("QuantumSCESS" PROPERTIES
    CXX_STANDARD 20
    CXX_STANDARD_REQUIRED YES
    CXX_EXTENSIONS NO
    POSITION_INDEPENDENT_CODE True
    INTERPROCEDURAL_OPTIMIZATION False
    OUTPUT_NAME "QuantumSCESS"
    ARCHIVE_OUTPUT_DIRECTORY "${SCESS_BIN_DIR}/Release-x64"
    LIBRARY_OUTPUT_DIRECTORY "${SCESS_BIN_DIR}/Release-x64"
    RUNTIME_OUTPUT_DIRECTORY "${SCESS_BIN_DIR}/Release-x64"
  )

  target_compile_options("QuantumSCESS" PRIVATE
      -mavx2
      -O2
  )
endif()

