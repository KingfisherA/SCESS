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
  "${SCESS_SRC_DIR}/Quantum/F12/F12.h"
  "${SCESS_SRC_DIR}/Quantum/F12/MP2_F12_Optimise.cpp"
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
  "${SCESS_SRC_DIR}/Quantum/PeriodicTable.cpp"
  "${SCESS_SRC_DIR}/Quantum/PeriodicTable.h"
  "${SCESS_SRC_DIR}/Quantum/QMC/CCQMC.cpp"
  "${SCESS_SRC_DIR}/Quantum/QMC/FCIQMC.cpp"
  "${SCESS_SRC_DIR}/Quantum/QMC/QMC.h"
  "${SCESS_SRC_DIR}/Quantum/F12/F12SinglesCorrection.cpp"
  "${SCESS_SRC_DIR}/Quantum/F12/CCQMC_F12.cpp"
  "${SCESS_SRC_DIR}/Tensor/TensorOperations.h"
  "${SCESS_SRC_DIR}/Tensor/TensorEinsum.h"
  "${SCESS_SRC_DIR}/Tensor/TensorSlice.h"
  "${SCESS_SRC_DIR}/Quantum/F12/MP2_F12_Fixed.cpp"
)

set_target_properties("QuantumSCESS" PROPERTIES PUBLIC_HEADER 
    "BasisSet.h"
    "Gaussian.h"
    "CoupledCluster.h"
    "ConfigurationInteraction.h"
    "DensityFunctional.h"
    "XCFunctional.h"
    "F12.h"
    "HartreeFock.h"
    "AtomicGrid.h"
    "GridPoint.h"
    "MolecularGrid.h"
    "Quadrature.h"
    "IntegralTransforms.h"
    "Integrals.h"
    "IntegralsF12.h"
    "IntegralsLibInt.h"
    "Molecule.h"
    "PeriodicTable.h"   
    "QMC.h"
    "Resources.h"
    "Assets.h"
    "Intrinsics.h"
    "Utilities.h"
)

find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
find_library(LIB_LAPACKE lapacke REQUIRED)

target_include_directories("QuantumSCESS" PRIVATE
  ${SCESS_SRC_DIR}
  ${GENG_INCLUDE_DIRS}
  ${IMGUI_INCLUDE_DIRS}
  ${SCESS_LIB_DIR}/cppitertools
  ${SCESS_LIB_DIR}/JSON/include
  ${SCESS_LIB_DIR}/Eigen
)

target_link_libraries("QuantumSCESS" PRIVATE
  LAPACK::LAPACK
  ${LIB_LAPACKE}
  BLAS::BLAS
  Boost::boost
  Boost::filesystem
  $<TARGET_NAME_IF_EXISTS:Boost::math>
  $<TARGET_NAME_IF_EXISTS:Boost::preprocessor>
  ${GENG_TOOLS_LIBRARIES}
  $<TARGET_NAME_IF_EXISTS:Libint2::int2>
  $<TARGET_NAME_IF_EXISTS:Libint2::libint2>
  Libxc::xc
)

target_compile_definitions("QuantumSCESS" PRIVATE
  $<$<PLATFORM_ID:Linux>:OS_LINUX>
  $<$<PLATFORM_ID:Windows>:OS_WINDOWS>
  EIGEN_USE_BLAS
  EIGEN_USE_LAPACKE
  LIBINT2_CPLUSPLUS_STD=2020
  $<$<CONFIG:DEBUG>:CONFIG_DEBUG>
  $<$<CONFIG:DEBUG>:_DEBUG>
  $<$<CONFIG:RELEASE>:CONFIG_RELEASE>
)

target_precompile_headers("QuantumSCESS" PRIVATE "${SCESS_SRC_DIR}/Quantum/QuantumPCH.h")

set(COMPILE_OPTIONS)
if((CMAKE_CXX_COMPILER_ID STREQUAL Clang) OR 
   (CMAKE_CXX_COMPILER_ID STREQUAL GNU))
    list(APPEND COMPILE_OPTIONS
      -m64
      -fPIC
      -std=c++20
    )
    if(CMAKE_BUILD_TYPE STREQUAL Debug)
      list(APPEND COMPILE_OPTIONS
        -g
        -O0
      )
    elseif(CMAKE_BUILD_TYPE STREQUAL Release)
      list(APPEND COMPILE_OPTIONS
        -O3
        -mavx2
      )
    endif()
elseif(CMAKE_CXX_COMPILER_ID STREQUAL MSVC)
  list(APPEND COMPILE_OPTIONS
      /std:c++20
    )
    if(CMAKE_BUILD_TYPE STREQUAL Debug)
      list(APPEND COMPILE_OPTIONS
        /Zi
        /Od
      )
    elseif(CMAKE_BUILD_TYPE STREQUAL Release)
      list(APPEND COMPILE_OPTIONS /O2)
    endif()
endif()
target_compile_options("QuantumSCESS" PRIVATE ${COMPILE_OPTIONS})

set_target_properties("QuantumSCESS" PROPERTIES
    CXX_STANDARD 20
    CXX_STANDARD_REQUIRED YES
    CXX_EXTENSIONS NO
    POSITION_INDEPENDENT_CODE True
    INTERPROCEDURAL_OPTIMIZATION False
    MSVC_RUNTIME_LIBRARY MultiThreaded$<$<CONFIG:Debug>:Debug>DLL
    OUTPUT_NAME "QuantumSCESS")

if(CMAKE_CONFIGURATION_TYPES)
  set(DEBUG_BIN_DIR "${SCESS_BIN_DIR}")
  set(RELEASE_BIN_DIR "${SCESS_BIN_DIR}")
else()
  set(DEBUG_BIN_DIR "${SCESS_BIN_DIR}/Debug")
  set(RELEASE_BIN_DIR "${SCESS_BIN_DIR}/Release")
endif()

if(CMAKE_BUILD_TYPE STREQUAL Debug)
  set_target_properties("QuantumSCESS" PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY "${DEBUG_BIN_DIR}"
    LIBRARY_OUTPUT_DIRECTORY "${DEBUG_BIN_DIR}"
    RUNTIME_OUTPUT_DIRECTORY "${DEBUG_BIN_DIR}"
  )
elseif(CMAKE_BUILD_TYPE STREQUAL Release)
  set_target_properties("QuantumSCESS" PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY "${RELEASE_BIN_DIR}"
    LIBRARY_OUTPUT_DIRECTORY "${RELEASE_BIN_DIR}"
    RUNTIME_OUTPUT_DIRECTORY "${RELEASE_BIN_DIR}"
  )
endif()
