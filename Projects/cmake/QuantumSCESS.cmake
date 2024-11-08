add_library("QuantumSCESS" SHARED
      "src/Quantum/Basis/BasisUtils.cpp"
      "src/Quantum/Basis/BasisUtils.h"
      "src/Quantum/Basis/Gaussian.cpp"
      "src/Quantum/Basis/Gaussian.h"
      "src/Quantum/CC/CCSD-DIIS.cpp"
      "src/Quantum/CC/CCSD.cpp"
      "src/Quantum/CC/CCSD_T.cpp"
      "src/Quantum/CC/CoupledCluster.h"
      "src/Quantum/CI/ConfigurationInteraction.h"
      "src/Quantum/CI/FullCI.cpp"
      "src/Quantum/DFT/DensityFunctional.cpp"
      "src/Quantum/DFT/DensityFunctional.h"
      "src/Quantum/DFT/XCFunctionals.cpp"
      "src/Quantum/DFT/XCFunctionals.h"
      "src/Quantum/HF/HartreeFock.h"
      "src/Quantum/HF/RestrictedHF.cpp"
      "src/Quantum/HF/UnrestrictedHF.cpp"
      "src/Quantum/Integrals/BoysFunc.cpp"
      "src/Quantum/Integrals/BoysFunc.h"
        "src/Quantum/Integrals/Grid/AtomicGrid.cpp"
        "src/Quantum/Integrals/Grid/AtomicGrid.h"
        "src/Quantum/Integrals/Grid/GridPoint.cpp"
        "src/Quantum/Integrals/Grid/GridPoint.h"
        "src/Quantum/Integrals/Grid/MolecularGrid.cpp"
        "src/Quantum/Integrals/Grid/MolecularGrid.h"
        "src/Quantum/Integrals/Grid/Quadrature.cpp"
        "src/Quantum/Integrals/Grid/Quadrature.h"
      "src/Quantum/Integrals/Integrals.cpp"
      "src/Quantum/Integrals/Integrals.h"
      "src/Quantum/Integrals/IntegralsF12.cpp"
      "src/Quantum/Integrals/IntegralsF12.h"
    "src/Quantum/Molecule.cpp"
    "src/Quantum/Molecule.h"
    "src/Quantum/PeriodicTable.cpp"
    "src/Quantum/PeriodicTable.h"
      "src/Quantum/QMC/CCQMC.cpp"
      "src/Quantum/QMC/FCIQMC.cpp"
      "src/Quantum/QMC/QMC.h"
    "src/Resources/AssetFolder.cpp"
    "src/Resources/Fonts.cpp"
    "src/Resources/Resources.h"
    "src/Resources/ShaderSources.cpp"
    "src/Utilities/Assets.h"
    "src/Utilities/Intrinsics.h"
    "src/Utilities/Utilities.cpp"
    "src/Utilities/Utilities.h"
)
if(CMAKE_BUILD_TYPE STREQUAL Debug)
  set_target_properties("QuantumSCESS" PROPERTIES
    OUTPUT_NAME "QuantumSCESS"
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin/linux/Debug-x64"
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin/linux/Debug-x64"
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin/linux/Debug-x64"
  )
endif()
target_include_directories("QuantumSCESS" PRIVATE
  $<$<CONFIG:Debug>:${CMAKE_CURRENT_SOURCE_DIR}/src>
  $<$<CONFIG:Debug>:${CMAKE_CURRENT_SOURCE_DIR}/Libraries/GraphicsEngine-1.0.0/include>
  $<$<CONFIG:Debug>:${CMAKE_CURRENT_SOURCE_DIR}/Libraries/glm>
  $<$<CONFIG:Debug>:${CMAKE_CURRENT_SOURCE_DIR}/Libraries/cppitertools>
  $<$<CONFIG:Debug>:${CMAKE_CURRENT_SOURCE_DIR}/Libraries/JSON/include>
  $<$<CONFIG:Debug>:${CMAKE_CURRENT_SOURCE_DIR}/Libraries/Eigen>
  $<$<CONFIG:Debug>:${CMAKE_CURRENT_SOURCE_DIR}/Libraries/LibInt2/include>
  $<$<CONFIG:Debug>:${CMAKE_CURRENT_SOURCE_DIR}/Libraries/LibXC/lib/include>
)
target_compile_definitions("QuantumSCESS" PRIVATE
  $<$<CONFIG:Debug>:OS_LINUX>
  $<$<CONFIG:Debug>:EIGEN_USE_BLAS>
  $<$<CONFIG:Debug>:EIGEN_USE_LAPACKE>
  $<$<CONFIG:Debug>:CONFIG_DEBUG>
  $<$<CONFIG:Debug>:_DEBUG>
)
target_link_directories("QuantumSCESS" PRIVATE
  $<$<CONFIG:Debug>:${CMAKE_CURRENT_SOURCE_DIR}/Libraries/GraphicsEngine-1.0.0>
  $<$<CONFIG:Debug>:${CMAKE_CURRENT_SOURCE_DIR}/Libraries/LibInt2>
  $<$<CONFIG:Debug>:${CMAKE_CURRENT_SOURCE_DIR}/Libraries/LibXC/lib>
)
target_link_libraries("QuantumSCESS"
  $<$<CONFIG:Debug>:lapack>
  $<$<CONFIG:Debug>:lapacke>
  $<$<CONFIG:Debug>:blas>
  $<$<CONFIG:Debug>:boost_filesystem>
  $<$<CONFIG:Debug>:GraphicsEngine>
  $<$<CONFIG:Debug>:libint2.a>
  $<$<CONFIG:Debug>:libxc.a>
)
target_compile_options("QuantumSCESS" PRIVATE
  $<$<AND:$<CONFIG:Debug>,$<COMPILE_LANGUAGE:C>>:-m64>
  $<$<AND:$<CONFIG:Debug>,$<COMPILE_LANGUAGE:C>>:-O0>
  $<$<AND:$<CONFIG:Debug>,$<COMPILE_LANGUAGE:C>>:-fPIC>
  $<$<AND:$<CONFIG:Debug>,$<COMPILE_LANGUAGE:C>>:-g>
  $<$<AND:$<CONFIG:Debug>,$<COMPILE_LANGUAGE:CXX>>:-m64>
  $<$<AND:$<CONFIG:Debug>,$<COMPILE_LANGUAGE:CXX>>:-O0>
  $<$<AND:$<CONFIG:Debug>,$<COMPILE_LANGUAGE:CXX>>:-fPIC>
  $<$<AND:$<CONFIG:Debug>,$<COMPILE_LANGUAGE:CXX>>:-g>
  $<$<AND:$<CONFIG:Debug>,$<COMPILE_LANGUAGE:CXX>>:-std=c++20>
)
if(CMAKE_BUILD_TYPE STREQUAL Debug)
  set_target_properties("QuantumSCESS" PROPERTIES
    CXX_STANDARD 20
    CXX_STANDARD_REQUIRED YES
    CXX_EXTENSIONS NO
    POSITION_INDEPENDENT_CODE True
    INTERPROCEDURAL_OPTIMIZATION False
  )
endif()
if(CMAKE_BUILD_TYPE STREQUAL Debug)
target_precompile_headers("QuantumSCESS" PUBLIC ./src/pch.h)
endif()
if(CMAKE_BUILD_TYPE STREQUAL Release)
  set_target_properties("QuantumSCESS" PROPERTIES
    OUTPUT_NAME "QuantumSCESS"
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin/linux/Release-x64"
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin/linux/Release-x64"
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin/linux/Release-x64"
  )
endif()
target_include_directories("QuantumSCESS" PRIVATE
  $<$<CONFIG:Release>:${CMAKE_CURRENT_SOURCE_DIR}/src>
  $<$<CONFIG:Release>:${CMAKE_CURRENT_SOURCE_DIR}/Libraries/GraphicsEngine-1.0.0/include>
  $<$<CONFIG:Release>:${CMAKE_CURRENT_SOURCE_DIR}/Libraries/glm>
  $<$<CONFIG:Release>:${CMAKE_CURRENT_SOURCE_DIR}/Libraries/cppitertools>
  $<$<CONFIG:Release>:${CMAKE_CURRENT_SOURCE_DIR}/Libraries/JSON/include>
  $<$<CONFIG:Release>:${CMAKE_CURRENT_SOURCE_DIR}/Libraries/Eigen>
  $<$<CONFIG:Release>:${CMAKE_CURRENT_SOURCE_DIR}/Libraries/LibInt2/include>
  $<$<CONFIG:Release>:${CMAKE_CURRENT_SOURCE_DIR}/Libraries/LibXC/lib/include>
)
target_compile_definitions("QuantumSCESS" PRIVATE
  $<$<CONFIG:Release>:OS_LINUX>
  $<$<CONFIG:Release>:EIGEN_USE_BLAS>
  $<$<CONFIG:Release>:EIGEN_USE_LAPACKE>
  $<$<CONFIG:Release>:CONFIG_RELEASE>
)
target_link_directories("QuantumSCESS" PRIVATE
  $<$<CONFIG:Release>:${CMAKE_CURRENT_SOURCE_DIR}/Libraries/GraphicsEngine-1.0.0>
  $<$<CONFIG:Release>:${CMAKE_CURRENT_SOURCE_DIR}/Libraries/LibInt2>
  $<$<CONFIG:Release>:${CMAKE_CURRENT_SOURCE_DIR}/Libraries/LibXC/lib>
)
target_link_libraries("QuantumSCESS"
  $<$<CONFIG:Release>:lapack>
  $<$<CONFIG:Release>:lapacke>
  $<$<CONFIG:Release>:blas>
  $<$<CONFIG:Release>:boost_filesystem>
  $<$<CONFIG:Release>:GraphicsEngine>
  $<$<CONFIG:Release>:libint2.a>
  $<$<CONFIG:Release>:libxc.a>
)
target_compile_options("QuantumSCESS" PRIVATE
  $<$<AND:$<CONFIG:Release>,$<COMPILE_LANGUAGE:C>>:-m64>
  $<$<AND:$<CONFIG:Release>,$<COMPILE_LANGUAGE:C>>:-O3>
  $<$<AND:$<CONFIG:Release>,$<COMPILE_LANGUAGE:C>>:-fPIC>
  $<$<AND:$<CONFIG:Release>,$<COMPILE_LANGUAGE:C>>:-mavx2>
  $<$<AND:$<CONFIG:Release>,$<COMPILE_LANGUAGE:CXX>>:-m64>
  $<$<AND:$<CONFIG:Release>,$<COMPILE_LANGUAGE:CXX>>:-O3>
  $<$<AND:$<CONFIG:Release>,$<COMPILE_LANGUAGE:CXX>>:-fPIC>
  $<$<AND:$<CONFIG:Release>,$<COMPILE_LANGUAGE:CXX>>:-mavx2>
  $<$<AND:$<CONFIG:Release>,$<COMPILE_LANGUAGE:CXX>>:-std=c++20>
)
if(CMAKE_BUILD_TYPE STREQUAL Release)
  set_target_properties("QuantumSCESS" PROPERTIES
    CXX_STANDARD 20
    CXX_STANDARD_REQUIRED YES
    CXX_EXTENSIONS NO
    POSITION_INDEPENDENT_CODE True
    INTERPROCEDURAL_OPTIMIZATION False
  )
endif()
if(CMAKE_BUILD_TYPE STREQUAL Release)
target_precompile_headers("QuantumSCESS" PUBLIC ./src/pch.h)
endif()