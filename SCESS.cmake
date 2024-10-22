add_executable("SCESS"
  "SCESS.APS"
  "SCESS.rc"
  "resource.h"
      "src/Application/Application.cpp"
      "src/Application/Application.h"
      "src/Chemistry/PeriodicTable.cpp"
      "src/Chemistry/PeriodicTable.h"
          "src/Chemistry/Quantum/Basis/BasisUtils.cpp"
          "src/Chemistry/Quantum/Basis/BasisUtils.h"
          "src/Chemistry/Quantum/Basis/Gaussian.cpp"
          "src/Chemistry/Quantum/Basis/Gaussian.h"
          "src/Chemistry/Quantum/CC/CCSD-DIIS.cpp"
          "src/Chemistry/Quantum/CC/CCSD.cpp"
          "src/Chemistry/Quantum/CC/CCSD_T.cpp"
          "src/Chemistry/Quantum/CC/CoupledCluster.h"
          "src/Chemistry/Quantum/CI/ConfigurationInteraction.h"
          "src/Chemistry/Quantum/CI/FullCI.cpp"
          "src/Chemistry/Quantum/DFT/DFT.cpp"
          "src/Chemistry/Quantum/DFT/DFT.h"
          "src/Chemistry/Quantum/DFT/XCFunctionals.cpp"
          "src/Chemistry/Quantum/DFT/XCFunctionals.h"
          "src/Chemistry/Quantum/HF/HartreeFock.h"
          "src/Chemistry/Quantum/HF/RestrictedHF.cpp"
          "src/Chemistry/Quantum/HF/UnrestrictedHF.cpp"
          "src/Chemistry/Quantum/Integrals/Gamma.cpp"
          "src/Chemistry/Quantum/Integrals/Gamma.h"
            "src/Chemistry/Quantum/Integrals/Grid/AtomicGrid.cpp"
            "src/Chemistry/Quantum/Integrals/Grid/AtomicGrid.h"
            "src/Chemistry/Quantum/Integrals/Grid/GridPoint.cpp"
            "src/Chemistry/Quantum/Integrals/Grid/GridPoint.h"
            "src/Chemistry/Quantum/Integrals/Grid/MolecularGrid.cpp"
            "src/Chemistry/Quantum/Integrals/Grid/MolecularGrid.h"
            "src/Chemistry/Quantum/Integrals/Grid/Quadrature.cpp"
            "src/Chemistry/Quantum/Integrals/Grid/Quadrature.h"
          "src/Chemistry/Quantum/Integrals/Integrals.cpp"
          "src/Chemistry/Quantum/Integrals/Integrals.h"
        "src/Chemistry/Quantum/Molecule.cpp"
        "src/Chemistry/Quantum/Molecule.h"
          "src/Chemistry/Quantum/QMC/CCQMC.cpp"
          "src/Chemistry/Quantum/QMC/FCIQMC.cpp"
          "src/Chemistry/Quantum/QMC/QMC.h"
          "src/Chemistry/Quantum/Tests/Tests.h"
        "src/Chemistry/Quantum/Tests.cpp"
        "src/Chemistry/Quantum/Tests.h"
      "src/GUI/CalculationLayer.cpp"
      "src/GUI/CalculationLayer.h"
      "src/GUI/ImSpinner.h"
      "src/GUI/ShaderSources.cpp"
      "src/GUI/ShaderSources.h"
      "src/Maths/Isosurfaces.cpp"
      "src/Maths/Isosurfaces.h"
      "src/Platform/Intrinsics.h"
      "src/Utilities/Utilities.cpp"
      "src/Utilities/Utilities.h"
    "src/pch.cpp"
    "src/pch.h"
)
if(CMAKE_BUILD_TYPE STREQUAL Debug)
  set_target_properties("SCESS" PROPERTIES
    OUTPUT_NAME "SCESS"
    ARCHIVE_OUTPUT_DIRECTORY "/media/charlie/Data/dev/C++/SCESS/bin/linux/Debug-x64"
    LIBRARY_OUTPUT_DIRECTORY "/media/charlie/Data/dev/C++/SCESS/bin/linux/Debug-x64"
    RUNTIME_OUTPUT_DIRECTORY "/media/charlie/Data/dev/C++/SCESS/bin/linux/Debug-x64"
  )
endif()
target_include_directories("SCESS" PRIVATE
  $<$<CONFIG:Debug>:/media/charlie/Data/dev/C++/SCESS/src>
  $<$<CONFIG:Debug>:/usr/include>
  $<$<CONFIG:Debug>:/usr/local/include>
  $<$<CONFIG:Debug>:/opt/etsf/include>
  $<$<CONFIG:Debug>:/usr/local/include/eigen3>
  $<$<CONFIG:Debug>:/media/charlie/Data/dev/C++/GraphicsEngine/src>
  $<$<CONFIG:Debug>:/media/charlie/Data/dev/Libraries/json-develop/include>
  $<$<CONFIG:Debug>:/media/charlie/Data/dev/Libraries/eigen-3.4.90>
  $<$<CONFIG:Debug>:/media/charlie/Data/dev/Libraries/glm>
  $<$<CONFIG:Debug>:/media/charlie/Data/dev/Libraries>
  $<$<CONFIG:Debug>:/media/charlie/Data/dev/Libraries/ImGui>
  $<$<CONFIG:Debug>:/media/charlie/Data/dev/Libraries/nativefiledialog/src/include>
  $<$<CONFIG:Debug>:/media/charlie/Data/dev/Libraries/cppitertools>
)
target_compile_definitions("SCESS" PRIVATE
  $<$<CONFIG:Debug>:OS_LINUX>
  $<$<CONFIG:Debug>:EIGEN_USE_BLAS>
  $<$<CONFIG:Debug>:EIGEN_USE_LAPACKE>
  $<$<CONFIG:Debug>:CONFIG_DEBUG>
  $<$<CONFIG:Debug>:_DEBUG>
)
target_link_directories("SCESS" PRIVATE
  $<$<CONFIG:Debug>:/media/charlie/Data/dev/Libraries/ImGui/bin/linux/Release-x64>
  $<$<CONFIG:Debug>:/media/charlie/Data/dev/Libraries/nativefiledialog/build/lib/Release/x64>
  $<$<CONFIG:Debug>:/media/charlie/Data/dev/C++/GraphicsEngine/bin/linux/Debug-x64/GraphicsEngine>
)
target_link_libraries("SCESS"
  $<$<CONFIG:Debug>:lapack>
  $<$<CONFIG:Debug>:lapacke>
  $<$<CONFIG:Debug>:blas>
  $<$<CONFIG:Debug>:GraphicsEngine>
  $<$<CONFIG:Debug>:glfw3>
  $<$<CONFIG:Debug>:ImGui>
  $<$<CONFIG:Debug>:libboost_filesystem.a>
  $<$<CONFIG:Debug>:GLEW.a>
  $<$<CONFIG:Debug>:GL>
  $<$<CONFIG:Debug>:m>
  $<$<CONFIG:Debug>:Xrandr>
  $<$<CONFIG:Debug>:Xi>
  $<$<CONFIG:Debug>:X11>
  $<$<CONFIG:Debug>:Xxf86vm>
  $<$<CONFIG:Debug>:pthread>
  $<$<CONFIG:Debug>:/opt/etsf/lib/libxc.a>
  $<$<CONFIG:Debug>:nfd>
  $<$<CONFIG:Debug>:gtk-3>
  $<$<CONFIG:Debug>:z>
  $<$<CONFIG:Debug>:gdk-3>
  $<$<CONFIG:Debug>:glib-2.0>
)
target_compile_options("SCESS" PRIVATE
  $<$<AND:$<CONFIG:Debug>,$<COMPILE_LANGUAGE:C>>:-m64>
  $<$<AND:$<CONFIG:Debug>,$<COMPILE_LANGUAGE:C>>:-g>
  $<$<AND:$<CONFIG:Debug>,$<COMPILE_LANGUAGE:CXX>>:-m64>
  $<$<AND:$<CONFIG:Debug>,$<COMPILE_LANGUAGE:CXX>>:-g>
  $<$<AND:$<CONFIG:Debug>,$<COMPILE_LANGUAGE:CXX>>:-std=c++20>
)
if(CMAKE_BUILD_TYPE STREQUAL Debug)
  set_target_properties("SCESS" PROPERTIES
    CXX_STANDARD 20
    CXX_STANDARD_REQUIRED YES
    CXX_EXTENSIONS NO
    POSITION_INDEPENDENT_CODE False
    INTERPROCEDURAL_OPTIMIZATION False
  )
endif()
if(CMAKE_BUILD_TYPE STREQUAL Release)
  set_target_properties("SCESS" PROPERTIES
    OUTPUT_NAME "SCESS"
    ARCHIVE_OUTPUT_DIRECTORY "/media/charlie/Data/dev/C++/SCESS/bin/linux/Release-x64"
    LIBRARY_OUTPUT_DIRECTORY "/media/charlie/Data/dev/C++/SCESS/bin/linux/Release-x64"
    RUNTIME_OUTPUT_DIRECTORY "/media/charlie/Data/dev/C++/SCESS/bin/linux/Release-x64"
  )
endif()
target_include_directories("SCESS" PRIVATE
  $<$<CONFIG:Release>:/media/charlie/Data/dev/C++/SCESS/src>
  $<$<CONFIG:Release>:/usr/include>
  $<$<CONFIG:Release>:/usr/local/include>
  $<$<CONFIG:Release>:/opt/etsf/include>
  $<$<CONFIG:Release>:/usr/local/include/eigen3>
  $<$<CONFIG:Release>:/media/charlie/Data/dev/C++/GraphicsEngine/src>
  $<$<CONFIG:Release>:/media/charlie/Data/dev/Libraries/json-develop/include>
  $<$<CONFIG:Release>:/media/charlie/Data/dev/Libraries/eigen-3.4.90>
  $<$<CONFIG:Release>:/media/charlie/Data/dev/Libraries/glm>
  $<$<CONFIG:Release>:/media/charlie/Data/dev/Libraries>
  $<$<CONFIG:Release>:/media/charlie/Data/dev/Libraries/ImGui>
  $<$<CONFIG:Release>:/media/charlie/Data/dev/Libraries/nativefiledialog/src/include>
  $<$<CONFIG:Release>:/media/charlie/Data/dev/Libraries/cppitertools>
)
target_compile_definitions("SCESS" PRIVATE
  $<$<CONFIG:Release>:OS_LINUX>
  $<$<CONFIG:Release>:EIGEN_USE_BLAS>
  $<$<CONFIG:Release>:EIGEN_USE_LAPACKE>
  $<$<CONFIG:Release>:CONFIG_RELEASE>
)
target_link_directories("SCESS" PRIVATE
  $<$<CONFIG:Release>:/media/charlie/Data/dev/Libraries/ImGui/bin/linux/Release-x64>
  $<$<CONFIG:Release>:/media/charlie/Data/dev/Libraries/nativefiledialog/build/lib/Release/x64>
  $<$<CONFIG:Release>:/media/charlie/Data/dev/C++/GraphicsEngine/bin/linux/Release-x64/GraphicsEngine>
)
target_link_libraries("SCESS"
  $<$<CONFIG:Release>:lapack>
  $<$<CONFIG:Release>:lapacke>
  $<$<CONFIG:Release>:blas>
  $<$<CONFIG:Release>:GraphicsEngine>
  $<$<CONFIG:Release>:glfw3>
  $<$<CONFIG:Release>:ImGui>
  $<$<CONFIG:Release>:libboost_filesystem.a>
  $<$<CONFIG:Release>:GLEW.a>
  $<$<CONFIG:Release>:GL>
  $<$<CONFIG:Release>:m>
  $<$<CONFIG:Release>:Xrandr>
  $<$<CONFIG:Release>:Xi>
  $<$<CONFIG:Release>:X11>
  $<$<CONFIG:Release>:Xxf86vm>
  $<$<CONFIG:Release>:pthread>
  $<$<CONFIG:Release>:/opt/etsf/lib/libxc.a>
  $<$<CONFIG:Release>:nfd>
  $<$<CONFIG:Release>:gtk-3>
  $<$<CONFIG:Release>:z>
  $<$<CONFIG:Release>:gdk-3>
  $<$<CONFIG:Release>:glib-2.0>
)
target_compile_options("SCESS" PRIVATE
  $<$<AND:$<CONFIG:Release>,$<COMPILE_LANGUAGE:C>>:-m64>
  $<$<AND:$<CONFIG:Release>,$<COMPILE_LANGUAGE:C>>:-O3>
  $<$<AND:$<CONFIG:Release>,$<COMPILE_LANGUAGE:C>>:-mavx2>
  $<$<AND:$<CONFIG:Release>,$<COMPILE_LANGUAGE:CXX>>:-m64>
  $<$<AND:$<CONFIG:Release>,$<COMPILE_LANGUAGE:CXX>>:-O3>
  $<$<AND:$<CONFIG:Release>,$<COMPILE_LANGUAGE:CXX>>:-mavx2>
  $<$<AND:$<CONFIG:Release>,$<COMPILE_LANGUAGE:CXX>>:-std=c++20>
)
if(CMAKE_BUILD_TYPE STREQUAL Release)
  set_target_properties("SCESS" PROPERTIES
    CXX_STANDARD 20
    CXX_STANDARD_REQUIRED YES
    CXX_EXTENSIONS NO
    POSITION_INDEPENDENT_CODE False
    INTERPROCEDURAL_OPTIMIZATION False
  )
endif()