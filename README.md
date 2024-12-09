
# SCESS
*Schrödinger's Cat's Electronic System Solver*

<div style="text-align:center">
<img src="./Assets/Icon/SCESS Small.png" width="150" height="150" /> 
</div>

# About
This is a basic implementation of different *ab-initio* quantum chemistry theories with a graphical editor to explore the solutions.

# Features
- Allows loading and displaying XYZ molecule files
- Can draw isosurface plots of the SCF orbital solutions
- Can calculate single-point energies with a variety of methods:
    - SCF Methods:
        - RHF (with DIIS)
        - UHF (with DIIS)
        - Kohn-Sham DFT
    - Correlation Methods:
        - MP2
        - Full-CI
        - CCSD (with DIIS)
        - CCSD(T)
        - CCQMC
        - FCIQMC


# Binaries
Release binaries are available for testing from the releases GitHub section.

# Build Instructions
*Currently only supports Linux build systems*
1. Clone the github repository *recursively* to get a copy of the submodules with the source: `git clone --recursive -j8 git@github.com:KingfisherA/SCESS.git`
2. Install the necessary dependency packages. This can be achieved by running `sudo Scripts\Get_Packages.sh`
3. Make both `Scripts\Build.sh` and `Libraries\Build_Libraries.sh` executable i.e. `chmod +x .\Scripts\Build.sh` and `chmod +x .\Libraries\Build_Libraries.sh`
4. Move into the `Scripts` folder; i.e. `cd .\Scripts`
5. Run `.\Build.sh`
6. Make a cup of tea while everything compiles ☕
7. Be happy when you have all your binaries 😊
8. The binaries are located in `..\bin\{OS}\{Release/Debug}-x64\`
9. Move into directory `cd ..\bin\linux\Release-x64`
10. You should be able to run `.\ViewerSCESS` straight away and the assets directory should be automatically found

# Dependencies

## Package Dependencies
- cmake
- clang++ (any version which supports C++20; normally compiled with clang >18)
- ninja
- LAPACK
- BLAS
- LAPACKE
- Boost
- libgtk-3

## Code Dependencies
These are automatically included and compiled in the project:
- Eigen
- cppitertools
- GLFW3
- GLEW
- Vulkan
- glm
- ImGui
- GraphicsEngine
- ImPlot
- ImSpinner
- Nlohmann JSON
- LibInt2
- LibXC
- nfd

# Warnings
- Currently this software is missing a lot of references to modified code it is using. This will be fixed soon.
- The basis sets assets do not have appropriate references either, however the basis sets are all derived from [basis set exchange](https://www.basissetexchange.org/). This will be updated soon.
- This software is unstable and may cause crashes. The coupled cluster theories are particularly bad at hogging large amounts of memory. Use at your own risk!
