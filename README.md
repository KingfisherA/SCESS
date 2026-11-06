<p align="center">
  <img src="Assets/Icon/SCESS%20Icon.svg" alt="SCESS icon" width="200" height="200" />
</p>

<h1 align="center">SCESS</h1>
<p align="center"><em>Schrödinger's Cat's Electronic System Solver</em></p>

# About
SCESS is an implementation of a range of *ab initio* quantum chemistry methods, with a graphical
viewer to explore the solutions. It is made up of:

- **QuantumSCESS**: a shared library containing the quantum chemistry code
- **SCESS**: a command line driver for running single-point calculations
- **ViewerSCESS**: a graphical viewer built on the [Dahlia](https://github.com/KingfisherA/GraphicsEngine)
  graphics engine (Metal by default on macOS, OpenGL 4.1 elsewhere; also Vulkan 1.3)

# Features
- Loads and displays XYZ molecule files
- Draws isosurface plots of the SCF orbital solutions
- Calculates single-point energies with a variety of methods:
    - SCF methods:
        - RHF (with DIIS)
        - UHF (with DIIS)
        - Kohn-Sham DFT
    - Correlation methods:
        - MP2 and MP2-F12 (fixed and optimised amplitudes)
        - Full-CI and truncated CI (Davidson diagonalisation; string-driven Hamiltonian products for full CI)
        - CCSD (closed-shell, spatial orbitals, with DIIS)
        - CCSD(T) (closed-shell triples over the unique virtual triples)
        - CCSD-F12b (3C(FIX) ansatz with the CABS approximation, closed-shell spatial-orbital CCSD)
        - CCQMC and FCIQMC
        - CCSD-QMC[MP2-F12]

# Building

## Prerequisites
SCESS builds the Dahlia engine from source as part of its own build, so clone both repositories
side by side (SCESS looks for `../Dahlia` by default):

```sh
git clone --recursive git@github.com:KingfisherA/SCESS.git
git clone --recursive git@github.com:KingfisherA/GraphicsEngine.git Dahlia
```

If you cloned without `--recursive`, run `git submodule update --init --recursive` in each repository.

You will also need:

- CMake 3.24 or newer, Ninja and a C++20 compiler (Clang 18+ or Apple Clang)
- Boost, fmt, LAPACK/LAPACKE, [Libxc](https://libxc.gitlab.io) and
  [Libint2](https://github.com/evaleev/libint) (built with C++ support and installed with its CMake package)
- Dahlia's dependencies: `glslc` (from shaderc or the Vulkan SDK), GLFW, GLEW, GLM and FreeType, plus
  optionally the Vulkan SDK and [metal-cpp](https://developer.apple.com/metal/cpp/) (see Dahlia's README)
- Linux only: GTK 3 (for the file dialogs)

On macOS, with Homebrew:

```sh
brew install cmake ninja boost fmt lapack libxc shaderc glfw glew glm freetype
```

On Debian/Ubuntu, `sudo Scripts/Get_Packages.sh` installs everything except Libint2.

## Configure and build
SCESS always builds into `./build`:

```sh
Scripts/Build.sh Release        # or Debug
```

This runs CMake through `cmake.sh`, which adds Homebrew's packages to the search path and uses
ccache when it is installed. The equivalent CMake commands are:

```sh
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build
```

In VS Code, the CMake Tools extension uses the same `build` folder (see `.vscode/settings.json`).
If Dahlia is somewhere other than `../Dahlia`, pass `-DSCESS_DAHLIA_DIR=/path/to/Dahlia`.

Each target is built into its own folder:

| Output                              | Contents                                    |
|-------------------------------------|---------------------------------------------|
| `build/QuantumSCESS/`               | `libQuantumSCESS`                           |
| `build/SCESS/SCESS`                 | the command line driver                     |
| `build/ViewerSCESS/ViewerSCESS`     | the graphical viewer                        |
| `build/ViewerSCESS/ViewerSCESS.app` | the macOS app bundle (see below)            |
| `build/Dahlia/`                     | `libDahlia`, built from `SCESS_DAHLIA_DIR`  |

## macOS app bundle
```sh
Scripts/Build.sh Release --target ViewerSCESSApp
```

This builds `build/ViewerSCESS/ViewerSCESS.app`, which has the SCESS icon and contains the `Assets`
folder and every non-system library it needs (including libDahlia and Homebrew's libraries), so it
can be moved anywhere and started from Finder. The bundle is ad-hoc signed.

# Running
```sh
./build/ViewerSCESS/ViewerSCESS              # Metal on macOS, OpenGL elsewhere
./build/ViewerSCESS/ViewerSCESS --opengl     # or --metal / --vulkan
./build/SCESS/SCESS                          # prints the available methods, molecules and basis sets
./build/SCESS/SCESS RHF H2O sto-3g
./build/SCESS/SCESS CCSD_F12b H2O cc-pVDZ-F12                      # CABS defaults to cc-pVDZ-F12-OptRI
./build/SCESS/SCESS CCSD_F12b H2O cc-pVDZ --CABS aug-cc-pVDZ-OptRI
```

## Basis sets
The basis sets are Basis Set Exchange JSON files in `Assets/Basis`, one folder per family. The folder is the family
shown in the viewer and the help text, and the file name (without `.json`) is the basis set name. Names are
case-insensitive, and `6-31G*`-style names are accepted for the `6-31G(d)`-style files.

| Folder      | Basis sets                                                                 |
|-------------|----------------------------------------------------------------------------|
| `Pople`     | STO-3G, 3-21G, 6-21G, 6-31G, 6-31G(d), 6-31G(d,p), 6-31+G(d), 6-311G(d,p), 6-311+G(d,p) |
| `Dunning`   | cc-pVXZ (X = D, T, Q, 5), aug-cc-pVXZ (X = D, T, Q)                        |
| `Karlsruhe` | def2-SVP, def2-SVPD, def2-TZVP, def2-TZVPP, def2-TZVPD, def2-QZVP, def2-QZVPD |
| `F12`       | cc-pVXZ-F12 (X = D, T, Q)                                                  |
| `CABS`      | cc-pVXZ-F12-OptRI and aug-cc-pVXZ-OptRI (X = D, T, Q)                      |

The new sets cover H–Kr where the family defines them (no ECPs are used). F12 methods need a complementary auxiliary
basis set: pass `--CABS` to the command line driver, or choose one in the viewer when an F12 method is selected. By
default the matching OptRI set is used (`cc-pVXZ-F12` → `cc-pVXZ-F12-OptRI`, `cc-pVXZ` and `aug-cc-pVXZ` →
`aug-cc-pVXZ-OptRI`). Any basis set can be used as the CABS; the OBS is projected out of the union of the two.
All integrals are evaluated with Libint2 on every hardware thread, so basis sets are limited to its maximum angular
momentum (4 for a default build), which excludes cc-pV5Z and the QZ OptRI sets. Shells use solid harmonic (spherical)
functions unless the basis set file marks them as Cartesian (`gto_cartesian`, as for the Pople sets' d shells), as in
PySCF; Scripts/verification/spherical_check.py gives the PySCF references used by the regression tests.

To add a basis set, download it from the Basis Set Exchange in JSON format and put it in one of these folders (or a
new one).

Both programs find the `Assets` folder (basis sets, molecules and periodic table data) in this order:

1. The `assetsPath` in a `settings.json` in the working directory or next to the executable
2. `Contents/Resources/Assets` inside the app bundle
3. An `Assets` folder next to the executable or in the working directory
4. The source tree's `Assets` folder

The viewer saves its window layout to `~/Library/Application Support/SCESS/imgui.ini` (macOS),
`~/.config/SCESS/imgui.ini` (Linux) or `%APPDATA%\SCESS\imgui.ini` (Windows). On the first run it
starts from the default layout in `Assets/imgui.ini`.

# Project layout
| Path                 | Contents                                                              |
|----------------------|-----------------------------------------------------------------------|
| `src/Quantum`        | The quantum chemistry methods, integrals and basis sets               |
| `src/SCESS`          | The command line driver                                               |
| `src/Viewer`         | The graphical viewer                                                  |
| `src/Resources`      | Asset lookup and the viewer's shaders (GLSL compiled to SPIR-V)       |
| `Assets`             | Basis sets, molecules, data, icons and the default viewer layout      |
| `cmake`              | The per-target CMake scripts and helpers                              |
| `Libraries`          | Git submodules: Eigen, cppitertools, nlohmann JSON, ImPlot and nfd    |
| `Scripts`            | Build, package and icon scripts, plus the verification scripts        |

The icon is drawn in `Assets/Icon/SCESS Icon.svg`. After changing it, run
`Scripts/Generate_Icons.sh` (needs [Inkscape](https://inkscape.org)) to regenerate the splash
screen image (`SCESS Icon.png`) and the macOS app icon (`SCESS Icon.icns`).

# Dependencies
Built as part of SCESS:
- [Dahlia](https://github.com/KingfisherA/GraphicsEngine) (with Dear ImGui and SPIRV-Cross)
- [Eigen](https://gitlab.com/libeigen/eigen)
- [cppitertools](https://github.com/ryanhaining/cppitertools)
- [nlohmann JSON](https://github.com/nlohmann/json)
- [ImPlot](https://github.com/epezent/implot)
- [nativefiledialog](https://github.com/mlabbe/nativefiledialog)
- ImSpinner

Installed separately: Boost, fmt, LAPACK/LAPACKE, Libxc, Libint2, GLFW, GLEW, GLM and FreeType.

# Warnings
- This software is missing references to some of the modified code it uses. This will be fixed soon.
- The basis set assets do not have individual references, but all of them are derived from the
  [Basis Set Exchange](https://www.basissetexchange.org/).
- This software is unstable and may crash. Every method stores the AO two-electron integrals (N⁴ doubles, 1.4 GB for
  N = 115), and full CI scales factorially with the number of orbitals. Use at your own risk! (H₂O/cc-pVTZ CCSD(T)
  needs about 650 MB, H₂O/6-31G full CI about 730 MB, and H₂O/cc-pVTZ-F12 CCSD-F12b about 2.5 GB.)
