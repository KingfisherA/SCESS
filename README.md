
# SCESS
*Schrödinger's Cat's Electronic System Solver*

<div style="text-align:center">
<img src="./Assets/Icon/SCESS Small.png" width="150" height="150" /> 
</div>

# About
This is a basic implementation of different *ab-initio* quantum chemistry theories with a graphical editor to explore the solutions.

# Features
- Allows loading and displaying XYZ molecule files
- Can perform RHF & UHF with DIIS acceleration
- Can perform Moller Plesset (up to order 2), CCSD and CCSD(T)
- Can draw isosurface plots of the Hartree-Fock solutions

# Binaries
Release binaries are available for testing from the releases GitHub section.

# Dependencies
- Eigen
- GraphicsEngine

# Warnings
- Currently this software is missing a lot of references to modified code it is using (especially regarding the integrator). This will be fixed soon.
- The basis sets and font assets do not have appropriate references either, however the basis sets are all derived from [basis set exchange](https://www.basissetexchange.org/). This will be updated soon.
- This software is extremely unstable and may cause crashes. The coupled cluster theories are particularly bad at hogging large amounts of memory. Use at your own risk!
- DFT is completely broken at the moment
