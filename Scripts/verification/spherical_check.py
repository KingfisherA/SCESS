"""PySCF references for SCESS with solid harmonic (spherical) basis functions.

The basis sets are read from SCESS's own Basis Set Exchange JSON files so that both codes use identical contractions.
Prints the H2O / cc-pVDZ RHF, UHF, MP2, CCSD, CCSD(T) and LDA (Slater + VWN5) energies and the H2O / 6-31G CISD and FCI
energies at the built-in water geometry.
"""
import json
import os
import sys

from pyscf import gto, scf, mp, cc, ci, dft, fci

ASSETS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "Assets", "Basis")
WATER = [["H", (0.7569520, -0.586349, 0.0)], ["O", (0, 0, 0)], ["H", (-0.7569520, -0.586349, 0.0)]]
SYMBOLS = {1: "H", 8: "O"}


def load_basis(family, name):
    with open(os.path.join(ASSETS, family, name + ".json")) as f:
        data = json.load(f)
    basis = {}
    for z, element in data["elements"].items():
        if int(z) not in SYMBOLS:
            continue
        shells = []
        for shell in element["electron_shells"]:
            exps = [float(e) for e in shell["exponents"]]
            ams = shell["angular_momentum"]
            for j, coeffs in enumerate(shell["coefficients"]):
                l = ams[j] if len(ams) > 1 else ams[0]
                prims = [[e, float(c)] for e, c in zip(exps, coeffs) if float(c) != 0.0]
                shells.append([l] + prims)
        basis[SYMBOLS[int(z)]] = shells
    return basis


def water(family, name):
    return gto.M(atom=WATER, basis=load_basis(family, name), unit="Angstrom", verbose=0, cart=False)


def main():
    mol = water("Dunning", "cc-pVDZ")
    rhf = scf.RHF(mol).run(conv_tol=1e-12)
    uhf = scf.UHF(mol).run(conv_tol=1e-12)
    mp2 = mp.MP2(rhf).run()
    ccsd = cc.CCSD(rhf).run(conv_tol=1e-12, conv_tol_normt=1e-10)
    et = ccsd.ccsd_t()
    lda = dft.RKS(mol)
    lda.xc = "slater,vwn5"
    lda.grids.level = 9
    lda.run(conv_tol=1e-12)

    print("H2O / cc-pVDZ (spherical)")
    print(f"  E(RHF)  = {rhf.e_tot:.10f}")
    print(f"  E(UHF)  = {uhf.e_tot:.10f}")
    print(f"  E(MP2)  = {mp2.e_corr:.10f}")
    print(f"  E(CCSD) = {ccsd.e_corr:.10f}")
    print(f"  E(T)    = {et:.10f}")
    print(f"  E(LDA)  = {lda.e_tot:.10f}")
    sys.stdout.flush()

    mol = water("Pople", "6-31G")
    rhf = scf.RHF(mol).run(conv_tol=1e-12)
    e_cisd = ci.CISD(rhf).run(conv_tol=1e-12).e_tot
    e_fci = fci.FCI(rhf).kernel()[0]
    print("H2O / 6-31G")
    print(f"  E(RHF)  = {rhf.e_tot:.10f}")
    print(f"  E(CISD) = {e_cisd:.10f}")
    print(f"  E(FCI)  = {e_fci:.10f}")


if __name__ == "__main__":
    main()
