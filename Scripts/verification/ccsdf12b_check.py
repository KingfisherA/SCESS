"""Independent check of SCESS CCSD-F12b (src/Quantum/F12/CCSD_F12b.cpp).

The F12 intermediates come from the NumPy port in mp2f12_port.py (PySCF integrals) and the CCSD equations are PySCF's
own spin-adapted RCCSD. The geminal coupling is added to PySCF's amplitude update and energy:
  doubles residual: += sum_kl c_ij^kl (C_kl^ab + V_ab^kl),  c_ij^kl = 3/8 d_ik d_jl + 1/8 d_il d_jk
  energy:           += sum_ijab (2 Gc_ij^ab - Gc_ij^ba) t_ij^ab,  Gc_ij^ab = sum_kl c_ij^kl C_kl^ab
with V_ab^kl = <ab|r12^-1 f12|kl> - sum_rs <ab|r12^-1|rs><rs|f12|kl> (no CABS sums, as in CCSD-F12x).

Usage: python3 ccsdf12b_check.py H2 [OBS CABS]
"""
import sys
import numpy as np
from pyscf import cc

from mp2f12_port import setup, transform_f12_basis, cabs_singles, intermediates, e_fixed, ein


def fixed_amplitudes(X_klab):
    """sum_kl c_ij^kl X_kl^ab = 3/8 X_ij^ab + 1/8 X_ji^ab"""
    return 3 / 8 * X_klab + 1 / 8 * X_klab.transpose(1, 0, 2, 3)


def ccsd_f12b(d, m):
    nocc, nO = m["nocc"], m["nO"]
    o, v, P = slice(0, nocc), slice(nocc, nO), slice(0, nO)
    I = intermediates(m, trace_bug=False)
    V_abkl = m["fg"][v, v, o, o] - ein("abrs,rskl->abkl", m["eri"][v, v, P, P], m["f"][P, P, o, o])
    G = fixed_amplitudes(I["C"] + V_abkl.transpose(2, 3, 0, 1))
    Gc = fixed_amplitudes(I["C"])

    mycc = cc.RCCSD(d["mf"])
    mycc.conv_tol, mycc.conv_tol_normt = 1e-11, 1e-9
    eps = d["mf"].mo_energy
    D = eps[o, None, None, None] + eps[None, o, None, None] - eps[None, None, v, None] - eps[None, None, None, v]

    base_update, base_energy = mycc.update_amps, mycc.energy

    def update_amps(t1, t2, eris):
        t1new, t2new = base_update(t1, t2, eris)
        # PySCF divides its residual by D last, so this adds G to the driving term
        return t1new, t2new + G / D

    def coupling(t2):
        return ein("ijab,ijab->", 2 * Gc - Gc.transpose(0, 1, 3, 2), t2)

    def energy(t1=None, t2=None, eris=None):
        t1 = mycc.t1 if t1 is None else t1
        t2 = mycc.t2 if t2 is None else t2
        return base_energy(t1, t2, eris) + coupling(t2)

    mycc.update_amps, mycc.energy = update_amps, energy
    mycc.verbose = 0
    mycc.kernel()
    assert mycc.converged
    return mycc.e_corr, coupling(mycc.t2), e_fixed(I, nocc)


if __name__ == "__main__":
    systems = {"H2": "H 0 0 0; H 0.7122300374 0 0",
               "H2O": "H 0.7569520 -0.586349 0; O 0 0 0; H -0.7569520 -0.586349 0"}
    name = sys.argv[1] if len(sys.argv) > 1 else "H2"
    obs, cabs = (sys.argv[2], sys.argv[3]) if len(sys.argv) > 3 else ("sto-3g", "cc-pvdz")
    d = setup(systems[name], obs, cabs)
    m = transform_f12_basis(d)
    e_cabs = cabs_singles(m)
    e_ccsd, e_coupling, e_f12 = ccsd_f12b(d, m)
    print(f"===== {name} {obs} / CABS {cabs}: nOBS={m['nO']} nCABS={m['N'] - m['nO']} nOcc={m['nocc']}")
    print(f"E(HF)                          = {m['E_HF']:.10f}")
    print(f"E(CABS singles)                = {e_cabs:.10f}")
    print(f"E(CCSD, with geminal coupling) = {e_ccsd:.10f}")
    print(f"    of which E(coupling)       = {e_coupling:.10f}")
    print(f"E(F12 geminals)                = {e_f12:.10f}")
    print(f"E(CCSD-F12b correlation)       = {e_ccsd + e_f12:.10f}")
    print(f"E(Total)                       = {m['E_HF'] + e_cabs + e_ccsd + e_f12:.10f}")
