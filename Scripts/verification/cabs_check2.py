"""Replicate SCESS TransformF12Basis + EvaluateCABSSinglesCorrection for H2O sto-3g / cc-pVDZ CABS."""
import numpy as np
from pyscf import gto, scf
geom = "H 0 0 0; H 0.7122300374 0 0"
obs = gto.M(atom=geom, basis="sto-3g", verbose=0)
cabs = gto.M(atom=geom, basis="cc-pvdz", verbose=0)
mf = scf.RHF(obs).run(conv_tol=1e-12)
both = gto.conc_mol(obs, cabs)
both._atm[obs.natm:, gto.CHARGE_OF] = 0  # remove duplicated nuclei
nO, N, nocc = obs.nao, obs.nao + cabs.nao, obs.nelectron // 2
S, h = both.intor("int1e_ovlp"), both.intor("int1e_kin") + both.intor("int1e_nuc")
eri = both.intor("int2e")                             # chemists'
Da = mf.mo_coeff[:, :nocc] @ mf.mo_coeff[:, :nocc].T  # D_Alpha (no factor 2)
J = np.einsum("ijmn,mn->ij", eri[:, :, :nO, :nO], Da)
K = np.einsum("inmj,mn->ij", eri[:, :nO, :nO, :], Da)
F = h + 2 * J - K

def scess_Q(sort_fix):
    P = np.zeros((N, N)); P[nO:, nO:] = np.eye(N - nO)
    P[:nO, :] -= np.linalg.solve(S[:nO, :nO], S[:nO, :])
    P[np.arange(nO), np.arange(nO)] = 0.0              # "if (i != mu)" skip -> diagonal stays 0
    s, U = np.linalg.eigh(P.T @ S @ P)
    X = U * np.where(np.abs(s) > 1e-8, 1 / np.sqrt(np.abs(s).clip(1e-300)), 0.0)
    Pp = P @ X
    e, V = np.linalg.eigh(Pp.T @ F @ Pp)
    Q = Pp @ V
    null = np.linalg.norm(Q, axis=0) < 1e-10
    print(f"  #null columns = {null.sum()} (nOBS = {nO}); Fock eigenvalues < 0 among non-null: {(e[~null] < 0).sum()}")
    print(f"  null-column positions: {np.where(null)[0].tolist()}")
    if sort_fix:
        Q = np.hstack([np.zeros((N, nO)), Q[:, ~null]])
    Q[:nO, :nO] = mf.mo_coeff
    return Q

def cabs_singles(Q):
    Fm = Q.T @ F @ Q
    nv = nO - nocc; nx = N - nO
    Dai = Fm.diagonal()[:nocc][None, :] - Fm.diagonal()[nocc:nO][:, None]
    Dxi = Fm.diagonal()[:nocc][None, :] - Fm.diagonal()[nO:][:, None]
    txi = Fm[nO:, :nocc].copy(); tai = np.zeros((nv, nocc))
    for _ in range(100):
        tai = Fm[nocc:nO, nO:] @ txi / Dai
        txi = (Fm[nO:, :nocc] + Fm[nO:, nocc:nO] @ tai) / Dxi
    return 2 * np.einsum("ix,xi", Fm[:nocc, nO:], txi)

print("As implemented:"); Q = scess_Q(False); print(f"  E(CABS singles) = {cabs_singles(Q):.6f} Eh")
print("Null space removed explicitly:"); Q = scess_Q(True); print(f"  E(CABS singles) = {cabs_singles(Q):.6f} Eh")

for mol_name, g in [("H2", "H 0 0 0; H 0.7122300374 0 0")]:
    pass
