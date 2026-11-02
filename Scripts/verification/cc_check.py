"""Reproduce SCESS's spin-orbital CCSD / (T) algorithms in numpy and compare to PySCF."""
import numpy as np
from pyscf import gto, scf, cc

A2B = 1.8897259886
mol = gto.M(atom=[["H", (0.7569520, -0.586349, 0.0)], ["O", (0, 0, 0)], ["H", (-0.7569520, -0.586349, 0.0)]],
            basis="sto-3g", unit="Angstrom", verbose=0)
mf = scf.RHF(mol).run(conv_tol=1e-12)
mycc = cc.CCSD(mf).run(conv_tol=1e-12)
et_ref = mycc.ccsd_t()

C, eps = mf.mo_coeff, mf.mo_energy
nmo, nocc_sp = C.shape[1], 2 * mol.nelectron // 2
eri = mol.ao2mo(C, compact=False).reshape([nmo] * 4)  # chemists (pq|rs)

# AntisymmetriseIntegrateOutSpinClosedShell (alternating alpha/beta)
n = 2 * nmo
p = np.arange(n)
sp = p % 2
sa = p // 2
coul = (sp[:, None, None, None] == sp[None, None, :, None]) & (sp[None, :, None, None] == sp[None, None, None, :])
exch = (sp[:, None, None, None] == sp[None, None, None, :]) & (sp[None, :, None, None] == sp[None, None, :, None])
g = (eri[np.ix_(sa, sa, sa, sa)].transpose(0, 2, 1, 3) * coul
     - eri[np.ix_(sa, sa, sa, sa)].transpose(0, 2, 3, 1) * exch)
fd = np.repeat(eps, 2)
o, v = slice(0, nocc_sp), slice(nocc_sp, n)
no, nv = nocc_sp, n - nocc_sp
Dia = fd[o, None] - fd[None, v]
Dijab = fd[o, None, None, None] + fd[None, o, None, None] - fd[None, None, v, None] - fd[None, None, None, v]


def ccsd(wmbej_factor):
    t1 = np.zeros((no, nv))
    t2 = g[o, o, v, v] / Dijab
    E_old = 0
    for _ in range(200):
        tt = t2 + 0.5 * (np.einsum('ia,jb->ijab', t1, t1) - np.einsum('ib,ja->ijab', t1, t1))
        tau = t2 + np.einsum('ia,jb->ijab', t1, t1) - np.einsum('ib,ja->ijab', t1, t1)
        Fae = np.einsum('mf,mafe->ae', t1, g[o, v, v, v]) - 0.5 * np.einsum('mnaf,mnef->ae', tt, g[o, o, v, v])
        Fmi = np.einsum('ne,mnie->mi', t1, g[o, o, o, v]) + 0.5 * np.einsum('inef,mnef->mi', tt, g[o, o, v, v])
        Fme = np.einsum('nf,mnef->me', t1, g[o, o, v, v])
        Wmnij = g[o, o, o, o].copy()
        X = np.einsum('je,mnie->mnij', t1, g[o, o, o, v]); Wmnij += X - X.transpose(0, 1, 3, 2)
        Wmnij += 0.25 * np.einsum('ijef,mnef->mnij', tau, g[o, o, v, v])
        Wabef = g[v, v, v, v].copy()
        X = np.einsum('mb,amef->abef', t1, g[v, o, v, v]); Wabef += -X + X.transpose(1, 0, 2, 3)
        Wabef += 0.25 * np.einsum('mnab,mnef->abef', tau, g[o, o, v, v])
        Wmbej = g[o, v, v, o].copy()
        Wmbej += np.einsum('jf,mbef->mbej', t1, g[o, v, v, v])
        Wmbej -= np.einsum('nb,mnej->mbej', t1, g[o, o, v, o])
        tmp = 0.5 * t2 + np.einsum('jf,nb->jnfb', t1, t1)
        Wmbej -= wmbej_factor * np.einsum('jnfb,mnef->mbej', tmp, g[o, o, v, v])
        r1 = (np.einsum('ie,ae->ia', t1, Fae) - np.einsum('ma,mi->ia', t1, Fmi)
              + np.einsum('imae,me->ia', t2, Fme) - np.einsum('nf,naif->ia', t1, g[o, v, o, v])
              - 0.5 * np.einsum('imef,maef->ia', t2, g[o, v, v, v])
              - 0.5 * np.einsum('mnae,nmei->ia', t2, g[o, o, v, o]))
        r2 = g[o, o, v, v].copy()
        X = np.einsum('ijae,be->ijab', t2, Fae - 0.5 * np.einsum('mb,me->be', t1, Fme)); r2 += X - X.transpose(0, 1, 3, 2)
        X = np.einsum('imab,mj->ijab', t2, Fmi + 0.5 * np.einsum('me,je->mj', Fme, t1)); r2 += -X + X.transpose(1, 0, 2, 3)
        r2 += 0.5 * np.einsum('mnab,mnij->ijab', tau, Wmnij) + 0.5 * np.einsum('ijef,abef->ijab', tau, Wabef)
        X = np.einsum('imae,mbej->ijab', t2, Wmbej) - np.einsum('ie,ma,mbej->ijab', t1, t1, g[o, v, v, o])
        r2 += X - X.transpose(0, 1, 3, 2) - X.transpose(1, 0, 2, 3) + X.transpose(1, 0, 3, 2)
        X = np.einsum('ie,abej->ijab', t1, g[v, v, v, o]); r2 += X - X.transpose(1, 0, 2, 3)
        X = np.einsum('ma,mbij->ijab', t1, g[o, v, o, o]); r2 += -X + X.transpose(0, 1, 3, 2)
        t1, t2 = r1 / Dia, r2 / Dijab
        E = 0.25 * np.einsum('ijab,ijab', g[o, o, v, v], t2) + 0.5 * np.einsum('ijab,ia,jb', g[o, o, v, v], t1, t1)
        if abs(E - E_old) < 1e-12:
            break
        E_old = E
    return E, t1, t2


def P_ijk(T, which):  # which: list of transpositions applied to occ axes (0,1,2)
    out = T.copy()
    for ax in which:
        perm = list(range(6)); perm[ax[0]], perm[ax[1]] = perm[ax[1]], perm[ax[0]]
        out -= T.transpose(perm)
    return out


def apply(T, occ_swaps, vir_swaps):
    return P_ijk(P_ijk(T, occ_swaps), vir_swaps)


D3 = (fd[o, None, None, None, None, None] + fd[None, o, None, None, None, None] + fd[None, None, o, None, None, None]
      - fd[None, None, None, v, None, None] - fd[None, None, None, None, v, None] - fd[None, None, None, None, None, v])
I_JK = [(0, 1), (0, 2)]   # P(i/jk) as written in SCESS: 1 - P_ij - P_ik
A_BC = [(3, 4), (3, 5)]   # P(a/bc): 1 - P_ab - P_ac
K_IJ = [(0, 2), (1, 2)]   # P(k/ij)
C_AB = [(3, 5), (4, 5)]   # P(c/ab)


def triples_scess(t1, t2, fixed):
    """SCESS CalculateTriplesCorrection; fixed=True uses the correct permutation operators."""
    X1 = np.einsum('dkbc,ijad->ijkabc', g[v, o, v, v], t2)
    X2 = np.einsum('jklc,ilab->ijkabc', g[o, o, o, v], t2)
    if fixed:
        t3 = (apply(X1, K_IJ, A_BC) - apply(X2, I_JK, C_AB)) / D3
    else:
        t3 = (apply(X1, I_JK, A_BC) - apply(X2, I_JK, A_BC)) / D3
    L = np.einsum('efbm,ijmaef->ijab', g[v, v, v, o], t3); L = L - L.transpose(0, 1, 3, 2)
    R = np.einsum('jemn,imnabe->ijab', g[o, v, o, o], t3); R = R - R.transpose(1, 0, 2, 3)
    VT3 = 0.5 * L - 0.5 * R
    E4_code = 0.25 * np.einsum('ijab,ijab', g[o, o, v, v], VT3 / Dijab)  # bra = first-order t2
    E4_std = 0.25 * np.einsum('ijab,ijab', t2, VT3)                       # bra = converged CCSD t2
    E5 = 0.25 * np.einsum('bcjk,ijkabc,ia', g[v, v, o, o], t3, t1)
    return E4_code + E5, E4_std + E5


print(f"PySCF   E(CCSD corr) = {mycc.e_corr:.10f}   E(T) = {et_ref:.10f}")
for fac in (0.25, 1.0):
    E, t1, t2 = ccsd(fac)
    print(f"Wmbej factor {fac:4}: E(CCSD corr) = {E:.10f}  (diff {E - mycc.e_corr:+.2e})")
E, t1, t2 = ccsd(1.0)
for fixed in (False, True):
    code, std = triples_scess(t1, t2, fixed)
    print(f"(T) perms {'fixed' if fixed else 'as-is'}: E(T) code-bra = {code:.10f} (diff {code - et_ref:+.2e})"
          f"   CCSD-t2 bra = {std:.10f} (diff {std - et_ref:+.2e})")
