"""Faithful numpy port of SCESS TransformF12Basis, EvaluateCABSSinglesCorrection, MP2_F12_Fixed,
MP2_F12_Optimise, plus independent cross-checks."""
import sys, numpy as np
from pyscf import scf
from f12_ints import build_mol, GeminalInts

ein = np.einsum


def setup(geom, obs="sto-3g", cabs="cc-pvdz"):
    mo, both = build_mol(geom, obs, cabs)
    mf = scf.RHF(mo).run(conv_tol=1e-12)
    nO, N, nocc = mo.nao, both.nao, mo.nelectron // 2
    S, T, V = both.intor("int1e_ovlp"), both.intor("int1e_kin"), both.intor("int1e_nuc")
    eri = both.intor("int2e")
    f12, f12sq, f12g, delsq = GeminalInts(both).all()
    return dict(mf=mf, nO=nO, N=N, nocc=nocc, S=S, T=T, V=V, eri=eri,
                f12=f12, f12sq=f12sq, f12g=f12g, delsq=delsq)


def transform_f12_basis(d):
    """Port of Transforms::TransformF12Basis (convertChemToPhys = true)."""
    nO, N, S, eri = d["nO"], d["N"], d["S"], d["eri"]
    C, nocc = d["mf"].mo_coeff, d["nocc"]
    SInv = np.linalg.inv(S[:nO, :nO])
    P = np.zeros((N, N))
    for i in range(N):
        if i >= nO: P[i, i] = 1.0
        for mu in range(nO):
            c = -SInv[mu, :] @ S[:nO, i]
            if i != mu: P[mu, i] = c
    s, U = np.linalg.eigh(P.T @ S @ P)
    X = U @ np.diag([1 / np.sqrt(v) if abs(v) > 1e-8 else 0.0 for v in s])
    Pp = P @ X
    D = C[:, :nocc] @ C[:, :nocc].T
    J = ein("ijmn,mn->ij", eri[:, :, :nO, :nO], D)
    K = ein("inmj,mn->ij", eri[:, :nO, :nO, :], D)
    F = d["T"] + d["V"] + 2 * J - K
    e, W = np.linalg.eigh(Pp.T @ F @ Pp)
    Q = Pp @ W
    Q[:nO, :nO] = C
    tr2 = lambda m: Q.T @ m @ Q
    tr4 = lambda t: ein("pqrs,pi,qj,rk,sl->ijkl", t, Q, Q, Q, Q, optimize=True).transpose(0, 2, 1, 3)
    out = dict(F=tr2(F), h=tr2(d["T"] + d["V"]), eri=tr4(eri), f=tr4(d["f12"]), f2=tr4(d["f12sq"]),
               fg=tr4(d["f12g"]), u=tr4(d["delsq"]))
    out.update(nO=nO, N=N, nocc=nocc, eps=d["mf"].mo_energy, E_HF=d["mf"].e_tot)
    return out


def cabs_singles(m):
    F, nocc, nO, N = m["F"], m["nocc"], m["nO"], m["N"]
    o, v, x = slice(0, nocc), slice(nocc, nO), slice(nO, N)
    Dai = F.diagonal()[o][None, :] - F.diagonal()[v][:, None]
    Dxi = F.diagonal()[o][None, :] - F.diagonal()[x][:, None]
    txi = F[x, o].copy(); tai = np.zeros_like(Dai); last = np.zeros_like(txi)
    for it in range(100):
        tai = F[v, x] @ txi / Dai
        txi = (F[x, o] + F[x, v] @ tai) / Dxi
        if abs((txi - last).sum()) < 1e-8: break
        last = txi.copy()
    return 2 * ein("ix,xi", F[o, x], txi)


def intermediates(m, trace_bug=True, include_xy_K=True):
    """Port of the V, X, Fock, B, C construction shared by MP2_F12_Fixed and MP2_F12_Optimise.
    trace_bug=True reproduces SCESS's Einsum reduction semantics (sum over BOTH repeated indices)."""
    nocc, nO, N = m["nocc"], m["nO"], m["N"]
    o, v, P, x = slice(0, nocc), slice(nocc, nO), slice(0, nO), slice(nO, N)
    g, f, f2, fg, u = m["eri"], m["f"], m["f2"], m["fg"], m["u"]
    V = fg[o, o, o, o] - ein("ijrs,klrs->ijkl", g[o, o, P, P], f[o, o, P, P]) \
        - ein("ijxm,klxm->ijkl", g[o, o, x, o], f[o, o, x, o]) - ein("ijmx,klmx->ijkl", g[o, o, o, x], f[o, o, o, x])
    X = f2[o, o, o, o] - ein("ijrs,klrs->ijkl", f[o, o, P, P], f[o, o, P, P]) \
        - ein("ijxm,klxm->ijkl", f[o, o, x, o], f[o, o, x, o]) - ein("ijmx,klmx->ijkl", f[o, o, o, x], f[o, o, o, x])
    Fm = m["h"].copy()
    if trace_bug:
        Jfun = lambda blk: blk.sum(axis=(1, 3))      # Einsum<P,I,Q,I> -> sum_{i,i'}
        Kfun = lambda blk: blk.sum(axis=(0, 1))      # Einsum<I,I,P,Q> -> sum_{i,i'}
    else:
        Jfun = lambda blk: ein("piqi->pq", blk)
        Kfun = lambda blk: ein("iipq->pq", blk)
    Fm[P, P] += 2 * Jfun(g[P, o, P, o]); Fm[P, x] += 2 * Jfun(g[o, P, o, x].transpose(1, 0, 3, 2))
    Fm[x, x] += 2 * Jfun(g[o, x, o, x].transpose(1, 0, 3, 2))
    Km = np.zeros_like(Fm)
    Km[P, P] = Kfun(g[o, o, P, P]); Km[P, x] = Kfun(g[o, o, P, x]); Km[x, x] = Kfun(g[o, o, x, x])
    Fm -= Km
    Fm[x, P] = Fm[P, x].T; Km[x, P] = Km[P, x].T
    FK = Fm + Km
    B = u[o, o, o, o].copy()
    S2 = lambda tmp: ein("mnpq,klpq->klmn", f[o, o, P, P], tmp)
    B -= S2(ein("pr,lkqr->klpq", FK[P, P], f[o, o, P, P]))
    B -= S2(ein("qr,klpr->klpq", FK[P, P], f[o, o, P, P]))
    B -= S2(ein("px,lkqx->klpq", FK[P, x], f[o, o, P, x]))
    B -= S2(ein("qx,klpx->klpq", FK[P, x], f[o, o, P, x]))
    A_iy = lambda tmp: ein("mniy,kliy->klmn", f[o, o, o, x], tmp)
    A_xj = lambda tmp: ein("nmjx,klxj->klmn", f[o, o, o, x], tmp)
    B -= A_iy(ein("yp,klip->kliy", FK[x, P], f[o, o, o, P]))
    B -= A_xj(ein("xp,lkjp->klxj", FK[x, P], f[o, o, o, P]))
    B -= A_iy(ein("ip,klpy->kliy", FK[o, P], f[o, o, P, x]))
    B -= A_iy(ein("yx,klix->kliy", FK[x, x], f[o, o, o, x]))
    B -= A_xj(ein("xy,lkjy->klxj", FK[x, x], f[o, o, o, x]))
    B -= A_xj(ein("jp,lkpx->klxj", FK[o, P], f[o, o, P, x]))
    B -= A_iy(ein("ix,klxy->kliy", FK[o, x], f[o, o, x, x]))
    B -= A_xj(ein("jy,lkyx->klxj", FK[o, x], f[o, o, x, x]))
    Y_xb = lambda tmp: ein("klxb,nmbx->klmn", tmp, f[o, o, v, x])
    Y_ay = lambda tmp: ein("klay,mnay->klmn", tmp, f[o, o, v, x])
    B -= Y_xb(ein("xp,lkbp->klxb", Km[x, P], f[o, o, v, P]))
    B -= Y_ay(ein("yp,klap->klay", Km[x, P], f[o, o, v, P]))
    B -= Y_xb(ein("xy,lkby->klxb", Km[x, x], f[o, o, v, x]))
    B -= Y_xb(ein("bp,lkpx->klxb", Km[v, P], f[o, o, P, x]))
    B -= Y_ay(ein("ap,klpy->klay", Km[v, P], f[o, o, P, x]))
    B -= Y_ay(ein("yx,klax->klay", Km[x, x], f[o, o, v, x]))
    B -= Y_xb(ein("by,lkyx->klxb", Km[v, x], f[o, o, x, x]))
    B -= Y_ay(ein("ax,klxy->klay", Km[v, x], f[o, o, x, x]))
    if include_xy_K:
        XY = lambda tmp: ein("klxy,mnxy->klmn", tmp, f[o, o, x, x])
        B -= XY(ein("xp,klpy->klxy", Km[x, P], f[o, o, P, x]))
        B -= XY(ein("yp,lkpx->klxy", Km[x, P], f[o, o, P, x]))
        B -= XY(ein("yz,lkzx->klxy", Km[x, x], f[o, o, x, x]))
        B -= XY(ein("xz,klzy->klxy", Km[x, x], f[o, o, x, x]))
    B += ein("kp,nmlp->klmn", FK[o, P], f2[o, o, o, P]) + ein("lp,mnkp->klmn", FK[o, P], f2[o, o, o, P])
    B += ein("kx,nmlx->klmn", FK[o, x], f2[o, o, o, x]) + ein("lx,mnkx->klmn", FK[o, x], f2[o, o, o, x])
    C = ein("bx,klax->klab", Fm[v, x], f[o, o, v, x]) + ein("ax,lkbx->klab", Fm[v, x], f[o, o, v, x])
    B -= ein("klab,mnab->klmn", f[o, o, v, v], C)
    B = 0.5 * (B + B.transpose(2, 3, 0, 1))
    return dict(V=V, X=X, B=B, C=C, F=Fm, K=Km, FK=FK)


def e_mp2(m):
    nocc, nO, g, eps = m["nocc"], m["nO"], m["eri"], m["eps"]
    o, v = slice(0, nocc), slice(nocc, nO)
    D = eps[o, None, None, None] + eps[None, o, None, None] - eps[None, None, v, None] - eps[None, None, None, v]
    gv = g[o, o, v, v]
    return ein("ijab,ijab->", gv, 2 * gv - gv.transpose(0, 1, 3, 2) ) * 1.0 / 1.0 if False else \
        ein("ijab,ijab,ijab->", gv, 2 * gv - gv.transpose(0, 1, 3, 2), 1 / D)


def e_fixed(I, nocc):
    V, X, B, F = I["V"], I["X"], I["B"], I["F"]
    E = 0.0
    for i in range(nocc):
        for j in range(nocc):
            ev = 5 / 8 * V[i, j, i, j] - 1 / 8 * V[j, i, i, j]
            ex = (F[i, i] + F[j, j]) * (7 / 32 * X[i, j, i, j] + 1 / 32 * X[i, j, j, i])
            eb = 7 / 32 * B[i, j, i, j] + 1 / 32 * B[i, j, j, i]
            E += 2 * ev + eb - ex
    return E


def e_opt(m, I, v_gside_ij=False):
    """Spin-orbital-correct optimised-amplitude MP2-F12 (what MP2_F12_Optimise intends).
    v_gside_ij=False reproduces SCESS's V_(k,l,i,j) indexing; True uses V_(i,j,k,l)."""
    nocc, nO, g, F = m["nocc"], m["nO"], m["eri"], I["F"]
    V, X, B, C = I["V"], I["X"], I["B"], I["C"]
    o, v = slice(0, nocc), slice(nocc, nO)
    ev = F.diagonal()
    E = 0.0
    pairs = [(k, l) for k in range(nocc) for l in range(nocc)]
    tri = [(k, l) for k in range(nocc) for l in range(k + 1, nocc)]
    for i in range(nocc):
        for j in range(nocc):
            D = ev[v][:, None] + ev[v][None, :] - ev[i] - ev[j]
            gij = g[i, j, v, v]
            Vk = (lambda k, l: V[i, j, k, l]) if v_gside_ij else (lambda k, l: V[k, l, i, j])
            # opposite spin
            Vt = np.array([Vk(k, l) - (C[k, l] * gij / D).sum() for k, l in pairs])
            Bt = np.array([[B[k, l, mm, n] - (ev[i] + ev[j]) * X[k, l, mm, n] - (C[k, l] * C[mm, n] / D).sum()
                            for mm, n in pairs] for k, l in pairs])
            E += -Vt @ np.linalg.solve(Bt, Vt)
            if i < j:
                ut = np.triu(np.ones_like(D), 1)
                Ca = lambda k, l: C[k, l] - C[k, l].T
                ga = gij - gij.T
                Vt = np.array([Vk(k, l) - Vk(l, k) - (Ca(k, l) * ga * ut / D).sum() for k, l in tri])
                Bt = np.array([[B[k, l, mm, n] - B[k, l, n, mm] - (ev[i] + ev[j]) * (X[k, l, mm, n] - X[k, l, n, mm])
                                - (Ca(k, l) * Ca(mm, n) * ut / D).sum() for mm, n in tri] for k, l in tri])
                E += 2 * (-Vt @ np.linalg.solve(Bt, Vt))
    return E


def independent_B(m, I, include_xy_K=True):
    """B from the thesis/Valeev 3C formula written directly over full index ranges (not mirroring the code)."""
    nocc, nO, N = m["nocc"], m["nO"], m["N"]
    o, v, x = slice(0, nocc), slice(nocc, nO), slice(nO, N)
    f, f2, u = m["f"], m["f2"], m["u"]
    FK, K, F = I["FK"], I["K"], I["F"]
    fkl = f[o, o]                                                 # f_kl^{alpha beta}, all alpha, beta
    rbar = ein("ag,klgb->klab", FK, fkl) + ein("bg,klag->klab", FK, fkl)
    rtil = ein("ag,klgb->klab", K, fkl) + ein("bg,klag->klab", K, fkl)
    Pm = np.zeros((N, N)); Pm[:nO, :nO] = 1; Pm[o, x] = 1; Pm[x, o] = 1   # P12 = pq + i a' + a' i
    Ym = np.zeros((N, N)); Ym[v, x] = 1; Ym[x, v] = 1                    # Y: a a' + a' a
    if include_xy_K: Ym[x, x] = 1
    A = u[o, o, o, o] + ein("kp,plmn->klmn", FK[o], f2[:, o, o, o]) + ein("lp,kpmn->klmn", FK[o], f2[o, :, o, o])
    A -= ein("mnab,klab,ab->klmn", fkl, rtil, Ym) + ein("mnab,klab,ab->klmn", fkl, rbar, Pm)
    Cm = ein("ax,klxb->klab", F[v, x], f[o, o, x, v]) + ein("bx,klax->klab", F[v, x], f[o, o, v, x])
    Bm = A - ein("klab,mnab->klmn", f[o, o, v, v], Cm)
    return 0.5 * (Bm + Bm.transpose(2, 3, 0, 1)), Cm


if __name__ == "__main__":
    systems = {"H2": "H 0 0 0; H 0.7122300374 0 0",
               "H2O": "H 0.7569520 -0.586349 0; O 0 0 0; H -0.7569520 -0.586349 0"}
    for name in sys.argv[1:]:
        d = setup(systems[name]); m = transform_f12_basis(d)
        print(f"===== {name}: nOBS={m['nO']} nCABS={m['N']-m['nO']} nOcc={m['nocc']}")
        print(f"E(HF)            = {m['E_HF']:.12f}")
        print(f"E(CABS singles)  = {cabs_singles(m):.12f}")
        print(f"E(MP2)           = {e_mp2(m):.12f}")
        Ib = intermediates(m, trace_bug=True); Ic = intermediates(m, trace_bug=False)
        print(f"Fock trace bug: max|F_code - F_correct| = {np.abs(Ib['F'] - Ic['F']).max():.3e};"
              f"  max|F_correct - F(TransformF12Basis)| = {np.abs(Ic['F'] - m['F']).max():.3e}")
        print(f"E(F12 fixed)  as coded = {e_fixed(Ib, m['nocc']):.12f}   with correct Fock = {e_fixed(Ic, m['nocc']):.12f}")
        Bi, Ci = independent_B(m, Ic)
        print(f"independent B vs port B: {np.abs(Bi - Ic['B']).max():.3e}   C: {np.abs(Ci - Ic['C']).max():.3e}")
        for lbl, I in (("as coded", Ib), ("correct Fock", Ic)):
            print(f"E(F12 opt, {lbl:12s}) V(k,l,i,j) = {e_opt(m, I):.12f}   V(i,j,k,l) = {e_opt(m, I, True):.12f}")
        I0 = intermediates(m, trace_bug=False, include_xy_K=False)
        print(f"without a'b' K-term (thesis formula): fixed = {e_fixed(I0, m['nocc']):.12f}  opt = {e_opt(m, I0):.12f}")
        V = Ic["V"]
        print(f"symmetry: |V_ijkl - V_klij| max = {np.abs(V - V.transpose(2,3,0,1)).max():.3e};"
              f" |B_klmn - B_lknm| max = {np.abs(Ic['B'] - Ic['B'].transpose(1,0,3,2)).max():.3e};"
              f" |X sym| = {np.abs(Ic['X'] - Ic['X'].transpose(2,3,0,1)).max():.3e}")
