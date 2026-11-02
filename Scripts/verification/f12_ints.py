"""Gaussian-geminal F12 integrals (libint cgtg / cgtg_x_coulomb / delcgtg2 equivalents) built from
PySCF range-separated ERIs:  d/dmu [erf(mu r)/r] = (2/sqrt(pi)) exp(-mu^2 r^2)."""
import numpy as np
from functools import lru_cache
from pyscf import gto

SQPI = np.sqrt(np.pi)
# SCESS CalculateSlaterContraction (gamma = 1): exponent, coefficient * (-1/gamma)
STG6G = [(0.2209, 0.3144), (1.004, 0.3037), (3.622, 0.1681),
         (12.16, 0.09811), (45.87, 0.06024), (254.4, 0.03726)]
CGTG = [(a, -c) for a, c in STG6G]


def build_mol(geom, obs, cabs, unit="Angstrom"):
    mo = gto.M(atom=geom, basis=obs, unit=unit, verbose=0)
    mc = gto.M(atom=geom, basis=cabs, unit=unit, verbose=0)
    both = gto.conc_mol(mo, mc)
    both._atm[mo.natm:, gto.CHARGE_OF] = 0  # do not double count nuclei
    return mo, both


class GeminalInts:
    def __init__(self, mol):
        self.mol = mol

        @lru_cache(maxsize=None)
        def E(mu):  # (pq|erf(mu r)/r|rs), chemists'
            with mol.with_range_coulomb(mu):
                return mol.intor("int2e")
        self.E = E

    def dE(self, mu):
        h = 1e-3 * mu
        E = self.E
        return (E(mu - 2 * h) - 8 * E(mu - h) + 8 * E(mu + h) - E(mu + 2 * h)) / (12 * h)

    def d2E(self, mu):
        h = 2e-3 * mu
        E = self.E
        return (-E(mu - 2 * h) + 16 * E(mu - h) - 30 * E(mu) + 16 * E(mu + h) - E(mu + 2 * h)) / (12 * h * h)

    def gauss(self, beta):            # exp(-beta r^2)
        return 0.5 * SQPI * self.dE(np.sqrt(beta))

    def r2gauss(self, beta):          # r^2 exp(-beta r^2)
        mu = np.sqrt(beta)
        return -SQPI / (4 * mu) * self.d2E(mu)

    def gauss_over_r(self, alpha, n=96):  # exp(-alpha r^2)/r = (2/sqrt(pi)) int_0^inf exp(-(alpha+t^2) r^2) dt
        x, w = np.polynomial.legendre.leggauss(n)
        s = 0.5 * (x + 1); w = 0.5 * w
        c = 4.0
        t = c * s / (1 - s); dt = c / (1 - s) ** 2
        out = 0.0
        for ti, wi, di in zip(t, w, dt):
            out = out + (2 / SQPI) * wi * di * self.gauss(alpha + ti * ti)
        return out

    def all(self, cgtg=CGTG):
        f12 = sum(c * self.gauss(a) for a, c in cgtg)
        f12sq = sum(ci * cj * self.gauss(ai + aj) for ai, ci in cgtg for aj, cj in cgtg)
        f12g = sum(c * self.gauss_over_r(a) for a, c in cgtg)
        delsq = sum(4 * ai * aj * ci * cj * self.r2gauss(ai + aj) for ai, ci in cgtg for aj, cj in cgtg)
        return f12, f12sq, f12g, delsq


if __name__ == "__main__":
    # self tests: Gaussian-over-r with alpha->0 must reproduce plain 1/r ERIs
    _, both = build_mol("H 0 0 0; H 0.7122300374 0 0", "sto-3g", "cc-pvdz")
    gi = GeminalInts(both)
    eri = both.intor("int2e")
    print("1/r via quadrature   max err:", np.abs(gi.gauss_over_r(1e-12) - eri).max())
    # exp(-b r^2) * exp(-c r^2) consistency: r^2 e^{-b r^2} = -(d/db) e^{-b r^2}
    b, h = 1.3, 1e-4
    fd = -(gi.gauss(b + h) - gi.gauss(b - h)) / (2 * h)
    print("r^2 gauss vs FD       max err:", np.abs(gi.r2gauss(b) - fd).max())
