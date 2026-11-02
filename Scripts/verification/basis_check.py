import json, os, numpy as np
from pyscf import gto, scf
REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
bse = json.load(open(os.path.join(REPO, 'Assets', 'Basis', 'sto-3g.json')))
def scess_basis(Z):
    """Mimic BasisSet::LoadAndAddBasisSet: every coefficient array is merged into each CGF."""
    out = []
    for sh in bse['elements'][str(Z)]['electron_shells']:
        exps = [float(e) for e in sh['exponents']]
        for l in sh['angular_momentum']:
            prims = []
            for coeffs in sh['coefficients']:
                for k, c in enumerate(coeffs):
                    if float(c) != 0.0: prims.append([exps[k], float(c)])
            out.append([l] + prims)
    return out
geom = "O 0 0 0; H 0.7569520 -0.586349 0; H -0.7569520 -0.586349 0"
ref = scf.RHF(gto.M(atom=geom, basis='sto-3g', verbose=0)).kernel()
bad = scf.RHF(gto.M(atom=geom, basis={'O': scess_basis(8), 'H': scess_basis(1)}, verbose=0)).kernel()
print(f"H2O STO-3G RHF correct {ref:.8f}   as-loaded-by-SCESS {bad:.8f}   diff {bad-ref:+.3e}")
