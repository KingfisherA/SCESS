import numpy as np
from pyscf import gto, scf, ao2mo
mol = gto.M(atom="H 0 0 0; H 0.7122300374 0 0", basis="sto-3g", verbose=0)
mf = scf.RHF(mol).run()
C = mf.mo_coeff; n = C.shape[1]
chem = ao2mo.kernel(mol, C, compact=False).reshape([n]*4)
phys = chem.transpose(0, 2, 1, 3)        # what TransformF12Basis(convertChemToPhys=true) hands over
def anti(t):                             # AntisymmetriseIntegrateOutSpinClosedShell (assumes chemists')
    N = 2*n; sp = np.arange(N) % 2; sa = np.arange(N)//2
    E = t[np.ix_(sa, sa, sa, sa)]
    c = (sp[:,None,None,None]==sp[None,None,:,None])&(sp[None,:,None,None]==sp[None,None,None,:])
    x = (sp[:,None,None,None]==sp[None,None,None,:])&(sp[None,:,None,None]==sp[None,None,:,None])
    return E.transpose(0,2,1,3)*c - E.transpose(0,2,3,1)*x
for name, t in [("correct (chem in)", chem), ("as called (phys in)", phys)]:
    g = anti(t)
    h1 = C.T @ mf.get_hcore() @ C
    H = np.array([[2*h1[0,0] + g[0,1,0,1], g[0,1,2,3]], [g[0,1,2,3], 2*h1[1,1] + g[2,3,2,3]]])
    e = np.linalg.eigvalsh(H)[0] + mol.energy_nuc()
    print(f"{name:20s} <D0|H|D(0a0b->1a1b)> = {g[0,1,2,3]:+.6f}   2x2 CI E = {e:.8f}")
