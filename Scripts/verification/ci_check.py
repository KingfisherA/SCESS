import itertools, numpy as np
from pyscf import gto, scf, fci, ao2mo
mol = gto.M(atom="Li 0 0 0; H 1.6 0 0", basis="sto-3g", verbose=0)
mf = scf.RHF(mol).run(conv_tol=1e-12)
C = mf.mo_coeff; nmo = C.shape[1]; ne = mol.nelectron
eri = ao2mo.kernel(mol, C, compact=False).reshape([nmo]*4)
h = C.T @ mf.get_hcore() @ C
n = 2*nmo; sp = np.arange(n)%2; sa = np.arange(n)//2
E = eri[np.ix_(sa,sa,sa,sa)]
coul = (sp[:,None,None,None]==sp[None,None,:,None])&(sp[None,:,None,None]==sp[None,None,None,:])
exch = (sp[:,None,None,None]==sp[None,None,None,:])&(sp[None,:,None,None]==sp[None,None,:,None])
g = E.transpose(0,2,1,3)*coul - E.transpose(0,2,3,1)*exch
H1 = np.where((np.arange(n)[:,None]+np.arange(n)[None,:])%2==0, h[np.ix_(sa,sa)], 0)

def cnt(x, nb): return sum((x>>i)&1 for i in range(nb))
def diffs(I, J, norb):
    nel = cnt(I[0],norb)+cnt(I[1],norb)
    exl = cnt(I[0]^J[0],norb)+cnt(I[1]^J[1],norb)
    shift = nel-exl; flips=0; ui=[]; uj=[]; cm=[]
    for s,(a,b) in enumerate(zip(I,J)):
        iel1=iel2=ie1=ie2=0
        for i in range(norb):
            if (b>>i)&1: iel2+=1
            if (a>>i)&1:
                iel1+=1
                if (b>>i)&1: cm.append(2*i+s)
                else: ie1+=1; flips+=shift-iel1+ie1; ui.append(2*i+s)
            elif (b>>i)&1: ie2+=1; flips+=shift-iel2+ie2; uj.append(2*i+s)
    return ui,uj,cm,(-1 if flips%2 else 1)
def Hel(I,J):
    ui,uj,cm,s = diffs(I,J,nmo)
    if len(ui)==0: return sum(H1[m,m] for m in cm)+sum(g[cm[a],cm[b],cm[a],cm[b]] for a in range(len(cm)) for b in range(a+1,len(cm)))
    if len(ui)==1: return s*(H1[ui[0],uj[0]]+sum(g[ui[0],k,uj[0],k] for k in cm))
    if len(ui)==2: return s*g[ui[0],ui[1],uj[0],uj[1]]
    return 0.0
dets=[]
for comb in itertools.combinations(range(n), ne):
    a=b=0
    for o in comb:
        if o%2==0: a|=1<<(o//2)
        else: b|=1<<(o//2)
    dets.append((a,b))
H=np.array([[Hel(I,J) for J in dets] for I in dets])
print("asymmetry", np.abs(H-H.T).max())
e_code = np.linalg.eigvalsh(H)[0]+mol.energy_nuc()
e_ref = fci.FCI(mf).kernel()[0]
print(f"SCESS-style FCI {e_code:.10f}  PySCF FCI {e_ref:.10f}  diff {e_code-e_ref:+.2e}")
