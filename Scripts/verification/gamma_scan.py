import numpy as np, mp2f12_port as P, f12_ints as FI
geom = "H 0 0 0; H 0.7122300374 0 0"
mo, both = FI.build_mol(geom, "sto-3g", "cc-pvdz")
base = P.setup(geom)
gi = FI.GeminalInts(both)
for gam in (0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4):
    cg = [(a * gam**2, -c / gam) for a, c in FI.STG6G]   # exp(-gam r) fit, f = -exp(-gam r)/gam
    d = dict(base); d["f12"], d["f12sq"], d["f12g"], d["delsq"] = gi.all(cg)
    m = P.transform_f12_basis(d); I = P.intermediates(m, trace_bug=False)
    print(f"gamma={gam:.1f}  fixed={P.e_fixed(I, 1):.9f}  opt={P.e_opt(m, I):.9f}")
