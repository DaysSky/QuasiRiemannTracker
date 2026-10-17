# direct check of twisted multiplicativity: compute g(n) by brute force for small composite n and compare
import math, cmath, numpy as np, sys
sys.argv = ['x', '400']
exec(open('gbias.py').read().split("sq = I.mu != 0")[0])
from eis import chi
def g_direct(n):
    Nn = norm(n); c0, c1 = tr_inv(n)
    # residues mod n: x0 + x1 w with 0<=x0,x1<Nn, reduce duplicates by canonical form via ideal lattice: brute force over box and dedupe
    seen = {}; s = 0
    for x0 in range(Nn):
        for x1 in range(Nn):
            # canonical rep: (x0,x1) mod n -> use mul by conj(n) mod Nn
            p = mul((x0, x1), conj(n)); key = (p[0] % Nn, p[1] % Nn)
            if key in seen: continue
            seen[key] = 1
            k = chi(n, (x0, x1), PR)
            if k == 6: continue
            s += cmath.exp(1j*math.pi*k/3) * cmath.exp(2j*math.pi*((x0*c0 + x1*c1) % 1.0))
    assert len(seen) == Nn, (len(seen), Nn)
    return s / math.sqrt(Nn)
cnt = 0
for i in range(1, len(I.N)):
    if I.mu[i] == 0 or I.cof[i] == 0: continue
    if I.N[i] > 400: break
    gd = g_direct(I.g[i]); err = abs(gd - gt[i])
    print(I.g[i], I.N[i], round(gd.real,4), round(gd.imag,4), round(gt[i].real,4), round(gt[i].imag,4), "ERR" if err > 1e-6 else "ok")
    cnt += 1
    if cnt > 12: break
