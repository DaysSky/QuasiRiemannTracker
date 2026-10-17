# Large-value structure of Mobius polynomials M_u(D) over the sextic family:
# count rows with |M_u|^2 >= V * mean, compare with Gaussian prediction exp(-V); check norm-residue / small-prime structure of large rows.
import sys, math, numpy as np
from mmdata import Ideals, rows, run_kernel, W
U = int(float(sys.argv[1])); r = float(sys.argv[2])
D = U**r; I = Ideals(int(D)+2); rx, ry = rows(U)
nn = rx*rx - rx*ry + ry*ry
def is3pow(v):
    while v % 3 == 0: v //= 3
    return v == 1
keep = np.array([not is3pow(int(v)) for v in nn]); rx, ry, nn = rx[keep], ry[keep], nn[keep]
lo, hi = I.rng(D/2, D)
specs = [(lo, hi, I.mu[lo:hi]*W(I.N[lo:hi]/D))]
V = run_kernel(I, rx, ry, specs, f"lv{U}")[:,0] / math.sqrt(D)
a2 = np.abs(V)**2; mu2 = a2.mean()
print(f"U={U:.0e} r={r} rows={len(a2)} mean|M|^2={mu2:.4f}")
for t in (2,4,6,8,10,12,14):
    obs = (a2 >= t*mu2).mean(); print(f"  P(|M|^2 >= {t:2d} mean) = {obs:.3e}   gaussian {math.exp(-t):.3e}  ratio {obs/math.exp(-t):.2f}")
# structure: fraction of large rows with small norm prime factors / smooth norms
big = a2 >= 8*mu2
def smallfac(v, B=50):
    c=0
    for p in (2,3,5,7,11,13,17,19,23,29,31,37,41,43,47):
        if v % p == 0: c+=1
    return c
sf_all = np.array([smallfac(int(v)) for v in nn[::max(1,len(nn)//20000)]]).mean()
sf_big = np.array([smallfac(int(v)) for v in nn[big]]).mean() if big.any() else float('nan')
print(f"  mean # small rational prime factors of N(u): all {sf_all:.3f}, large rows {sf_big:.3f} (n_big={big.sum()})")
