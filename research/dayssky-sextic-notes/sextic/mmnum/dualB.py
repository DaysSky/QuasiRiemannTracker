# Dual side of the LM involution: B_h(D) = sum_n mu(n) gtilde(n) chi_n(h) W(N(n)/D) over ALL rows h with N(h) <= K = D^2/U.
# Distribution of |B_h|^2/sum|w|^2; which h carry the excess (sixth-power-ish h, small h)?
import sys, math, cmath, numpy as np
from mmdata import Ideals, rows, run_kernel, W
from eis import chi_prime, norm
U = int(float(sys.argv[1])); s = float(sys.argv[2])
D = U**s; K = D*D/U
I = Ideals(int(D) + 2)
sys.argv = ['g', str(int(D) + 2)]
src = open('gbias.py').read().split("sq = I.mu != 0")[0]
ns = {}; exec(src, ns); gt = ns['gt']; I2 = ns['I']
assert len(I2.N) == len(I.N)
lo, hi = I.rng(D/2, D)
w = I.mu[lo:hi]*gt[lo:hi]*W(I.N[lo:hi]/D)
hx, hy = rows(int(K), sixth_power_free=False)
V = run_kernel(I, hx, hy, [(lo, hi, w.real), (lo, hi, w.imag)], f"dB{U}")
B = V[:, 0] + 1j*V[:, 1]
a2 = np.abs(B)**2 / float((np.abs(w)**2).sum())
Nh = hx*hx - hx*hy + hy*hy
print(f"U={U:.0e} s={s} D={D:.0f} K={K:.0f} rows h={len(hx)} mean|B|^2/diag = {a2.mean():.4f}  max = {a2.max():.1f}")
order = np.argsort(-a2)
for j in order[:8]:
    print(f"   h=({hx[j]},{hy[j]}) N(h)={Nh[j]} |B|^2/diag={a2[j]:.1f}")
for lo_, hi_ in ((0, K**0.25), (K**0.25, K**0.5), (K**0.5, K)):
    m = (Nh > lo_) & (Nh <= hi_)
    print(f"   N(h) in ({lo_:.0f},{hi_:.0f}]: count {m.sum()} mean {a2[m].mean():.3f}")
