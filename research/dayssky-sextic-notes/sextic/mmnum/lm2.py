# LM for 1/L^2: coefficients mu2 = mu*mu (multiplicative: p -> -2, p^2 -> 1, p^k -> 0 for k>=3). Ratio to diagonal.
import sys, math, numpy as np
from mmdata import Ideals, rows, run_kernel, W
U = int(float(sys.argv[1])); smax = float(sys.argv[2])
I = Ideals(int(U**smax) + 2); rx, ry = rows(U)
n = len(I.N); m2 = np.zeros(n); m2[0] = 1
# exponent of spf in each ideal: track via cof chain
for i in range(1, n):
    c = I.cof[i]; p = I.spf[i]
    # count multiplicity of p in i
    k = 1; cc = c
    while cc > 0 and I.spf[cc] == p: k += 1; cc = I.cof[cc]
    val = {1: -2.0, 2: 1.0}.get(k, 0.0)
    m2[i] = val * m2[cc]
rng = np.random.default_rng(13); fp = rng.choice([-1, 1], size=len(I.primes))
# random multiplicative analogue of mu*mu: f(p) = 2*eps_p ... use (eps*eps): p -> 2 eps_p, p^2 -> 1
r2 = np.zeros(n); r2[0] = 1
for i in range(1, n):
    c = I.cof[i]; p = I.spf[i]; k = 1; cc = c
    while cc > 0 and I.spf[cc] == p: k += 1; cc = I.cof[cc]
    val = {1: 2.0*fp[p], 2: 1.0}.get(k, 0.0)
    r2[i] = val * r2[cc]
SS = [s for s in (0.5, 0.75, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5) if s <= smax + 1e-9]
specs = []; key = []
for s in SS:
    D = U**s; lo, hi = I.rng(D/2, D)
    for lab, c in (('mu*mu', m2), ('rmul2', r2)):
        w = c[lo:hi]*W(I.N[lo:hi]/D); specs.append((lo, hi, w)); key.append((s, lab, float((w**2).sum())))
V = run_kernel(I, rx, ry, specs, f"lm2{U}")
for j, (s, lab, d) in enumerate(key):
    a2 = np.abs(V[:, j])**2
    print(f"U={U:.0e} s={s:.2f} {lab:6s} mean|.|^2/diag = {a2.mean()/d:.4f}   max/mean = {a2.max()/a2.mean():.1f}", flush=True)
