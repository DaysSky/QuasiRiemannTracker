# Moments of the centrally normalised plain block S_u(X) = X^{-1/2} sum_{n} chi_n(u) W(N(n)/X), X = c*U^{1/2}, over all sixth-power-free rows N(u) <= U.
# Input 4 (HIGH-side sixth moment) needs mean|S|^6 << U^{lambda} * mean|S|^2^3 with lambda < 0.20; copies give 0.417. Random model: m4/m2^2 = 2, m6/m2^3 = 6.
import sys, time, numpy as np
from mmdata import Ideals, rows, run_kernel, W
U = int(float(sys.argv[1])); c = float(sys.argv[2]) if len(sys.argv) > 2 else 1.0
t0 = time.time(); X = c*U**0.5
I = Ideals(int(X) + 2); rx, ry = rows(U)
lo, hi = I.rng(X/2, X)
rng = np.random.default_rng(3)
fp = rng.choice([-1.0, 1.0], size=len(I.primes))   # random completely multiplicative +-1 control
cm = np.ones(len(I.N)); 
for i in range(1, len(I.N)): cm[i] = fp[I.spf[i]]*cm[I.cof[i]]
w = W(I.N[lo:hi]/X)
specs = [(lo, hi, w), (lo, hi, w*cm[lo:hi])]
V = run_kernel(I, rx, ry, specs, f"six{U}") / np.sqrt(X)
for j, name in enumerate(['plain', 'rand-mult']):
    a = np.abs(V[:, j])**2
    m2, m4, m6 = a.mean(), (a**2).mean(), (a**3).mean()
    print(f"U={U:.0e} X={X:.0f} rows={len(rx)} {name:9s} m2={m2:.4f} m4/m2^2={m4/m2**2:.3f} m6/m2^3={m6/m2**3:.3f} max|S|^2/m2={a.max()/m2:.1f}", flush=True)
print(f"time {time.time()-t0:.1f}s")
