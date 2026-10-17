# cross-check C kernel against the reference symbol implementation
import numpy as np, random
from mmdata import *
from eis import chi, zeta6, ZERO
I = Ideals(3000)
print("ideals", len(I.N), "primes", len(I.primes), "count check N<=3000 primary odd ideals ~0.302X:", len(I.N)/3000)
rx, ry = rows(3000)
print("rows", len(rx), "ratio to 2pi/sqrt3 U:", len(rx)/(2*np.pi/np.sqrt(3)*3000))
# spec: each single ideal i -> weight 1  (to read off chi_{n_i}(u) individually) for a sample of ideals
random.seed(2)
samp = sorted(random.sample(range(1, len(I.N)), 60))
specs = [(i, i+1, np.ones(1)) for i in samp]
sel = np.random.default_rng(0).choice(len(rx), 300, replace=False)
V = run_kernel(I, rx[sel], ry[sel], specs, 'test')
bad = 0
for a, r in enumerate(sel):
    u = (int(rx[r]), int(ry[r]))
    for b, i in enumerate(samp):
        ref = zeta6(chi(I.g[i], u, I.primes))
        if abs(V[a, b] - ref) > 1e-5: bad += 1
print("kernel vs reference mismatches:", bad, "of", len(sel)*len(samp))
# mu check: compare to factorization
mm = 0
for i in samp:
    f = I.factor(i); sq = len(set(f)) == len(f)
    mm += (I.mu[i] != ((-1)**len(f) if sq else 0))
print("mu mismatches:", mm)
