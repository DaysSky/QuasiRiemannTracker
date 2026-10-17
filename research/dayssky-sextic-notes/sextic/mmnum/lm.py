# Long Mobius mean square: MS(s) = (1/#rows) sum_u |D^{-1/2} sum_n mu(n) chi_n(u) W(N(n)/D)|^2, D = U^s, s in (1, 1.5].
# Copies predict growth U^{5(s-1)/6} relative to s=1; LM(alpha'=0) predicts flat. Controls: random +-1 (indep, multiplicative).
import sys, json, time, math
import numpy as np
from mmdata import Ideals, rows, run_kernel, W
U = int(float(sys.argv[1])); smax = float(sys.argv[2]); seed = 7
t0 = time.time()
Dmax = U**smax
I = Ideals(int(Dmax) + 2); rx, ry = rows(U); nrows = len(rx)
print(f"ideals {len(I.N)} rows {nrows} build {time.time()-t0:.1f}s", flush=True)
rng = np.random.default_rng(seed)
sqf = (I.mu != 0)
fp = rng.choice([-1, 1], size=len(I.primes))
eps_mul = np.zeros(len(I.N), np.int64); eps_mul[0] = 1
for i in range(1, len(I.N)):
    if sqf[i]: eps_mul[i] = fp[I.spf[i]] * eps_mul[I.cof[i]]
eps_mul *= sqf
COEF = {'mu': I.mu.astype(float), 'rmul': eps_mul.astype(float), 'sqf': sqf.astype(float)}
SS = [s for s in (0.5, 0.75, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5) if s <= smax + 1e-9]
specs = []; key = []
for s in SS:
    D = U**s; lo, hi = I.rng(D/2, D)
    for v, c in COEF.items():
        specs.append((lo, hi, c[lo:hi] * W(I.N[lo:hi] / D))); key.append((s, v, float((c[lo:hi]**2 * W(I.N[lo:hi]/D)**2).sum()/D)))
V = run_kernel(I, rx, ry, specs, f"lmU{U}")
out = []
for j, (s, v, diag) in enumerate(key):
    a2 = np.abs(V[:, j])**2 / U**s
    rec = dict(s=s, var=v, mean=float(a2.mean()), diag=diag, ratio=float(a2.mean()/diag) if diag > 0 else None,
               max=float(a2.max()))
    out.append(rec)
    print(f"U={U:.0e} s={s:.2f} {v:4s} mean|M|^2/D={rec['mean']:.4f} diag={diag:.4f} ratio={rec['ratio']:.4f} max={rec['max']:.1f}", flush=True)
json.dump(dict(U=U, rows=nrows, res=out), open(f"lm_U{U}.json", "w"), indent=1)
print(f"total {time.time()-t0:.1f}s")
