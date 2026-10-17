# Mobius fourth moment numerics: M4(r) = U^{-1} sum_{0<N(u)<=U, u 6pf} |M_u(U^r)|^4, versus Gaussian prediction 2*(M2)^2.
import sys, json, time, math
import numpy as np
from mmdata import Ideals, rows, run_kernel, W
U = int(float(sys.argv[1])); seed = int(sys.argv[2]) if len(sys.argv) > 2 else 7
t0 = time.time()
I = Ideals(U); rx, ry = rows(U); nrows = len(rx)
rng = np.random.default_rng(seed)
sqf = (I.mu != 0)
eps_ind = rng.choice([-1, 1], size=len(I.N)) * sqf
fp = rng.choice([-1, 1], size=len(I.primes))
eps_mul = np.zeros(len(I.N), np.int64); eps_mul[0] = 1
for i in range(1, len(I.N)):
    if sqf[i]: eps_mul[i] = fp[I.spf[i]] * eps_mul[I.cof[i]]
eps_mul *= sqf
COEF = {'mu': I.mu.astype(float), 'rnd': eps_ind.astype(float), 'rmul': eps_mul.astype(float)}
RS = [0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
specs = []; key = []
for r in RS:
    D = U**r; lo, hi = I.rng(D/2, D)
    for v, c in COEF.items():
        specs.append((lo, hi, c[lo:hi] * W(I.N[lo:hi] / D))); key.append((r, v))
V = run_kernel(I, rx, ry, specs, f"m4U{U}")
out = []
for j, (r, v) in enumerate(key):
    M = V[:, j] / math.sqrt(U**r)
    a2 = np.abs(M)**2
    m2 = a2.sum() / U; m4 = (a2**2).sum() / U
    # Gaussian prediction uses the per-row mean square (rows/U * E|M|^2)
    e2 = a2.mean(); g4 = (nrows / U) * 2 * e2**2
    top = np.sort(a2**2)[::-1]
    rec = dict(r=r, var=v, M2=m2, M4=m4, gauss=g4, ratio=m4/g4, top1_share=float(top[0]/a2.__pow__(2).sum()),
               top10_share=float(top[:10].sum()/(a2**2).sum()), max_abs=float(math.sqrt(a2.max())))
    out.append(rec)
    print(f"U={U:.0e} r={r:.2f} {v:4s} M2={m2:.4f} M4={m4:.4f} 2*M2^2(rows)={g4:.4f} ratio={m4/g4:.4f} "
          f"top1={rec['top1_share']:.3f} top10={rec['top10_share']:.3f} max|M|={rec['max_abs']:.2f}", flush=True)
json.dump(dict(U=U, rows=nrows, res=out), open(f"m4_U{U}.json", "w"), indent=1)
print(f"total {time.time()-t0:.1f}s rows={nrows}")
