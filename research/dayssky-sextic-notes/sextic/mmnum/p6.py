# Plain sixth moment numerics for the sextic family: S_u(N) = N^{-1/2} sum_n chi_n(u) W(N(n)/N),
# rows u sixth-power-free with 0 < N(u) <= U.  Reports M2, M4/(2 M2^2), M6/(6 M2^3) (complex-Gaussian ratios),
# and mixed sixth moments |S(N1)S(N2)S(N3)|^2.  Controls: random multiplicative and independent random signs.
import sys, json, time, math
import numpy as np
from mmdata import Ideals, rows, run_kernel, W
U = int(float(sys.argv[1])); seed = int(sys.argv[2]) if len(sys.argv) > 2 else 7
MS = [float(x) for x in (sys.argv[3].split(',') if len(sys.argv) > 3 else ['0.3','0.35','0.4','0.45','0.5'])]
t0 = time.time()
Xmax = max(U**m for m in MS)
I = Ideals(int(Xmax) + 2); rx, ry = rows(U); nrows = len(rx)
rng = np.random.default_rng(seed)
n = len(I.N)
fp = rng.choice([-1, 1], size=len(I.primes)).astype(float)
rmul = np.zeros(n); rmul[0] = 1
for i in range(1, n):
    rmul[i] = fp[I.spf[i]] * rmul[I.cof[i]]          # completely multiplicative random signs
rind = rng.choice([-1, 1], size=n).astype(float)
COEF = {'one': np.ones(n), 'rmul': rmul, 'rnd': rind}
specs = []; key = []
for m in MS:
    N = U**m; lo, hi = I.rng(N/2, N)
    for v, c in COEF.items():
        specs.append((lo, hi, c[lo:hi] * W(I.N[lo:hi] / N))); key.append((m, v))
V = run_kernel(I, rx, ry, specs, f"p6U{U}")
res = []
cols = {}
for j, (m, v) in enumerate(key):
    S = V[:, j] / math.sqrt(U**m); cols[(m, v)] = S
    a2 = np.abs(S)**2
    M2 = a2.mean(); M4 = (a2**2).mean(); M6 = (a2**3).mean()
    rec = dict(m=m, var=v, M2=M2, r4=M4/(2*M2**2), r6=M6/(6*M2**3), max=float(np.sqrt(a2.max())),
               top1_6=float((a2**3).max()/(a2**3).sum()))
    res.append(rec)
    print(f"U={U:.0e} m={m:.2f} {v:4s} M2={M2:.4f} r4={rec['r4']:.4f} r6={rec['r6']:.4f} max={rec['max']:.2f} top1(6)={rec['top1_6']:.4f}", flush=True)
# mixed sixth moments for (m,m,m') combos of the 'one' family
for v in COEF:
    for (m1, m2, m3) in [(0.3, 0.4, 0.5), (0.4, 0.45, 0.5), (0.45, 0.45, 0.45), (0.5, 0.5, 0.5), (0.35, 0.45, 0.5)]:
        if all((mm, v) in cols for mm in (m1, m2, m3)):
            P = np.abs(cols[(m1, v)] * cols[(m2, v)] * cols[(m3, v)])**2
            g = np.prod([np.mean(np.abs(cols[(mm, v)])**2) for mm in (m1, m2, m3)])
            print(f"   mixed {v:4s} ({m1},{m2},{m3}): mean|S1S2S3|^2 / prod M2 = {P.mean()/g:.4f}", flush=True)
json.dump(dict(U=U, rows=nrows, res=res), open(f"p6_U{U}.json", "w"), indent=1)
print(f"total {time.time()-t0:.1f}s rows={nrows}")
