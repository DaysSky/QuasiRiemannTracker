# Thin-family plain fourth moment (input 2, DH_SEXTIC §98): mean of |S1 S2(u)|^2 over rows u = 0 mod m^2, against the mean over all rows.
# Natural size (ratio ~ 1) is cost 2; Lemma 18.1's twist rule corresponds to ratio ~ N(m); positivity to ratio ~ N(m)^2.
import sys, numpy as np
from mmdata import Ideals, rows, run_kernel, W, mul, norm
U = int(float(sys.argv[1])); n = float(sys.argv[2])
X = U**n
I = Ideals(int(X) + 2); rx, ry = rows(U)
lo, hi = I.rng(X/2, X); w = W(I.N[lo:hi]/X)
V = run_kernel(I, rx, ry, [(lo, hi, w)], f"thin{U}_{n}")
S = V[:, 0] / np.sqrt(X); A = np.abs(S)**4
meanA = A.mean(); rng = np.random.default_rng(1)
print(f"U={U:.0e} X=U^{n}: rows={len(rx)} mean|S1S2|^2={meanA:.4f}")
for pr in I.primes:
    P = pr['P']
    if P > 200: break
    g = pr['g']; g2 = mul(g, g); Ng = norm(g2); cg = (g2[0] - g2[1], -g2[1])
    pa = rx*cg[0] - ry*cg[1]; pb = rx*cg[1] + ry*cg[0] - ry*cg[1]
    sel = (pa % Ng == 0) & (pb % Ng == 0); k = int(sel.sum())
    if k < 30: continue
    ratio = A[sel].mean()/meanA
    ctrl = np.median([A[rng.choice(len(A), k, replace=False)].mean()/meanA for _ in range(50)])
    sd = np.std([A[rng.choice(len(A), k, replace=False)].mean()/meanA for _ in range(50)])
    print(f"  N(m)={P:4d}  rows with m^2|u: {k:6d}  ratio={ratio:.3f}   random-subset control {ctrl:.3f} +- {sd:.3f}   (twist-rule loss would be ~{P})")
