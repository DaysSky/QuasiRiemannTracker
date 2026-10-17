# Twisted fourth moments T(m,m') = sum_u |S1 S2(u)|^2 chi_m(u) conj(chi_m'(u)) for prime m != m' (slot primes),
# with S1 = S2 the centrally normalised plain block at length X = U^n (2n + z <= 1 region). Compare |T| / sum|S1S2|^2 with N(mm')^{-1/2}.
import sys, numpy as np
from mmdata import Ideals, rows, run_kernel, W
U = int(float(sys.argv[1])); n = float(sys.argv[2]); zlen = float(sys.argv[3])
X = U**n; Zp = U**zlen
I = Ideals(int(max(X, Zp)) + 2); rx, ry = rows(U)
lo, hi = I.rng(X/2, X); w = W(I.N[lo:hi]/X)
# prime ideals with norm in [Zp/2, Zp]: pick split primes from the Ideals table

pidx = [i for i in range(*I.rng(Zp/2, Zp)) if i > 0 and I.cof[i] == 0]
print("prime slots found", len(pidx), file=sys.stderr)
specs = [(lo, hi, w)] + [(p, p+1, np.ones(1)) for p in pidx[:40]]
V = run_kernel(I, rx, ry, specs, f"tw{U}")
S = V[:, 0] / np.sqrt(X); A = np.abs(S)**4          # |S1 S2|^2 with S1 = S2
chis = V[:, 1:]                                     # chi_p(u) values (unit modulus or 0)
base = A.sum(); res = []
for i in range(chis.shape[1]):
    for j in range(i+1, chis.shape[1]):
        T = (A * chis[:, i] * np.conj(chis[:, j])).sum()
        Nmm = I.N[pidx[i]] * I.N[pidx[j]]
        res.append((abs(T)/base, Nmm**-0.5))
r = np.array(res)
print(f"U={U:.0e} X=U^{n} slots~U^{zlen} pairs={len(r)}  median |T|/base={np.median(r[:,0]):.4f}  median N(mm')^-1/2={np.median(r[:,1]):.4f}  "
      f"fitted theta'={np.median(np.log(r[:,0])/np.log(1/r[:,1])) if len(r) else float('nan'):.3f}  rows-noise 1/sqrt(rows)~{1/np.sqrt(len(rx)):.4f}")
