# Cross-scale correlation of Mobius large values: rows' |M(U^r)|^2/mean at r in {0.745, 1.0, 1.162}; random multiplicative control.
import sys, math, numpy as np
from mmdata import Ideals, rows, run_kernel, W
U = int(float(sys.argv[1])); RS = [0.5, 0.745, 1.0, 1.162]
I = Ideals(int(U**max(RS)) + 2); rx, ry = rows(U)
rng = np.random.default_rng(11); sqf = I.mu != 0
fp = rng.choice([-1, 1], size=len(I.primes)); em = np.zeros(len(I.N), np.int64); em[0] = 1
for i in range(1, len(I.N)):
    if sqf[i]: em[i] = fp[I.spf[i]]*em[I.cof[i]]
specs = []
for c in (I.mu.astype(float), em.astype(float)):
    for r in RS:
        D = U**r; lo, hi = I.rng(D/2, D); specs.append((lo, hi, c[lo:hi]*W(I.N[lo:hi]/D)))
V = run_kernel(I, rx, ry, specs, f"xs{U}")
for k, lab in enumerate(('mu', 'rmul')):
    A = np.abs(V[:, 4*k:4*k+4])**2; A = A/A.mean(axis=0)
    L = np.log(A + 1e-12)
    C = np.corrcoef(L.T)
    print(f"{lab}: corr(log|M|^2) between scales {RS}:\n{np.array2string(C, precision=3)}")
    # adversary-type rows: top 1% at r*=0.745 and t*=1.162, what is their |M(1)|^2/mean?
    thr = np.quantile(A, 0.99, axis=0)
    both = (A[:, 1] >= thr[1]) & (A[:, 3] >= thr[3])
    print(f"   rows top1% at both 0.745 and 1.162: {both.sum()} (indep. expect {len(A)*1e-4:.1f}); their mean |M(U^1)|^2/mean = {A[both,2].mean() if both.any() else float('nan'):.2f}; frac with |M(U^1)|^2 < median: {(A[both,2] < np.median(A[:,2])).mean() if both.any() else float('nan'):.2f}")
