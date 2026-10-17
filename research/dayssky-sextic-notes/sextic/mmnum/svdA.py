# Singular values of A = (chi_n(u)) for sixth-power-free rows N(u)<=U and squarefree columns N(n) in (D/2, D].
# Sharp sieve <=> s1^2 ~ max(#rows, #cols)*c. BGL allows (UD)^{2/3}. Compare mu's projection.
import sys, math, numpy as np
from mmdata import Ideals, rows, run_kernel
U = int(float(sys.argv[1])); s = float(sys.argv[2])
D = U**s; I = Ideals(int(D) + 2); rx, ry = rows(U)
lo, hi = I.rng(D/2, D); idx = np.arange(lo, hi)[I.mu[lo:hi] != 0]
specs = [(int(i), int(i)+1, np.array([1.0])) for i in idx]
A = run_kernel(I, rx, ry, specs, f"svd{U}")          # rows x cols
nr, nc = A.shape
sv = np.linalg.svd(A, compute_uv=False)
mu = I.mu[idx].astype(float)
Amu = A @ mu
print(f"U={U} s={s} rows={nr} cols={nc}  s1^2={sv[0]**2:.0f} s10^2={sv[9]**2:.0f} median s^2={np.median(sv**2):.0f}  "
      f"max(rows,cols)={max(nr,nc)}  (UD)^(2/3)={(nr*nc)**(2/3):.0f}")
print(f"   ||A mu||^2/(rows*cols) = {np.sum(np.abs(Amu)**2)/(nr*nc):.4f} ;  ||A e||^2/(rows*cols) for random +-1 e: ", end='')
rng = np.random.default_rng(0); e = rng.choice([-1.0, 1.0], nc); print(f"{np.sum(np.abs(A@e)**2)/(nr*nc):.4f}")
print("   top 8 s^2/cols:", np.round(sv[:8]**2/nc, 2))
