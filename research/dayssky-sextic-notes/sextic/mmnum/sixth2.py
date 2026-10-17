import sys, time, numpy as np
from mmdata import Ideals, rows, run_kernel, W
U = int(float(sys.argv[1])); c = float(sys.argv[2]) if len(sys.argv) > 2 else 1.0
X = c*U**0.5; I = Ideals(int(X) + 2); rx, ry = rows(U)
lo, hi = I.rng(X/2, X); w = W(I.N[lo:hi]/X)
V = run_kernel(I, rx, ry, [(lo, hi, w)], f"sx{U}")[:, 0] / np.sqrt(X)
a = np.abs(V)**2; Nu = rx*rx - rx*ry + ry*ry
def smooth23(n):
    n = int(n)
    for p in (2, 3):
        while n % p == 0: n //= p
    return n == 1
top = np.argsort(-a)[:12]
print("largest rows: (N(u), |S|^2/m2, 2-3-smooth norm?)")
print([(int(Nu[i]), round(a[i]/a.mean(), 1), smooth23(Nu[i])) for i in top])
theta = np.array([smooth23(n) for n in Nu])
for name, mask in (('all', np.ones(len(a), bool)), ('drop {2,3}-smooth norms', ~theta), ('also drop N(u)<1000', (~theta) & (Nu >= 1000))):
    b = a[mask]; m2 = b.mean()
    print(f"U={U:.0e} X={X:.0f} {name:24s} rows={mask.sum()} m4/m2^2={(b**2).mean()/m2**2:.3f} m6/m2^3={(b**3).mean()/m2**3:.3f} max/m2={b.max()/m2:.1f}")
