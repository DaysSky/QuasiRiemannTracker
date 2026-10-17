import sys, math, numpy as np
from mmdata import Ideals, rows, run_kernel, W
U = int(float(sys.argv[1])); m = float(sys.argv[2])
N = U**m; I = Ideals(int(N)+2); rx, ry = rows(U)
lo, hi = I.rng(N/2, N)
V = run_kernel(I, rx, ry, [(lo, hi, W(I.N[lo:hi]/N))], f"p6top{U}")[:, 0] / math.sqrt(N)
a = np.abs(V); o = np.argsort(-a)[:25]
for i in o:
    x, y = int(rx[i]), int(ry[i]); n = x*x - x*y + y*y
    print(f"u={x:+d}{y:+d}w N={n} |S|={a[i]:.3f}")
