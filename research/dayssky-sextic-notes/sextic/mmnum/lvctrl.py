import sys, math, numpy as np
from mmdata import Ideals, rows, run_kernel, W
U = int(float(sys.argv[1])); r = float(sys.argv[2])
D = U**r; I = Ideals(int(D)+2); rx, ry = rows(U)
nn = rx*rx - rx*ry + ry*ry
def is3pow(v):
    while v % 3 == 0: v //= 3
    return v == 1
keep = np.array([not is3pow(int(v)) for v in nn]); rx, ry = rx[keep], ry[keep]
rng = np.random.default_rng(3); n = len(I.N)
fp = rng.choice([-1,1], size=len(I.primes)).astype(float)
rm = np.zeros(n); rm[0]=1
for i in range(1,n): rm[i] = fp[I.spf[i]]*rm[I.cof[i]]
rm *= (I.mu != 0)       # random multiplicative supported on squarefree, like mu
lo, hi = I.rng(D/2, D)
specs = [(lo,hi,I.mu[lo:hi]*W(I.N[lo:hi]/D)), (lo,hi,rm[lo:hi]*W(I.N[lo:hi]/D))]
V = run_kernel(I, rx, ry, specs, f"lvc{U}") / math.sqrt(D)
for j,name in enumerate(['mu','rmul']):
    a2 = np.abs(V[:,j])**2; m = a2.mean()
    print(name, ' '.join(f"{t}:{(a2>=t*m).mean()/math.exp(-t):.2f}" for t in (4,8,12,16,20)), f"max/mean={a2.max()/m:.1f}")
