# Off-diagonal of the LM dual sieve split by the square class of the row h (h = h0 * k^2, h0 squarefree up to units).
# Rows with h0 a unit (pure squares) have no quadratic twist; RS applies there. Question: where does the off-diagonal live?
import sys, math, numpy as np, collections
exec(open('ds2.py').read().split('rng = np.random.default_rng(2)')[0].replace('P = int(sys.argv[1]); Ks = [int(float(x)) for x in sys.argv[2].split(\',\')]','P = int(sys.argv[1]); Ks=[int(sys.argv[2])]'))
Kmax = Ks[0]; R = int(math.isqrt(4*Kmax//3)) + 3
A, B = np.meshgrid(np.arange(-R, R+1), np.arange(-R, R+1), indexing='ij'); A = A.ravel(); B = B.ravel()
nr = A*A - A*B + B*B
m = (nr > 0) & (nr <= Kmax) & (A % 3 == 1) & (B % 3 == 0)
A, B, nr = A[m], B[m], nr[m]
c = mu*g6
T = np.zeros(len(A), complex); diag = np.zeros(len(A))
for j, (n, N, facs) in enumerate(cols):
    x = np.ones(len(A), complex)
    for data in facs: x *= np.conj(chi6_rows(data, A, B))
    T += c[j]*x; diag += np.abs(c[j]*x)**2
off = np.abs(T)**2 - diag
# square class: squarefree kernel of the rational norm is a proxy; exact: h is a square times unit iff all prime valuations even
def sqclass(a, b):
    nrm = a*a - a*b + b*b; f = 1; d = 2; x = nrm
    while d*d <= x:
        e = 0
        while x % d == 0: x//=d; e+=1
        if e % 2: f *= d
        d += 1
    if x > 1: f *= x
    return f
cls = np.array([sqclass(int(a), int(b)) for a, b in zip(A, B)])
print("rows", len(A), "sum diag", diag.sum(), "sum off", off.sum(), "rel", off.sum()/diag.sum())
sq = cls == 1
print(f"norm-square rows: {sq.sum()}  off/diag there {off[sq].sum()/diag[sq].sum():+.4f}   elsewhere {off[~sq].sum()/diag[~sq].sum():+.4f}")
small = cls <= 50
print(f"rows with norm square-class <= 50: {small.sum()}  off/diag {off[small].sum()/diag[small].sum():+.4f}")
# dispersion of class means
g = collections.defaultdict(list)
for o, d, k in zip(off, diag, cls): g[k].append((o, d))
big = [(k, sum(o for o,_ in v)/sum(d for _,d in v), len(v)) for k, v in g.items() if len(v) >= 8]
vals = np.array([b[1] for b in big]); print("classes with >=8 rows:", len(big), "mean of class off/diag", vals.mean().round(4), "sd", vals.std().round(4))
