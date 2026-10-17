# Rows-dominant sextic large sieve test: S(K) = sum_{N(h)<=K} |sum_{pi: N(pi)~P} c_pi conj((h/pi)_6)|^2 vs (#h) sum|c|^2
import sys, math, cmath, numpy as np
P = int(sys.argv[1]); Ks = [int(float(x)) for x in sys.argv[2].split(',')]; seed = 1
W6 = np.exp(2j*np.pi*np.arange(6)/6); W = cmath.exp(2j*math.pi/3)
def mul(x, y):
    a, b = x; c, d = y
    return (a*c - b*d, a*d + b*c - b*d)
U6 = [(1,0),(-1,0),(0,1),(0,-1),(1,1),(-1,-1)]
def primary(x):
    for u in U6:
        y = mul(u, x)
        if y[0] % 3 == 1 and y[1] % 3 == 0: return y
def sieve(n):
    s = np.ones(n+1, bool); s[:2] = False
    for i in range(2, int(n**.5)+1):
        if s[i]: s[i*i::i] = False
    return np.nonzero(s)[0]
def two_sq(p):
    for b in range(1, int(math.isqrt(4*p//3))+2):
        d = 4*p - 3*b*b
        if d < 0: break
        r = math.isqrt(d)
        if r*r == d and (b + r) % 2 == 0: return ((b+r)//2, b)
def factor(n):
    f=set(); d=2
    while d*d<=n:
        while n%d==0: f.add(d); n//=d
        d+=1
    if n>1: f.add(n)
    return f
cols = []
for p in sieve(2*P):
    p = int(p)
    if p < P or p % 3 != 1: continue
    ab = two_sq(p)
    for x in (ab, (ab[0]-ab[1], -ab[1])):
        pi = primary(x); a, b = pi; wr = (-a*pow(b, -1, p)) % p
        fs = factor(p-1); g = 2
        while any(pow(g, (p-1)//q, p) == 1 for q in fs): g += 1
        ind = np.zeros(p, np.int64); v = 1
        for k in range(p-1): ind[v] = k; v = v*g % p
        # identify zeta0 = image of g^{(p-1)/6} among sixth roots: sixth roots in F_p are +-1, +-wr, +-wr^2 <-> W6 powers
        gv = pow(g, (p-1)//6, p)
        # W6[j] = e^{2 pi i j/6}; w = e^{2pi i/3} = W6[2]; -1 = W6[3]
        img = {1: 0, wr % p: 2, wr*wr % p: 4, p-1: 3, (-wr) % p: 5, (-wr*wr) % p: 1}
        j0 = img[gv]
        chi6 = lambda r: np.where(r % p == 0, 0, W6[(j0*(ind[r % p] % 6)) % 6])
        # sextic Gauss sum g6(pi) = sum_d chi6(d) e(Tr(d/pi)), Tr(d/pi) = d*T/p, T = Tr(pi) = 2a - b
        d = np.arange(1, p); T = (2*a - b) % p
        g6 = np.sum(chi6(d)*np.exp(2j*np.pi*(d*T % p)/p))/math.sqrt(p)
        chi3 = lambda r: chi6(r)**2
        g3 = np.sum(chi3(d)*np.exp(2j*np.pi*(d*T % p)/p))/math.sqrt(p)
        cols.append((pi, p, wr, j0, ind, g6, g3))
print("columns", len(cols), "P", P, file=sys.stderr)
rng = np.random.default_rng(seed)
seqs = {'mu*g6 (LM dual)': np.array([-c[5] for c in cols]),
        'conj g3 (biased)': np.array([np.conj(c[6]) for c in cols]),
        'random': np.exp(2j*np.pi*rng.random(len(cols)))}
Kmax = max(Ks)
# rows: primary h = a + b w with N(h) <= Kmax
R = int(math.isqrt(4*Kmax//3)) + 3
A, Bv = np.meshgrid(np.arange(-R, R+1), np.arange(-R, R+1), indexing='ij'); A = A.ravel(); Bv = Bv.ravel()
nr = A*A - A*Bv + Bv*Bv
m = (nr > 0) & (nr <= Kmax) & (A % 3 == 1) & (Bv % 3 == 0)
A, Bv, nr = A[m], Bv[m], nr[m]
order = np.argsort(nr); A, Bv, nr = A[order], Bv[order], nr[order]
print("rows", len(A), file=sys.stderr)
X = np.zeros((len(A), len(cols)), complex)
for j, (pi, p, wr, j0, ind, g6, g3) in enumerate(cols):
    r = (A + Bv*wr) % p
    X[:, j] = np.where(r == 0, 0, np.conj(W6[(j0*(ind[r] % 6)) % 6]))
for name, c in seqs.items():
    T = X @ c; T2 = np.abs(T)**2; cum = np.cumsum(T2)
    out = []
    for K in Ks:
        k = np.searchsorted(nr, K, side='right')
        out.append(f"K={K}: rows={k} ratio={cum[k-1]/(k*np.sum(np.abs(c)**2)):.4f}")
    print(name, '|', '  '.join(out), flush=True)
