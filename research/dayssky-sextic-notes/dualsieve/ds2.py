# Rows-dominant sextic sieve with squarefree composite columns n (split primes, squarefree rational norm), direct Gauss sums.
import sys, math, numpy as np
P = int(sys.argv[1]); Ks = [int(float(x)) for x in sys.argv[2].split(',')]
W6 = np.exp(2j*np.pi*np.arange(6)/6)
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
# split primes up to 2P with their data
PR = {}
for p in sieve(2*P):
    p = int(p)
    if p % 3 != 1: continue
    ab = two_sq(p); lst = []
    for x in (ab, (ab[0]-ab[1], -ab[1])):
        pi = primary(x); a, b = pi; wr = (-a*pow(b, -1, p)) % p
        fs = factor(p-1); g = 2
        while any(pow(g, (p-1)//q, p) == 1 for q in fs): g += 1
        ind = np.zeros(p, np.int64); v = 1
        for k in range(p-1): ind[v] = k; v = v*g % p
        gv = pow(g, (p-1)//6, p)
        img = {1: 0, wr % p: 2, wr*wr % p: 4, p-1: 3, (-wr) % p: 5, (-wr*wr) % p: 1}
        lst.append((pi, p, wr, img[gv], ind))
    PR[p] = lst
plist = sorted(PR)
def chi6_rows(data, A, B):
    pi, p, wr, j0, ind = data
    r = (A + B*wr) % p
    return np.where(r == 0, 0, W6[(j0*(ind[r] % 6)) % 6])
# columns: products of primes from distinct rational primes, norm in [P, 2P]
cols = []
def rec(i, prod_el, norm, facs):
    if norm >= P and norm <= 2*P and facs: cols.append((prod_el, norm, list(facs)))
    for k in range(i, len(plist)):
        p = plist[k]
        if norm*p > 2*P: break
        for data in PR[p]:
            facs.append(data); rec(k+1, mul(prod_el, data[0]), norm*p, facs); facs.pop()
rec(0, (1,0), 1, [])
print("columns", len(cols), file=sys.stderr)
# Gauss sums g6(n) = sum_{d mod N} chi6_n(d) e(d*Tr(n)/N), chi6_n(d) = prod (d/pi_i)_6 for rational d
cvals = []
for (n, N, facs) in cols:
    d = np.arange(N); ch = np.ones(N, complex)
    for data in facs:
        pi, p, wr, j0, ind = data
        r = d % p
        ch *= np.where(r == 0, 0, W6[(j0*(ind[r] % 6)) % 6])
    T = (2*n[0] - n[1]) % N
    g = np.sum(ch*np.exp(2j*np.pi*(d*T % N)/N))/math.sqrt(N)
    mu = (-1)**len(facs)
    cvals.append((mu, g))
mu = np.array([c[0] for c in cvals]); g6 = np.array([c[1] for c in cvals])
print("mean|g6|", np.mean(np.abs(g6)), file=sys.stderr)
rng = np.random.default_rng(2)
seqs = {'mu*g6 (LM dual)': mu*g6, 'g6 (no Mobius)': g6, 'random': np.exp(2j*np.pi*rng.random(len(cols)))}
Kmax = max(Ks); R = int(math.isqrt(4*Kmax//3)) + 3
A, B = np.meshgrid(np.arange(-R, R+1), np.arange(-R, R+1), indexing='ij'); A = A.ravel(); B = B.ravel()
nr = A*A - A*B + B*B
m = (nr > 0) & (nr <= Kmax) & (A % 3 == 1) & (B % 3 == 0)
A, B, nr = A[m], B[m], nr[m]; o = np.argsort(nr); A, B, nr = A[o], B[o], nr[o]
print("rows", len(A), file=sys.stderr)
Tacc = {k: np.zeros(len(A), complex) for k in seqs}
for j, (n, N, facs) in enumerate(cols):
    x = np.ones(len(A), complex)
    for data in facs: x *= np.conj(chi6_rows(data, A, B))
    for k, c in seqs.items(): Tacc[k] += c[j]*x
for k, c in seqs.items():
    cum = np.cumsum(np.abs(Tacc[k])**2); nc = np.sum(np.abs(c)**2)
    out = []
    for K in Ks:
        kk = np.searchsorted(nr, K, side='right'); out.append(f"{K}:{cum[kk-1]/(kk*nc):.4f}")
    print(f"{k:18s}", ' '.join(out), flush=True)
