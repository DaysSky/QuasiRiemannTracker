# Is the normalised quadratic Gauss sum g2(n) over Z[w] a function of n mod 4 (n squarefree, primary, split primes)?
# Also test eps(n) = mu(n) g6(n) / ((n/|n|) conj g3(n)) for periodicity mod 4*9 etc.
import sys, math, numpy as np, itertools, collections
sys.path.insert(0, '../dualsieve')
P = int(sys.argv[1]) if len(sys.argv) > 1 else 400
W6 = np.exp(2j*np.pi*np.arange(6)/6)
def mul(x, y):
    a, b = x; c, d = y
    return (a*c - b*d, a*d + b*c - b*d)
U6 = [(1,0),(-1,0),(0,1),(0,-1),(1,1),(-1,-1)]
def primary(x):
    for u in U6:
        y = mul(u, x)
        if y[0] % 3 == 1 and y[1] % 3 == 0: return y
def isprime(n):
    if n < 2: return False
    i = 2
    while i*i <= n:
        if n % i == 0: return False
        i += 1
    return True
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
PR = []
for p in range(7, P+1):
    if not isprime(p) or p % 3 != 1: continue
    ab = two_sq(p)
    for x in (ab, (ab[0]-ab[1], -ab[1])):
        pi = primary(x); a, b = pi; wr = (-a*pow(b, -1, p)) % p
        fs = factor(p-1); g = 2
        while any(pow(g, (p-1)//q, p) == 1 for q in fs): g += 1
        ind = np.zeros(p, np.int64); v = 1
        for k in range(p-1): ind[v] = k; v = v*g % p
        gv = pow(g, (p-1)//6, p)
        img = {1: 0, wr % p: 2, wr*wr % p: 4, p-1: 3, (-wr) % p: 5, (-wr*wr) % p: 1}
        PR.append((pi, p, wr, img[gv], ind))
def gauss(n, N, facs, k):
    # chi_n^{6/k}: k=6 sextic, 3 cubic (chi^2), 2 quadratic (chi^3)
    e = 6 // k
    d = np.arange(N); ch = np.ones(N, complex)
    for (pi, p, wr, j0, ind) in facs:
        r = d % p
        ch *= np.where(r == 0, 0, W6[(e*j0*(ind[r] % 6)) % 6])
    T = (2*n[0] - n[1]) % N
    return np.sum(ch*np.exp(2j*np.pi*(d*T % N)/N))/math.sqrt(N)
rows = []
for i in range(len(PR)):
    rows.append(((PR[i][0], PR[i][1], [PR[i]])))
for i, j in itertools.combinations(range(len(PR)), 2):
    if PR[i][1] == PR[j][1]: continue
    N = PR[i][1]*PR[j][1]
    if N > 40000: continue
    rows.append((mul(PR[i][0], PR[j][0]), N, [PR[i], PR[j]]))
for i, j, k in itertools.combinations(range(len(PR)), 3):
    if len({PR[i][1], PR[j][1], PR[k][1]}) < 3: continue
    N = PR[i][1]*PR[j][1]*PR[k][1]
    if N > 60000: continue
    rows.append((mul(mul(PR[i][0], PR[j][0]), PR[k][0]), N, [PR[i], PR[j], PR[k]]))
print("elements", len(rows), file=sys.stderr)
def ph(z, m):  # phase index in multiples of 2pi/m
    return int(round((np.angle(z) / (2*np.pi)) * m)) % m
for mod in (4, 12, 36):
    tab2 = collections.defaultdict(set); tabe = collections.defaultdict(set)
    for (n, N, facs) in rows:
        g2 = gauss(n, N, facs, 2); g3 = gauss(n, N, facs, 3); g6 = gauss(n, N, facs, 6)
        mu = (-1)**len(facs)
        lam = complex(n[0] - n[1]/2, n[1]*math.sqrt(3)/2); lam /= abs(lam)
        eps = mu*g6/(lam*np.conj(g3))
        key = (n[0] % mod, n[1] % mod)
        tab2[key].add(ph(g2, 12)); tabe[key].add(ph(eps, 12))
    amb2 = sum(1 for v in tab2.values() if len(v) > 1); ambe = sum(1 for v in tabe.values() if len(v) > 1)
    print(f"mod {mod}: classes {len(tab2)}, g2 ambiguous {amb2}, eps ambiguous {ambe}")
