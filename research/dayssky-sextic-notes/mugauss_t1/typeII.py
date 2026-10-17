# Sums of normalised cubic Gauss sums over squarefree primary c in Z[w] (w = e^{2pi i/3}), Dunn-Radziwill convention
# g(c) = sum_{d mod c} (d/c)_3 exp(2 pi i Tr(d/c)),  |g~(c)| = mu^2(c),  g~(pi)^3 = -pi/|pi|.
# Outputs partial sums of  g~(c), mu(c) g~(c), and prime g~ at norms X/2^k.
import sys, math, cmath, numpy as np
X = int(sys.argv[1])
W = cmath.exp(2j*math.pi/3)
def mul(x, y):
    a, b = x; c, d = y
    return (a*c - b*d, a*d + b*c - b*d)
def nrm(x): a, b = x; return a*a - a*b + b*b
def cj(x): a, b = x; return (a - b, -b)
def prim(x): return x[0] % 3 == 1 and x[1] % 3 == 0
U6 = [(1,0),(-1,0),(0,1),(0,-1),(1,1),(-1,-1)]
def primary(x):
    for u in U6:
        y = mul(u, x)
        if prim(y): return y
def sieve(n):
    s = np.ones(n+1, bool); s[:2] = False
    for i in range(2, int(n**.5)+1):
        if s[i]: s[i*i::i] = False
    return np.nonzero(s)[0]
def factor(n):
    f = set(); d = 2
    while d*d <= n:
        while n % d == 0: f.add(d); n //= d
        d += 1
    if n > 1: f.add(n)
    return f
def primroot(p):
    fs = factor(p-1); g = 2
    while any(pow(g, (p-1)//q, p) == 1 for q in fs): g += 1
    return g
def two_squares(p):  # a^2 - ab + b^2 = p
    for b in range(1, int(math.isqrt(4*p//3))+2):
        disc = 4*p - 3*b*b
        if disc < 0: break
        r = math.isqrt(disc)
        if r*r == disc and (b + r) % 2 == 0: return ((b + r)//2, b)
# --- F_{q^2} arithmetic for inert q: elements (a,b) = a + b w mod q
def fq2_mul(x, y, q):
    a, b = x; c, d = y
    return ((a*c - b*d) % q, (a*d + b*c - b*d) % q)
def fq2_pow(x, e, q):
    r = (1, 0)
    while e:
        if e & 1: r = fq2_mul(r, x, q)
        x = fq2_mul(x, x, q); e >>= 1
    return r
def cubic_sym(c, P):
    """(c/P)_3 as a complex cube root of unity (0 if P | c). P = (pi, kind, p, extra)."""
    pi, kind, p, ex = P
    if kind == 's':
        wr = ex; r = (c[0] + c[1]*wr) % p
        if r == 0: return 0
        e = pow(r, (p-1)//3, p)
        return 1 if e == 1 else (W if e == wr else W*W)
    else:
        q = p; r = (c[0] % q, c[1] % q)
        if r == (0, 0): return 0
        e = fq2_pow(r, (q*q-1)//3, q)
        return {(1,0): 1, (0,1): W, ((-1) % q, (-1) % q): W*W}[e]
ps = sieve(X)
Pl = []
for p in ps:
    p = int(p)
    if p == 3: continue
    if p % 3 == 1:
        ab = two_squares(p)
        for x in (ab, cj(ab)):
            pi = primary(x); a, b = pi
            wr = (-a * pow(b, -1, p)) % p
            Pl.append((pi, 's', p, wr))
    elif p*p <= X:
        Pl.append((primary((p, 0)), 'i', p, None))
Pl.sort(key=lambda P: nrm(P[0]))
gt = {}
for P in Pl:
    pi, kind, p, wr = P
    if kind == 's':
        g = primroot(p); gv = pow(g, (p-1)//3, p); j = 1 if gv == wr else 2
        y = np.arange(1, p, dtype=np.int64); y3 = (y*y % p)*y % p
        eta = [np.exp(2j*np.pi*((pow(g, r, p)*y3) % p)/p).sum()/3 for r in range(3)]
        G = sum((W**((j*r) % 3))*eta[r] for r in range(3))
        T = (2*pi[0] - pi[1]) % p
        chiT = cubic_sym((T, 0), P)
        gt[pi] = (chiT.conjugate() if chiT != 0 else 0)*G/math.sqrt(p)
    else:
        q = p; s = 0
        for a in range(q):
            for b in range(q):
                if a == 0 and b == 0: continue
                chi = cubic_sym((a, b), P)
                tr = (2*a - b) * 1.0
                s += chi*cmath.exp(2j*math.pi*(-tr/q))   # Tr(d/(-q)) = -Tr(d)/q
        gt[pi] = s/q
print("primes", len(Pl), "check cube:", max(abs(gt[P[0]]**3 + complex(P[0][0] + P[0][1]*W.real, P[0][1]*W.imag)/math.sqrt(nrm(P[0]))) for P in Pl), file=sys.stderr, flush=True)

# Balanced Type II weight: w(c) = #{d | c : X^0.4 < N(d) <= X^0.6}. Compare sum mu g~ w (Type II for Mobius) with sum g~ w (theta-biased).
import itertools
lo, hi = X**0.4, X**0.6
K = 24
Smw = np.zeros(K, complex); Sw = np.zeros(K, complex)
def bucket(nm):
    return int(math.floor(math.log2(X/nm))) if nm <= X else -1
def wcount(norms):
    t = 0
    for r in range(1, len(norms)):
        for S in itertools.combinations(norms, r):
            pr = 1
            for v in S: pr *= v
            if lo < pr <= hi: t += 1
    return t
def rec(start, c, nm, gval, norms):
    for i in range(start, len(Pl)):
        P = Pl[i]; pi = P[0]; np_ = nrm(pi)
        if nm*np_ > X: break
        s = cubic_sym(c, P) if norms else 1
        g2 = gval*gt[pi]*(s*s)
        c2 = mul(c, pi); n2 = nm*np_; nl = norms + [np_]
        k = bucket(n2)
        if k < K and len(nl) >= 2:
            w = wcount(nl)
            if w:
                Sw[k] += w*g2; Smw[k] += w*((-1)**len(nl))*g2
        rec(i+1, c2, n2, g2, nl)
sys.setrecursionlimit(10000)
rec(0, (1, 0), 1, 1.0, [])
print("X", X, " balanced Type II weight d in (X^.4, X^.6]  (cumulative over N(c) <= Xk, with the d-range fixed by X)")
cw = np.cumsum(Sw[::-1])[::-1]; cmw = np.cumsum(Smw[::-1])[::-1]
for k in range(0, 4):
    Xk = X/2**k
    print(f"{Xk:9.0f}  sum g~ w / Xk^(5/6) = {cw[k].real/Xk**(5/6):+.4f}   |sum mu g~ w| / Xk^(1/2) = {abs(cmw[k])/Xk**0.5:.3f}   / Xk^(5/6) = {abs(cmw[k])/Xk**(5/6):.4f}")
