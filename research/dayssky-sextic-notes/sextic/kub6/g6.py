# Sextic Gauss sums over O = Z[w], w = e^{2 pi i/3}.  Elements are pairs (a, b) = a + b w.
# Residue symbol (x/pi)_6 for split primes pi (N = p = 1 mod 6): x^{(p-1)/6} mod pi, read through w -> r mod p.
# Additive character on O/a: e(Tr(h x / (a sqrt(-3)))), well defined because sqrt(-3) is the different.
# Gauss sum g(h, a) = sum_{x mod a} chi_a(x) e(Tr(h x/(a sqrt(-3)))), chi_a = prod over primes of (./pi)_6.
import math, cmath, numpy as np
W = complex(-0.5, math.sqrt(3)/2)
def mul(x, y): a, b = x; c, d = y; return (a*c - b*d, a*d + b*c - b*d)
def nrm(x): a, b = x; return a*a - a*b + b*b
def cj(x): a, b = x; return (a - b, -b)
def tocx(x): return x[0] + x[1]*W
def primary(x):
    # unit multiple that is = 2 mod 3 in the sense a = 2 mod 3, b = 0 mod 3  (pi = -1 mod 3 convention)
    u = (1, 0)
    for k in range(6):
        a, b = x
        if a % 3 == 2 and b % 3 == 0: return x
        x = mul(x, (1, 1))  # multiply by -w^2 = 1 + w (a unit of order 6)
    raise ValueError
def rat_primes(n):
    s = np.ones(n+1, bool); s[:2] = False
    for i in range(2, int(n**0.5)+1):
        if s[i]: s[i*i::i] = False
    return np.nonzero(s)[0]
def split_primes(P):
    out = []
    for p in rat_primes(P):
        p = int(p)
        if p % 3 != 1: continue
        m = int((4*p/3)**0.5) + 2; found = []
        for b in range(-m, m+1):
            disc = b*b - 4*(b*b - p)
            if disc < 0: continue
            s = int(round(disc**0.5))
            if s*s != disc or (b + s) % 2: continue
            for a in ((b+s)//2, (b-s)//2):
                if nrm((a, b)) == p:
                    g = primary((a, b))
                    if g not in found: found.append(g)
        for g in found:
            a, b = g; r = (-a*pow(b, -1, p)) % p   # w -> r mod pi  (a + b r = 0 mod p)
            assert (r*r + r + 1) % p == 0
            out.append(dict(g=g, p=p, r=r))
    out.sort(key=lambda d: (d['p'], d['g'])); return out
ZETA6 = [cmath.exp(2j*math.pi*k/6) for k in range(6)]
class Prime:
    def __init__(self, d):
        self.g, self.p, self.r = d['g'], d['p'], d['r']
        p, r = self.p, self.r
        # discrete log table for (x/pi)_6 : x^{(p-1)/6} is a 6th root of unity in F_p; identify with zeta6^k via w -> r
        # zeta6 = -w^2 ; in F_p: w = r, zeta6 = -r^2
        z = (-r*r) % p; self.k = {}
        acc = 1
        for k in range(6): self.k[acc] = k; acc = acc*z % p
        self.e = (p-1)//6
    def red(self, x):  # O -> F_p
        return (x[0] + x[1]*self.r) % self.p
    def sym(self, x):  # (x/pi)_6 as exponent k (zeta6^k), or None if pi | x
        v = self.red(x)
        if v == 0: return None
        return self.k[pow(v, self.e, self.p)]
def addchar(h, x, a):
    # e(Tr(h x / (a sqrt(-3))))
    z = tocx(mul(h, x)) / (tocx(a) * complex(0, math.sqrt(3)))
    return cmath.exp(2j*math.pi*2*z.real)
def gauss_brute(h, a, primes_of_a):
    # sum over x mod a (representatives: x = s + t w, 0 <= s, t < N(a) suffices only for prime a with N(a) = p; general: use box mod a)
    N = nrm(a); tot = 0
    seen = set()
    # enumerate residues mod a by x = s + t w with 0 <= s,t < N, reduced via canonical key
    for s in range(N):
        for t in range(N):
            x = (s, t)
            key = tuple(P.red(x) for P in primes_of_a)
            if key in seen: continue
            seen.add(key)
            k = 0; zero = False
            for P in primes_of_a:
                kk = P.sym(x)
                if kk is None: zero = True; break
                k += kk
            if zero: continue
            tot += ZETA6[k % 6]*addchar(h, x, a)
        if len(seen) == N: break
    return tot
