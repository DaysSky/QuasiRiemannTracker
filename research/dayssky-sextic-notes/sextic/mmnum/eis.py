# Eisenstein integers O = Z[w], w = exp(2 pi i/3); element (a,b) <-> a + b w.
# Conventions follow OpenAI paper.tex Sec. "Arithmetic and analytic preliminaries":
#   primary = congruent to 1 mod 3 (a = 1 mod 3, b = 0 mod 3); S = primes above 2,3;
#   chi_p(u) = sixth root of unity == u^{(Np-1)/6} mod p ; chi_c = prod chi_p^{v_p(c)}.
# Sixth roots of unity are encoded by an exponent e in 0..5 meaning zeta6^e, zeta6 = 1 + w = exp(i pi/3);
# ZERO = 6 encodes the value 0.
import numpy as np

ZERO = 6

def mul(x, y):
    a, b = x; c, d = y
    return (a*c - b*d, a*d + b*c - b*d)

def conj(x):
    a, b = x            # conj(w) = w^2 = -1 - w
    return (a - b, -b)

def norm(x):
    a, b = x
    return a*a - a*b + b*b

def divexact(x, y):
    """x / y assuming y | x."""
    n = norm(y); p = mul(x, conj(y))
    assert p[0] % n == 0 and p[1] % n == 0, (x, y)
    return (p[0]//n, p[1]//n)

def divides(y, x):
    n = norm(y); p = mul(x, conj(y))
    return p[0] % n == 0 and p[1] % n == 0

UNITS = [(1,0),(0,1),(-1,-1),(-1,0),(0,-1),(1,1)]   # zeta6^k? no: list of all six units
# zeta6 = 1 + w ; powers zeta6^k for k=0..5
Z6 = [(1,0)]
for _k in range(5):
    Z6.append(mul(Z6[-1], (1,1)))
assert mul(Z6[-1], (1,1)) == (1,0)

def is_primary(x):
    return x[0] % 3 == 1 and x[1] % 3 == 0

def primary_assoc(x):
    """unique associate congruent to 1 mod 3 (x coprime to 3)."""
    for u in Z6:
        y = mul(u, x)
        if is_primary(y):
            return y
    raise ValueError("not coprime to 3: %r" % (x,))

def is_odd(x):
    return not (x[0] % 2 == 0 and x[1] % 2 == 0)

def rat_primes(n):
    s = np.ones(n+1, bool); s[:2] = False
    for i in range(2, int(n**0.5)+1):
        if s[i]: s[i*i::i] = False
    return np.nonzero(s)[0]

def prime_ideals(X):
    """primary generators of prime ideals outside S={2,3}-primes with norm <= X, sorted by norm.
       returns list of dicts: g=(a,b) primary generator, P=norm, kind='split'/'inert', r (w -> r mod P) or q."""
    out = []
    for p in rat_primes(X):
        p = int(p)
        if p in (2, 3): continue
        if p % 3 == 1:
            # find a,b with a^2-ab+b^2 = p
            found = []
            m = int((4*p/3)**0.5) + 2
            for b in range(-m, m+1):
                # a^2 - a b + b^2 - p = 0
                disc = b*b - 4*(b*b - p)
                if disc < 0: continue
                s = int(round(disc**0.5))
                if s*s != disc: continue
                for a in ((b + s)//2, (b - s)//2):
                    if (b + s) % 2 == 0 and norm((a, b)) == p:
                        g = primary_assoc((a, b))
                        if g not in found: found.append(g)
            assert len(found) == 2, (p, found)
            for g in sorted(found):
                a, b = g
                r = (-a * pow(b, -1, p)) % p
                assert (r*r + r + 1) % p == 0
                out.append(dict(g=g, P=p, kind='split', r=r))
        else:
            if p*p <= X:
                out.append(dict(g=(-p, 0), P=p*p, kind='inert', q=p))
    out.sort(key=lambda d: (d['P'], d['g']))
    return out

# ---------- reference sextic symbol (slow, exact) ----------
def _fq2_mul(x, y, q):
    a, b = x; c, d = y
    return ((a*c - b*d) % q, (a*d + b*c - b*d) % q)

def _fq2_pow(x, e, q):
    r = (1, 0); x = (x[0] % q, x[1] % q)
    while e:
        if e & 1: r = _fq2_mul(r, x, q)
        x = _fq2_mul(x, x, q); e >>= 1
    return r

def chi_prime(pr, u):
    """sextic symbol (u/p)_6 as exponent in 0..5, or ZERO."""
    x, y = u
    if pr['kind'] == 'split':
        P, r = pr['P'], pr['r']
        t = (x + y*r) % P
        if t == 0: return ZERO
        h = pow(t, (P-1)//6, P)
        for k in range(6):
            zk = (Z6[k][0] + Z6[k][1]*r) % P
            if zk == h: return k
        raise RuntimeError
    else:
        q = pr['q']
        if x % q == 0 and y % q == 0: return ZERO
        h = _fq2_pow((x, y), (q*q-1)//6, q)
        for k in range(6):
            if (Z6[k][0] % q, Z6[k][1] % q) == h: return k
        raise RuntimeError

def factor_primary(n, primes):
    """factor a primary odd element n over prime list (trial division)."""
    res = []
    m = n
    for pr in primes:
        if norm(m) == 1: break
        if norm(m) % pr['P'] != 0 and pr['kind'] == 'split': continue
        v = 0
        while divides(pr['g'], m):
            m = divexact(m, pr['g']); v += 1
        if v: res.append((pr, v))
    assert norm(m) == 1, (n, m)
    return res

def chi(n, u, primes):
    """chi_n(u) for primary odd n (exponent or ZERO)."""
    e = 0
    for pr, v in factor_primary(n, primes):
        c = chi_prime(pr, u)
        if c == ZERO: return ZERO
        e = (e + v*c) % 6
    return e

# ---------- quadratic four-term phase and reciprocity factor (paper Lemma quadratic-four-term) ----------
def Gamma(c):
    a, b = c
    ipow = lambda k: [1, 1j, -1, -1j][k % 4]
    return (1 + ipow(-b) + ipow(a) + ipow(b - a)) / 2

def Rfac(a, b):
    return Gamma(mul(a, b)) / (Gamma(a) * Gamma(b))

def zeta6(e):
    return 0 if e == ZERO else np.exp(1j*np.pi*e/3)
