# Möbius-weighted cubic Gauss sums over Z[omega]: is there a Patterson X^{5/6} bias?
# Elements a+b*w, w = e^{2 pi i/3}, w^2 = -1-w. Norm a^2 - ab + b^2.
# g(c) = sum_{d mod c} (d/c)_3 e(Tr(d/(c*sqrt(-3)))) -- we use Dunn-Radziwill convention e~(z) = exp(2 pi i (z + zbar)) with
# argument d/c (as written in their (2.1) up to the lambda=sqrt(-3) normalisation; the bias exponent does not depend on it).
import sys, math, cmath, numpy as np
X = int(sys.argv[1])
def mul(x, y):
    a, b = x; c, d = y
    return (a*c - b*d, a*d + b*c - b*d)
def norm(x): a, b = x; return a*a - a*b + b*b
def conj(x): a, b = x; return (a - b, -b)
def prim(x):  # x == 1 mod 3
    a, b = x; return a % 3 == 1 and b % 3 == 0
# enumerate primes: split p = 1 mod 3 -> two conjugate primary primes of norm p; inert q = 2 mod 3 -> primary -q, norm q^2
def sieve(n):
    s = np.ones(n+1, bool); s[:2] = False
    for i in range(2, int(n**.5)+1):
        if s[i]: s[i*i::i] = False
    return np.nonzero(s)[0]
ps = sieve(X)
units = [(1,0),(-1,0),(0,1),(0,-1),(1,1),(-1,-1)]
def make_primary(x):
    for u in units:
        y = mul(u, x)
        if prim(y): return y
    return None
primes = []  # (pi, norm, rational p)
for p in ps:
    p = int(p)
    if p == 3: continue
    if p % 3 == 1:
        # find a,b with a^2 - ab + b^2 = p
        found = None
        for b in range(1, int(math.isqrt(4*p//3))+2):
            disc = 4*p - 3*b*b
            if disc < 0: break
            r = math.isqrt(disc)
            if r*r == disc:
                for a in [(b + r)//2, (b - r)//2]:
                    if (b + r) % 2 == 0 and a*a - a*b + b*b == p: found = (a, b); break
            if found: break
        pi = make_primary(found); pib = make_primary(conj(found))
        primes.append((pi, p, p)); primes.append((pib, p, p))
    elif p*p <= X:
        primes.append((make_primary((p,0)), p*p, p))
print("primes:", len(primes), file=sys.stderr)
# residue map for split primes: Z[w]/pi ~ F_p, w -> wr with pi | (w - wr)
def wroot(pi, p):
    a, b = pi
    # pi = a + b w = 0 mod pi -> w = -a/b mod p (b invertible mod p since p splits and pi not rational)
    return (-a * pow(b, -1, p)) % p
def cubic_char_table(p):
    # discrete log mod p via primitive root
    g = 2
    while True:
        ok = all(pow(g, (p-1)//q, p) != 1 for q in set(factor(p-1)))
        if ok: break
        g += 1
    ind = np.zeros(p, np.int64); x = 1
    for k in range(p-1): ind[x] = k; x = x*g % p
    return g, ind
def factor(n):
    f = []; d = 2
    while d*d <= n:
        while n % d == 0: f.append(d); n //= d
        d += 1
    if n > 1: f.append(n)
    return f
W = cmath.exp(2j*math.pi/3)
info = {}
def symbol_split(x, pi, p, wr, g, ind):
    # (x/pi)_3 for x in Z[w]: reduce x = a + b w -> a + b wr mod p, then x^{(p-1)/3} = wr^k
    a, b = x; r = (a + b*wr) % p
    if r == 0: return 0
    e = pow(r, (p-1)//3, p)
    # identify e with w^k: w -> wr
    if e == 1: return 1
    if e == wr: return W
    if e == wr*wr % p: return W*W
    raise ValueError
gauss = {}
for (pi, nm, p) in primes:
    if nm != p: continue
    wr = wroot(pi, p); g, ind = cubic_char_table(p)
    info[pi] = (p, wr)
    # (d/pi)_3 for rational d = 0..p-1 ; additive char: e~(d/pi) = exp(2 pi i Tr(d/pi)) with Tr(d/pi) = d*Tr(pibar)/p
    d = np.arange(1, p)
    e = np.array([pow(int(x), (p-1)//3, p) for x in d]) if p < 2000 else None
    if e is None:
        k = ind[d] % 3  # chi(d) = (g^{ind})^{(p-1)/3} -> need value of g^{(p-1)/3} as w^j
        gv = pow(g, (p-1)//3, p); j = 1 if gv == wr else 2
        chi = W**((j*k) % 3)
    else:
        chi = np.where(e == 1, 1, np.where(e == wr, W, W*W))
    trpb = (2*pi[0] - pi[1])  # Tr(pibar) = Tr(pi) = 2a - b
    # Dunn-Radziwill: e~(z) = exp(2 pi i (z + zbar)); z = d/pi -> z + zbar = d*(pibar + pi)/p = d*Tr(pi)/p
    add = np.exp(2j*np.pi*(d*trpb % p)/p)
    gauss[pi] = complex(np.sum(chi*add))/math.sqrt(p)
print("gauss done", file=sys.stderr)
np.save('gauss_%d.npy' % X, np.array([(pi[0], pi[1], gauss[pi].real, gauss[pi].imag) for pi in gauss]))
