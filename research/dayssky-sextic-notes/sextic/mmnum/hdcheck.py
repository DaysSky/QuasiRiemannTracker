# Check Hasse-Davenport: g6(p) = chibar(4) chi^2(-1) g3(p)^2 g2(p)/Q, and whether mu(p) g~2(p) is a fixed character of p.
import math, cmath, numpy as np, sys
from eis import prime_ideals, mul, conj, norm, Z6
X = int(float(sys.argv[1]))
LAM = (1, 2)
def tr_inv(n):
    d = mul(LAM, n); Nd = norm(d); e0, e1 = conj(d)
    return (2*e0 - e1)/Nd, (-2*e1 - (e0 - e1))/Nd
PR = [p for p in prime_ideals(X) if p['kind'] == 'split']
rows = []
for pr in PR:
    g = pr['g']; P = pr['P']; r = pr['r']; c0, c1 = tr_inv(g)
    x = np.arange(1, P, dtype=np.int64)
    def vpow(b, ex, mod):
        res = np.ones_like(b); b = b % mod
        while ex:
            if ex & 1: res = (res*b) % mod
            b = (b*b) % mod; ex >>= 1
        return res
    h = vpow(x.copy(), (P-1)//6, P)
    zk = [(Z6[k][0] + Z6[k][1]*r) % P for k in range(6)]
    e = np.full(len(x), -1, np.int64)
    for k, v in enumerate(zk): e[h == v] = k
    phase = np.exp(2j*math.pi*((x*c0) % 1.0))
    g6 = np.sum(np.exp(1j*math.pi*e/3)*phase)/math.sqrt(P)
    g3 = np.sum(np.exp(2j*math.pi*e/3)*phase)/math.sqrt(P)
    g2 = np.sum(np.exp(1j*math.pi*e)*phase)/math.sqrt(P)
    # chi_p(4) and chi_p(-1): exponents of 4 and -1 (=P-1) as residues
    e4 = int(e[4-1]) if P > 4 else 0; em1 = int(e[P-2])
    pred = cmath.exp(-1j*math.pi*e4/3) * cmath.exp(1j*math.pi*2*em1/3)
    ratio = g6/(g3**2*g2)
    rows.append((P, g[0], g[1], ratio, pred, g2, em1))
bad = sum(1 for t in rows if abs(t[3]-t[4]) > 1e-6)
print("split primes", len(rows), " Hasse-Davenport mismatches:", bad)
# quadratic Gauss sum: is -g2 (= mu(p) g2) determined by (a,b) mod small moduli?
from collections import defaultdict
for mod in (4, 6, 12, 36):
    tab = defaultdict(set)
    for P, a, b, ratio, pred, g2, em1 in rows:
        tab[(a % mod, b % mod)].add((round(g2.real, 6), round(g2.imag, 6)))
    amb = sum(1 for k, v in tab.items() if len(v) > 1)
    print(f"g2 determined by (a,b) mod {mod}? classes {len(tab)}, ambiguous {amb}; values {sorted(set(v for s in tab.values() for v in s))}")
# cubic Gauss sum cubed vs the Groessencharacter g/|g| (or its conjugate): ratio a fixed root of unity?
W3 = cmath.exp(2j*math.pi/3)
vals_g, vals_gb = set(), set()
for pr in PR[:600]:
    pass
