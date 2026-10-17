# Task 1: verification of the sextic residue symbol implementation (reference, exact arithmetic).
import random, itertools, cmath, math
import numpy as np
from eis import *

random.seed(1)
X = 4000
PR = prime_ideals(X)
print("prime ideals outside S with norm <=", X, ":", len(PR))
W = cmath.exp(2j*math.pi/3)
def cplx(x): return x[0] + x[1]*W
LAM = cplx((1, 2))

def rand_elt(B=60):
    while True:
        u = (random.randint(-B, B), random.randint(-B, B))
        if u != (0, 0): return u

def rand_primary_odd(maxN):
    while True:
        a = 3*random.randint(-40, 40) + 1; b = 3*random.randint(-40, 40)
        n = (a, b)
        if is_odd(n) and 1 < norm(n) <= maxN:
            return n

fails = {}
def check(name, ok):
    fails.setdefault(name, [0, 0]); fails[name][0] += 1
    if not ok: fails[name][1] += 1

# 1. multiplicativity in u and in n, chi_n(a^6) = 1 on coprime a, supplementary laws
for _ in range(3000):
    pr = random.choice(PR)
    u1, u2 = rand_elt(), rand_elt()
    c1, c2, c12 = chi_prime(pr, u1), chi_prime(pr, u2), chi_prime(pr, mul(u1, u2))
    if ZERO in (c1, c2): check("mult_u(prime)", c12 == ZERO)
    else: check("mult_u(prime)", c12 == (c1 + c2) % 6)
    # Euler criterion cross-check by brute force for small norm: chi(u) == u^{(P-1)/6} built into chi_prime
    a = rand_elt(15); a6 = a
    for _k in range(5): a6 = mul(a6, a)
    ca = chi_prime(pr, a)
    check("chi_p(a^6)=1 (coprime)", chi_prime(pr, a6) == (ZERO if ca == ZERO else 0))
    P = pr['P']
    check("chi_p(w)=w^{(P-1)/6}", chi_prime(pr, (0, 1)) == (2*((P-1)//6)) % 6)
    check("chi_p(-1)=(-1)^{(P-1)/6}", chi_prime(pr, (-1, 0)) == (3*((P-1)//6)) % 6)
    # Euler criterion sanity: chi_p(u)^6 == 1 and chi_p(u) == 1 for u a 6th power residue
for _ in range(400):
    n1, n2 = rand_primary_odd(300), rand_primary_odd(300)
    u = rand_elt()
    c1, c2, c12 = chi(n1, u, PR), chi(n2, u, PR), chi(mul(n1, n2), u, PR)
    if ZERO in (c1, c2): check("mult_n", c12 == ZERO)
    else: check("mult_n", c12 == (c1 + c2) % 6)

# 2. sextic reciprocity chi_b(a) = r(a,b) chi_a(b) for coprime primary odd a,b (paper Lemma fixed-gauss-phase)
Rvals = {}
for _ in range(1500):
    a, b = rand_primary_odd(600), rand_primary_odd(600)
    cab, cba = chi(b, a, PR), chi(a, b, PR)   # chi_b(a), chi_a(b)
    if cab == ZERO:
        check("coprime<->nonzero", cba == ZERO); continue
    r = Rfac(a, b)
    check("R(a,b) is +-1", abs(r - round(r.real)) < 1e-12 and abs(abs(r) - 1) < 1e-12)
    lhs = zeta6(cab); rhs = r * zeta6(cba)
    check("sextic reciprocity chi_b(a)=R(a,b)chi_a(b)", abs(lhs - rhs) < 1e-9)
    check("cubic reciprocity (a/b)_3=(b/a)_3", (2*cab) % 6 == (2*cba) % 6)
    Rvals[round(r.real)] = Rvals.get(round(r.real), 0) + 1
# 3. diagonal R(n,n) = chi_n(-1)
for _ in range(500):
    n = rand_primary_odd(2000)
    check("R(n,n)=chi_n(-1)", abs(Rfac(n, n) - zeta6(chi(n, (-1, 0), PR))) < 1e-9)

# 4. Gauss-sum orientation: gamma_2(p)^3 = -alpha(p) and gamma_1 gamma_2 = mu alpha G (paper eq. (signal))
def e_add(z):   # e(z) = exp(2 pi i Tr(z/lambda))
    return cmath.exp(2j*math.pi*2*(z/LAM).real)
def gamma(j, c, residues):
    cc = cplx(c); s = 0
    for v in residues:
        e = chi(c, v, PR)
        if e == ZERO: continue
        s += cmath.exp(1j*math.pi*((j*e) % 6)/3) * e_add(cplx(v)/cc)
    return s / math.sqrt(norm(c))
split = [p for p in PR if p['kind'] == 'split' and p['P'] < 400]
for pr in split:
    c = pr['g']; P = pr['P']
    res = [(v, 0) for v in range(P)]
    g1, g2, g3 = gamma(1, c, res), gamma(2, c, res), gamma(3, c, res)
    alpha = cplx(c)/abs(cplx(c))
    check("gamma_2(p)^3 = -alpha(p)", abs(g2**3 + alpha) < 1e-8)
    G = np.conj(zeta6(chi(c, (4, 0), PR))) * Gamma(c)
    check("gamma_1 gamma_2 = mu alpha G (prime)", abs(g1*g2 - (-1)*alpha*G) < 1e-8)
    check("Gamma(p)=gamma_3(p)", abs(Gamma(c) - g3) < 1e-8)
# composite c = p1 p2 with distinct split rational norms (O/c = Z/N(c))
pairs = [(p, q) for p, q in itertools.combinations(split[:14], 2) if p['P'] != q['P']][:25]
for p, q in pairs:
    c = mul(p['g'], q['g']); N = norm(c)
    res = [(v, 0) for v in range(N)]
    g1, g2 = gamma(1, c, res), gamma(2, c, res)
    alpha = cplx(c)/abs(cplx(c))
    G = np.conj(zeta6(chi(c, (4, 0), PR))) * Gamma(c)
    check("gamma_2(c)^3 = mu(c) alpha(c) (c=p1p2)", abs(g2**3 - alpha) < 1e-7)
    check("gamma_1 gamma_2 = mu alpha G (c=p1p2)", abs(g1*g2 - alpha*G) < 1e-7)

for k, (n, f) in fails.items():
    print(f"{k:48s} tests={n:5d} failures={f}")
print("R(a,b) value counts:", Rvals)
