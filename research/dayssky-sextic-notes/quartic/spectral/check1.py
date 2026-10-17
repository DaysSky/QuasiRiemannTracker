# (1) validate Z/N model against gauss4 (chi, psi, Kloosterman), (2) Weil bound |S(r,r';pi)|<=2 sqrt(q),
# (3) identity S(x^3,x^3;c) = chi_c(x) S(1,x^6;c), (4) Mellin identity S(1,n;pi)=(1/(q-1)) sum_lambda conj(lambda)(n) g(lambda) g(lambda chi)
import numpy as np, math, cmath, sys, os, random
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from zmod import ZmodSplit, FieldInert, split_primes
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from gauss4 import chi as chi4, gmul, gmod, gnorm, gconj, residues
random.seed(1)
# (1) validation on c = pi (N=13) and c = pi1*pi2
for c in [split_primes(13)[1], gmul(split_primes(13)[1], split_primes(17)[2])]:
    Z = ZmodSplit(c); N = Z.N; a, b = c
    bad = 0
    for x in range(N):
        z = (x, 0)
        if abs(Z.chi[x] - chi4(z, c)) > 1e-9: bad += 1
        zc = gmul(z, gconj(c)); ps = cmath.exp(2j*math.pi*zc[1]/N)
        if abs(Z.e(np.array([x]))[0] - ps) > 1e-9: bad += 1
    # Kloosterman directly in Z[i]
    r, rp = (2, 1), (3, -1)
    S = 0
    for d in residues(c):
        cd = chi4(d, c)
        if cd == 0: continue
        # inverse of d mod c
        for e in residues(c):
            if gmod((gmul(d, e)[0]-1, gmul(d, e)[1]), c) == (0, 0): break
        z = gmul(r, d); w = gmul(rp, e); s = (z[0]+w[0], z[1]+w[1]); sc = gmul(s, gconj(c))
        S += cd*cmath.exp(2j*math.pi*sc[1]/N)
    ri = (r[0]+r[1]*Z.iota) % N; rpi = (rp[0]+rp[1]*Z.iota) % N
    print("validate c=", c, "N=", N, "mismatches:", bad, " Kloost direct vs model:", abs(S-Z.kloost(ri, rpi)))
# (2)+(3)
worst = 0; worst_id = 0
for pi in split_primes(400, 13):
    Z = ZmodSplit(pi); N = Z.N
    for _ in range(20):
        r = random.randrange(1, N); rp = random.randrange(1, N)
        worst = max(worst, abs(Z.kloost(r, rp))/(2*math.sqrt(N)))
    for x in random.sample(range(1, N), 5):
        lhs = Z.kloost(pow(x, 3, N), pow(x, 3, N)); rhs = Z.chi[x]*Z.kloost(1, pow(x, 6, N))
        worst_id = max(worst_id, abs(lhs-rhs))
print("split primes N<=400: max |S|/(2 sqrt q) =", round(worst, 4), " max |S(x^3,x^3)-chi(x)S(1,x^6)| =", worst_id)
for p in [7, 11, 19, 23]:
    F = FieldInert(p); q = F.N; w = 0
    for _ in range(20):
        r = random.randrange(1, q); rp = random.randrange(1, q); w = max(w, abs(F.kloost(r, rp))/(2*math.sqrt(q)))
    print(f"inert p={p} q={q}: max |S|/(2 sqrt q) = {w:.4f}")
# (4) Mellin identity on one prime (q=29)
pi = split_primes(29, 29)[0]; Z = ZmodSplit(pi); N = Z.N
g = 2
while any(pow(g, (N-1)//t, N) == 1 for t in [2, 7]): g += 1   # generator mod 29 (29-1=28=4*7)
dlog = {pow(g, k, N): k for k in range(N-1)}
def lam(k, x): return cmath.exp(2j*math.pi*k*dlog[x]/(N-1))
def gs(k, extra_chi=False):
    s = 0
    for x in range(1, N):
        s += lam(k, x)*(Z.chi[x] if extra_chi else 1)*Z.e(np.array([x]))[0]
    return s
G = [gs(k)*gs(k, True) for k in range(N-1)]
err = 0
for n in [1, 2, 5, 17]:
    rhs = sum(np.conj(lam(k, n))*G[k] for k in range(N-1))/(N-1)
    err = max(err, abs(rhs - Z.kloost(1, n)))
print("Mellin identity S(1,n)=(1/(q-1)) sum_l conj(l(n)) g(l)g(l chi): max err", err)
