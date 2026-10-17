# composite check: eps(n) = mu(n) g6(n) / ((n/|n|) conj(g3(n))) for squarefree primary n (split primes only, for speed)
import math, cmath, numpy as np, sys
from collections import Counter
from eis import prime_ideals, mul, conj, norm, Z6, chi_prime as chp
X = int(float(sys.argv[1]))
LAM = (1, 2)
def tr_inv(n):
    d = mul(LAM, n); Nd = norm(d); e0, e1 = conj(d)
    return (2*e0 - e1)/Nd, (-2*e1 - (e0 - e1))/Nd
PR = [p for p in prime_ideals(X) if p['kind'] == 'split']
W = cmath.exp(2j*math.pi/3)
g6p, g3p = [], []
for pr in PR:
    g = pr['g']; P = pr['P']; r = pr['r']; c0, c1 = tr_inv(g)
    x = np.arange(1, P, dtype=np.int64)
    def vpow(b, ex, mod):
        res = np.ones_like(b); b = b % mod
        while ex:
            if ex & 1: res = (res*b) % mod
            b = (b*b) % mod; ex >>= 1
        return res
    h = vpow(x.copy(), (P-1)//6, P); zk = [(Z6[k][0] + Z6[k][1]*r) % P for k in range(6)]
    e = np.full(len(x), -1, np.int64)
    for k, v in enumerate(zk): e[h == v] = k
    phase = np.exp(2j*math.pi*((x*c0) % 1.0))
    g6p.append(np.sum(np.exp(1j*math.pi*e/3)*phase)/math.sqrt(P))
    g3p.append(np.sum(np.exp(2j*math.pi*e/3)*phase)/math.sqrt(P))
# products of two distinct split primes with norm <= X: twisted multiplicativity
# g6(pq) = chi_p(q) chi_q(p) g6(p) g6(q);  g3(pq) = chi_p(q)^2 chi_q(p)^2 g3(p) g3(q)   (chi = sextic symbol, exponent k -> zeta6^k)
cnt = Counter(); tot = 0
idx = list(range(len(PR)))
for i in idx:
    for j in idx:
        if j <= i: continue
        if PR[i]['P']*PR[j]['P'] > X: break
        e_pq = chp(PR[i], PR[j]['g']); e_qp = chp(PR[j], PR[i]['g'])
        cr = cmath.exp(1j*math.pi*((e_pq + e_qp) % 6)/3)
        g6 = cr*g6p[i]*g6p[j]
        g3 = cr**2*g3p[i]*g3p[j]
        nn = mul(PR[i]['g'], PR[j]['g']); z = complex(nn[0] + nn[1]*W); u = z/abs(z)
        eps = (+1)*g6/(u*np.conj(g3))        # mu(pq) = +1
        k = cmath.phase(eps)/(math.pi/6)
        cnt[round(k, 3)] += 1; tot += 1
print("pairs pq:", tot, " eps phases (units of pi/6):", dict(cnt))
cnt1 = Counter()
for i in idx:
    pr = PR[i]; z = complex(pr['g'][0] + pr['g'][1]*W); u = z/abs(z)
    eps = (-1)*g6p[i]/(u*np.conj(g3p[i])); cnt1[round(cmath.phase(eps)/(math.pi/6), 3)] += 1
print("primes p:", len(idx), " eps phases (units of pi/6):", dict(cnt1))
