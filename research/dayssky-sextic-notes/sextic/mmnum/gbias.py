# Sextic Gauss sums over Z[w]: g(n) = sum_{x mod n} chi_n(x) e(Tr(x/(lam n))), lam = sqrt(-3) = 1+2w.
# Test for Patterson-type bias: S(X) = sum_{N(n)<=X, n sf primary odd} c(n) gtilde(n), c = 1 or mu, and over primes.
import sys, math, cmath, numpy as np
from eis import prime_ideals, mul, conj, norm, Z6
from mmdata import Ideals
X = int(float(sys.argv[1]))
W = cmath.exp(2j*math.pi/3)
def cx(x): return x[0] + x[1]*W
LAM = (1, 2)
def tr_inv(n):
    # Tr(y/(lam n)) for y in Z[w] is Z-linear in (y0,y1): returns coefficients (c0,c1) rational
    d = mul(LAM, n); Nd = norm(d); cd = conj(d)          # 1/d = conj(d)/N(d)
    # y*conj(d) = (y0 + y1 w)(e0 + e1 w); Tr(a + b w) = 2a - b
    e0, e1 = cd
    # (y0 + y1 w)(e0 + e1 w) = y0e0 - y1e1 + (y0e1 + y1e0 - y1e1) w
    c0 = (2*e0 - e1) / Nd           # coefficient of y0 in Tr
    c1 = (-2*e1 - (e0 - e1)) / Nd   # coefficient of y1: 2(-e1) - (e0 - e1)
    return c0, c1
PR = prime_ideals(X)
I = Ideals(X)
# Gauss sums at primes
gp = np.zeros(len(PR), complex)
for j, pr in enumerate(PR):
    g = pr['g']; P = pr['P']
    c0, c1 = tr_inv(g)
    if pr['kind'] == 'split':
        r = pr['r']           # w -> r mod P ; residues represented by integers x in 0..P-1 (y1 = 0)
        x = np.arange(1, P, dtype=np.int64)
        def vpow(b, ex, mod):
            res = np.ones_like(b); b = b % mod
            while ex:
                if ex & 1: res = (res * b) % mod
                b = (b * b) % mod; ex >>= 1
            return res
        h = vpow(x.copy(), (P-1)//6, P)
        zk = [(Z6[k][0] + Z6[k][1]*r) % P for k in range(6)]
        e = np.full(len(x), -1, np.int64)
        for k, v in enumerate(zk): e[h == v] = k
        assert (e >= 0).all()
        phase = np.exp(2j*math.pi*((x*c0) % 1.0))
        gp[j] = np.sum(np.exp(1j*math.pi*e/3) * phase)
    else:
        q = pr['q']   # norm q^2; residues x0 + x1 w, 0<=x0,x1<q
        from eis import chi_prime
        s = 0
        for x0 in range(q):
            for x1 in range(q):
                if x0 == 0 and x1 == 0: continue
                k = chi_prime(pr, (x0, x1))
                s += cmath.exp(1j*math.pi*k/3) * cmath.exp(2j*math.pi*((x0*c0 + x1*c1) % 1.0))
        gp[j] = s
gt_p = gp / np.sqrt(np.array([p['P'] for p in PR], float))
print("primes", len(PR), " |gtilde| mean", np.abs(gt_p).mean(), flush=True)
# squarefree n via twisted multiplicativity: g(n p) = chi_n(p) chi_p(n) g(n) g(p) for coprime primary n, p
from eis import chi_prime as chp
n_ = len(I.N); gt = np.zeros(n_, complex); gt[0] = 1.0
for i in range(1, n_):
    if I.mu[i] == 0: continue
    j = int(I.spf[i]); c = int(I.cof[i])
    if c == 0: gt[i] = gt_p[j]; continue
    pr = PR[j]; ncof = I.g[c]
    # chi_n(p): symbol of p modulo n = prod over primes of n; chi_p(n) = chp(pr, n)
    e_pn = chp(pr, ncof)               # chi_p(n)
    # chi_n(p) = prod_{q | n} chi_q(p)
    e_np = 0; k = c
    while k > 0:
        jj = int(I.spf[k]); e_np = (e_np + chp(PR[jj], pr['g'])) % 6; k = int(I.cof[k])
    gt[i] = cmath.exp(1j*math.pi*((e_pn + e_np) % 6)/3) * gt[c] * gt_p[j]
sq = I.mu != 0
for lab, c in (("all sf", sq.astype(float)), ("mu", I.mu.astype(float))):
    cs = np.cumsum(c * gt)
    for Xt in (X/16, X/4, X):
        k = int(np.searchsorted(I.N, Xt, 'right')) - 1
        S = cs[k]
        print(f"{lab:7s} X={Xt:9.0f} S={S.real:+10.2f}{S.imag:+10.2f}i |S|={abs(S):9.2f} |S|/X^(2/3)={abs(S)/Xt**(2/3):.4f} |S|/X^(1/2)={abs(S)/Xt**0.5:.4f}", flush=True)
Pn = np.array([p['P'] for p in PR], float)
for Xt in (X/16, X/4, X):
    m = Pn <= Xt
    S = np.sum(gt_p[m] * np.log(Pn[m]))
    print(f"primes  X={Xt:9.0f} |S_Lambda|={abs(S):9.2f} /X^(2/3)={abs(S)/Xt**(2/3):.4f} /X^(1/2)={abs(S)/Xt**0.5:.4f} /(X^(1/2)logX)={abs(S)/(Xt**0.5*math.log(Xt)):.4f}")
