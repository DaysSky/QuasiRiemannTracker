# Fast sextic Kubota sums K(h) = sum_{a sqfree, split-prime support, N(a) in [Y/2,Y]} W(N(a)/Y) gt(h,a),  gt = g/sqrt(N).
# gt(h, a) = conj(chi_a(h)) * gt(1, a),  gt(1, a) = prod_i gt(1,pi_i) * prod_{i<j} (pi_j/pi_i)(pi_i/pi_j).
import sys, math, cmath, numpy as np
sys.path.insert(0, '.')
from g6 import split_primes, Prime, mul, nrm, tocx, ZETA6
def g1_prime(P):
    p = P.p; z = 1/(tocx(P.g)*complex(0, math.sqrt(3))); c = round(2*z.real*p) % p
    assert abs(2*z.real*p - round(2*z.real*p)) < 1e-6
    s = np.arange(1, p); ks = np.array([P.k[pow(int(x), P.e, p)] for x in s])
    return np.sum(np.exp(2j*np.pi*ks/6) * np.exp(2j*np.pi*c*s/p)) / math.sqrt(p)
def build(Y):
    ps = [Prime(d) for d in split_primes(Y)]
    for P in ps: P.g1 = g1_prime(P)
    # squarefree products of split primes with norm <= Y (DFS)
    mods = []  # list of (norm, prime index list)
    def dfs(start, nm, idx):
        for i in range(start, len(ps)):
            n2 = nm*ps[i].p
            if n2 > Y: break
            mods.append((n2, idx + [i])); dfs(i+1, n2, idx + [i])
    dfs(0, 1, [])
    return ps, mods
def W(t):
    t = np.asarray(t, float); s = 4*t - 3; out = np.zeros_like(t); m = np.abs(s) < 1
    out[m] = np.exp(1 - 1/(1 - s[m]**2)); return out
def coeffs(ps, mods, Y):
    """for moduli with N in [Y/2, Y]: weight W * gt(1,a), and the prime index lists"""
    sel = [(n, idx) for (n, idx) in mods if Y/2 <= n <= Y]
    c = np.zeros(len(sel), complex); lists = []
    for j, (n, idx) in enumerate(sel):
        v = complex(1)
        for i in idx: v *= ps[i].g1
        k = 0
        for a in range(len(idx)):
            for b in range(a+1, len(idx)):
                Pa, Pb = ps[idx[a]], ps[idx[b]]
                k += Pa.sym(Pb.g) + Pb.sym(Pa.g)
        c[j] = W(n/Y)*v*ZETA6[k % 6]; lists.append(idx)
    return c, lists
