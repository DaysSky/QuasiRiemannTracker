# Task 3 numerics: principal rows u = a^6 and fixed-character rows u = k0 a^6 (masks (n,a)=1).
#   chi_n(k0 a^6) = chi_n(k0) 1_{(n,a)=1}.
import sys, math, json
import numpy as np
from mmdata import Ideals, W
from eis import chi_prime, ZERO, mul

X = int(float(sys.argv[1])) if len(sys.argv) > 1 else 100000
I = Ideals(X)
n = len(I.N)
# prime-divisibility masks for small primes (for coprimality to a)
A0 = 300
small = [j for j, p in enumerate(I.primes) if p['P'] <= A0]
divmask = {j: np.zeros(n, bool) for j in small}
for i in range(1, n):
    for j in set(I.factor(i)):
        if j in divmask: divmask[j][i] = True
avals = [i for i in range(n) if I.N[i] <= A0]          # ideals a (primary generators), N(a) <= A0
def cop(a_idx):
    m = np.ones(n, bool)
    for j in set(I.factor(a_idx)):
        m &= ~divmask[j]
    return m
COP = {a: cop(a) for a in avals}

def chivec(k0):
    """chi_n(k0) for all ideals n as complex array."""
    ep = np.array([chi_prime(p, k0) for p in I.primes])
    e = np.zeros(n, np.int64)
    for i in range(1, n):
        a, b = ep[I.spf[i]], e[I.cof[i]]
        e[i] = ZERO if (a == ZERO or b == ZERO) else (a + b) % 6
    z = np.exp(1j*np.pi*e/3); z[e == ZERO] = 0
    return z

def poly(coef, L, mask=None):
    lo, hi = I.rng(L/2, L)
    c = coef[lo:hi] * W(I.N[lo:hi]/L)
    if mask is not None: c = c * mask[lo:hi]
    return c.sum() / math.sqrt(L)

mu = I.mu.astype(float)
out = {}
Ds = [10**k for k in (3, 3.5, 4, 4.5, 5)] + ([3e5] if X >= 3e5 else [])
print(f"X={X}, #a (N(a)<= {A0}) = {len(avals)}")
print("principal rows a^6: M_{a^6}(D) = D^{-1/2} sum_{(n,a)=1} mu(n) W(Nn/D)")
for D in Ds:
    Ms = np.array([poly(mu, D, COP[a]) for a in avals])
    print(f"  D={D:9.0f}: mean_a |M|^2 = {np.mean(abs(Ms)**2):.4f}  max_a |M|^2 = {np.max(abs(Ms)**2):.4f}"
          f"  |mean_a M|^2 = {abs(Ms.mean())**2:.2e}   (trivial bound ~ D; zero at beta would give ~D^(2beta-1))")
one = np.ones(n)
print("principal rows: |S_{a^6}(N)|^2 / N  (expected ~ const * (phi(a)/N(a))^2 )")
for N in [1e2, 1e3, 1e4, 1e5]:
    Ss = np.array([poly(one, N, COP[a]) for a in avals])
    print(f"  N={N:8.0f}: mean_a |S|^2/N = {np.mean(abs(Ss)**2)/N:.4f}   a=1: {abs(Ss[0])**2/N:.4f}")
print("fixed nonprincipal rows u = k0 (and k0 a^6): |S_{k0}(N)|^2 and |M_{k0}(D)|^2")
K0 = {'w': (0, 1), '-1': (-1, 0), 'lambda': (1, 2), '2': (2, 0), 'pi7': I.primes[0]['g'],
      'pi7*pi13': mul(I.primes[0]['g'], I.primes[2]['g'])}
for name, k0 in K0.items():
    cv = chivec(k0)
    Sline = "  ".join(f"N=1e{math.log10(N):.1f}:{abs(poly(cv, N))**2:.1e}" for N in [1e2, 1e3, 1e4, 1e5])
    Mline = "  ".join(f"D=1e{math.log10(D):.1f}:{abs(poly(mu*cv, D))**2:.3f}" for D in [1e3, 1e4, 1e5])
    Smask = max(abs(poly(cv, 1e4, COP[a]))**2 for a in avals)
    print(f"  k0={name:9s} |S|^2: {Sline}   max_a|S_(k0 a^6)(1e4)|^2={Smask:.1e}   |M|^2: {Mline}")
# barrier table: all-row MM at r+m=1 restricted to principal rows a^6 forces zeta_F zero-free for Re s > theta(m)
print("barrier exponent theta(m) = (11-12m)/(12(1-m)):")
for m in (0, 0.15, 0.25, 1/3, 0.375, 0.4, 0.5):
    print(f"  m={m:.4f}  theta={(11-12*m)/(12*(1-m)):.5f}  principal-row excess over U given beta*=7/8: U^{(3*m-1)/12:+.4f}")
