# Check the hypergeometric (rank-6) Mellin form of the dual sum:
#   Ssum(j;pi) = (1/(q-1)) sum_lambda (conj(chi) lambda^6)(j) g(lambda) g(lambda chi) g(chi lambda^{-6}),  g(1) := -1,
# and that Ssum is REAL (S(x^3,x^3) at -x is the conjugate).
import numpy as np, math, cmath, sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from zmod import ZmodSplit, split_primes
for pi in [split_primes(37, 37)[0], split_primes(61, 61)[0], split_primes(109, 109)[0]]:
    Z = ZmodSplit(pi); N = Z.N
    g0 = next(g for g in range(2, N) if all(pow(g, (N-1)//t, N) != 1 for t in set(f for f in range(2, N) if (N-1) % f == 0 and all(f % s for s in range(2, f)))))
    dl = np.zeros(N, np.int64); v = 1
    for k in range(N-1): dl[v] = k; v = v*g0 % N
    x = np.arange(1, N); ex = Z.e(x)
    kchi = int(round(np.angle(Z.chi[g0])/(2*np.pi)*(N-1))) % (N-1)    # chi = lambda_{kchi}
    def g(k):
        k %= (N-1)
        return -1.0 if k == 0 else np.sum(np.exp(2j*np.pi*k*dl[x]/(N-1))*ex)
    G = {k: g(k) for k in range(N-1)}
    Sx = Z.kloost_cubes(); Ss = Z.dual(Sx)
    err = 0
    for j in [1, 2, 3, 7, N-1]:
        tot = 0
        for k in range(N-1):
            eta = (kchi - 6*k) % (N-1)            # chi * lambda^{-6}
            tot += np.exp(-2j*np.pi*eta*dl[j]/(N-1))*G[k]*G[(k+kchi) % (N-1)]*G[eta]
        err = max(err, abs(tot/(N-1) - Ss[j]))
    print(f"q={N}: max |Mellin formula - Ssum| = {err:.2e}   max |Im Ssum|/q = {np.max(np.abs(Ss.imag))/N:.2e}")
