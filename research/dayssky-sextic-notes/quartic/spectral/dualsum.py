# Complete dual sums after Poisson in m:  Ssum(j;c) = sum_{x mod c} S(x^3,x^3;c) psi(j x/c),
# S(r,r';c) = sum_{d mod c}^* chi_c(d) psi((r d + r' dbar)/c)   (odd part of the Kubota Kloosterman sum).
# Prediction (Deligne/Katz, rank-6 Fourier transform sheaf): |Ssum(j;pi)| <= 6 q for every j; moments of
# Ssum/q over j: E|.|^2 -> 1, E|.|^4 -> 2 if G_geom^0 = SL_6 (3 if orthogonal/symplectic).
import numpy as np, math, sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from zmod import ZmodSplit, FieldInert, split_primes
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from gauss4 import gmul
maxN = int(sys.argv[1]) if len(sys.argv) > 1 else 600
rows = []
def stats(name, q, Ss, Sx):
    a = np.abs(Ss)/q
    jn0 = np.arange(len(Ss)) != 0
    m2 = np.mean(a[jn0]**2); m4 = np.mean(a[jn0]**4); m6 = np.mean(a[jn0]**6)
    kmax = np.max(np.abs(Sx[1:]))/math.sqrt(q)
    print(f"{name:>22s} q={q:5d} q%12={q%12:2d}  max|S(x^3,x^3)|/sqrtq={kmax:5.3f}  max_j|Ssum|/q={a.max():6.3f} (j=0: {a[0]:5.3f})"
          f"  E|.|^2={m2:5.3f} E|.|^4={m4:5.3f} E|.|^6={m6:6.3f}  E[Ssum]/q={np.mean(Ss[jn0]).real/q:+.3f}")
    rows.append((q, a.max(), m2, m4, m6))
for pi in split_primes(maxN, 13):
    Z = ZmodSplit(pi); Sx = Z.kloost_cubes(); Ss = Z.dual(Sx)
    stats(f"pi={pi}", Z.N, Ss, Sx)
for p in [7, 11, 19, 23]:
    if p*p > maxN: break
    F = FieldInert(p); Sx = F.kloost_cubes(); Ss = F.dual(Sx)
    stats(f"inert p={p}", F.N, Ss, Sx)
R = np.array(rows)
big = R[R[:, 0] > maxN/3]
print("primes with q>maxN/3: mean E|.|^2=%.3f  mean E|.|^4=%.3f  mean E|.|^6=%.3f   overall max |Ssum|/q=%.3f (Deligne bound 6)"
      % (big[:, 2].mean(), big[:, 3].mean(), big[:, 4].mean(), R[:, 1].max()))
