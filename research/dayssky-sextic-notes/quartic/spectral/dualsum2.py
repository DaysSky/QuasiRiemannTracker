# Dual complete sums for composite squarefree c = pi1*pi2 and for c = pi^2 (Kuznetsov's c-sum runs over all c).
# Also: is Ssum(j;pi) real up to a fixed phase?  (structure test)
import numpy as np, math, sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from zmod import ZmodSplit, split_primes
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from gauss4 import gmul
P = split_primes(120, 13)
print("composite c = pi1 pi2:  Deligne+CRT bound |Ssum| <= 36 N(c);  also |S(x^3,x^3;c)| <= 4 sqrt(N c)")
for i in range(len(P)):
    for k in range(i+1, len(P)):
        c = gmul(P[i], P[k]); N = c[0]**2+c[1]**2
        if N > 2600 or math.gcd(*c) != 1: continue
        Z = ZmodSplit(c); Sx = Z.kloost_cubes(); Ss = Z.dual(Sx)
        units = [x for x in range(N) if math.gcd(x, N) == 1]
        print(f"  c={c} N={N:5d}: max|S(x^3,x^3)|/sqrtN={np.max(np.abs(Sx))/math.sqrt(N):5.3f}  max|Ssum|/N={np.max(np.abs(Ss))/N:6.3f}"
              f"  E|Ssum/N|^2={np.mean(np.abs(Ss/N)**2):5.3f}")
        if len(P) > 6 and k > i+3: break
print("prime squares c = pi^2 (kappa = quadratic char mod pi: Salie-type):")
for pi in split_primes(45, 5):
    c = gmul(pi, pi); N = c[0]**2+c[1]**2
    Z = ZmodSplit(c); Sx = Z.kloost_cubes(); Ss = Z.dual(Sx)
    print(f"  pi={pi} N(c)={N:5d}: max|S(x^3,x^3)|/sqrtN={np.max(np.abs(Sx))/math.sqrt(N):5.3f}  max|Ssum|/N={np.max(np.abs(Ss))/N:6.3f}"
          f"  E|Ssum/N|^2={np.mean(np.abs(Ss/N)**2):5.3f}  #nonzero j / N = {np.mean(np.abs(Ss)>1e-6*N):.3f}")
print("phase structure of Ssum(j;pi): spread of arg(Ssum^2) over j (0 => real up to fixed phase)")
for pi in split_primes(400, 300):
    Z = ZmodSplit(pi); Sx = Z.kloost_cubes(); Ss = Z.dual(Sx)
    w = Ss[1:]**2; w = w[np.abs(w) > 1e-6*Z.N**2]
    print(f"  q={Z.N}: |mean(Ssum^2/|Ssum|^2)| = {abs(np.mean(w/np.abs(w))):.3f}   real-Gaussian check E|.|^4/(E|.|^2)^2 = "
          f"{np.mean(np.abs(Ss[1:])**4)/np.mean(np.abs(Ss[1:])**2)**2:.3f}")
