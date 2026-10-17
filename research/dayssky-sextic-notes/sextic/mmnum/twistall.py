# Twisted plain fourth moments T(psi) = mean_v |sum_{n1,n2} w w psi(n1 n2) chi_{n1 n2}(v)|^2 / X^2 over ALL characters psi of (O/m)^*,
# for small split primes m. Sextic-type psi (chi_n(m)^j, orders dividing 6) are the thin families (row shifts); the others are "foreign" twists.
# Parseval: average over all psi of T(psi) = (p-1) * sum_c mean|P_c|^2 / ... ; the question is whether T(psi) is natural for every psi.
import sys, numpy as np, cmath
from mmdata import Ideals, rows, run_kernel, W
U = int(float(sys.argv[1])); n = float(sys.argv[2]); X = U**n
I = Ideals(int(X) + 2); rx, ry = rows(U)
lo, hi = I.rng(X/2, X); w = W(I.N[lo:hi]/X); gens = I.g[lo:hi]
ms = [pr for pr in I.primes if pr['kind'] == 'split' and pr['P'] in (7, 13, 19)][::2]
specs = []; meta = []
for pr in ms:
    p, r = pr['P'], pr['r']
    cls = np.array([(a + b*r) % p for (a, b) in gens])
    for c in range(1, p):
        specs.append((lo, hi, w*(cls == c))); meta.append((p, c))
S = run_kernel(I, rx, ry, specs, f"twall{U}_{n}") / np.sqrt(X)
full = run_kernel(I, rx, ry, [(lo, hi, w)], f"twallf{U}_{n}")[:, 0] / np.sqrt(X)
base = np.mean(np.abs(full)**4)
print(f"U={U:.0e} X=U^{n} rows={len(rx)} columns={hi-lo} mean|S^2|^2={base:.4f}")
k0 = 0
for pr in ms:
    p = pr['P']; Sc = S[:, k0:k0+p-1]; k0 += p-1
    g = next(x for x in range(2, p) if all(pow(x, (p-1)//q, p) != 1 for q in range(2, p) if (p-1) % q == 0 and all(q % d for d in range(2, q))))
    dlog = {pow(g, j, p): j for j in range(p-1)}
    # product sums P_c = sum_{c1 c2 = c} S_c1 S_c2
    P = np.zeros_like(Sc)
    for c1 in range(1, p):
        for c2 in range(1, p):
            P[:, (c1*c2 % p) - 1] += Sc[:, c1-1]*Sc[:, c2-1]
    res = []
    for k in range(p-1):
        psi = np.array([cmath.exp(2j*cmath.pi*k*dlog[c]/(p-1)) for c in range(1, p)])
        T = np.mean(np.abs(P @ psi)**2)
        order = (p-1)//np.gcd(k, p-1)
        res.append((k, order, T/base))
    avg = np.mean([t for _, _, t in res])
    line = "  ".join(f"k={k}(ord {o}):{t:.3f}{'*' if 6 % o == 0 else ''}" for k, o, t in res)
    print(f"N(m)={p}: T(psi)/T(1-free)  [* = sextic-type, i.e. thin family]  average {avg:.3f}\n   {line}", flush=True)
