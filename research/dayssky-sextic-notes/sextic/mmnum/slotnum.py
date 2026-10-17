# Marked moment numerics: Sum_u |M_u(N) P_u(Z)|^2 vs (Sum|M|^2)(Sum|P|^2)/rows, with slot coefficients a(p):
#  'one' a=1, 'gconj' a = conj(gtilde(p)) (Gauss-correlated: dual-coherent), 'rand' random phases.
# Cost-2 lemma is uniform in a(p); cost-1 conjecture predicts ratio ~ 1 when N Z <= U < N Z^2 (for non-Gauss a).
import sys, math, json, numpy as np
from mmdata import Ideals, rows, run_kernel, W
U = int(float(sys.argv[1])); r = float(sys.argv[2]); zs = [float(x) for x in sys.argv[3].split(',')]
N = U**r; Zmax = U**max(zs)
I = Ideals(int(max(N, Zmax)) + 2); rx, ry = rows(U); nr = len(rx)
# Gauss sums at primes (reuse gbias machinery, small primes only)
sys.argv = ['g', str(int(Zmax) + 2)]
src = open('gbias.py').read().split("# squarefree n via twisted")[0]
ns = {}; exec(src, ns)
PR = ns['PR']; gt_p = ns['gt_p']
pidx = {pr['g']: j for j, pr in enumerate(PR)}
isprime = np.zeros(len(I.N), bool); gt_col = np.zeros(len(I.N), complex)
for i in range(1, len(I.N)):
    if I.cof[i] == 0:
        isprime[i] = True
        j = pidx.get(I.g[i]);
        if j is not None: gt_col[i] = gt_p[j]
rng = np.random.default_rng(3)
specs = []; key = []
lo, hi = I.rng(N/2, N); specs.append((lo, hi, I.mu[lo:hi]*W(I.N[lo:hi]/N))); key.append(('M', None, None))
for z in zs:
    Z = U**z; lo, hi = I.rng(Z/2, Z)
    base = isprime[lo:hi]*W(I.N[lo:hi]/Z)
    for lab, a in (('one', np.ones(hi-lo)), ('gconj', np.conj(gt_col[lo:hi])), ('rand', np.exp(2j*np.pi*rng.random(hi-lo)))):
        # kernel takes real weights: split complex into re/im specs
        specs.append((lo, hi, (base*a).real)); specs.append((lo, hi, (base*a).imag)); key.append(('P', z, lab))
V = run_kernel(I, rx, ry, specs, f"slot{U}")
M = V[:, 0]; m2 = np.mean(np.abs(M)**2)
out = []; c = 1
for k in key[1:]:
    P = V[:, c] + 1j*V[:, c+1]; c += 2
    p2 = np.mean(np.abs(P)**2); mp = np.mean(np.abs(M*P)**2)
    rec = dict(z=k[1], a=k[2], ratio=float(mp/(m2*p2)), NZ=float(r+k[1]), NZ2=float(r+2*k[1]))
    out.append(rec)
    print(f"U={U:.0e} r={r} z={k[1]:.3f} a={k[2]:5s} r+z={r+k[1]:.3f} r+2z={r+2*k[1]:.3f}  mean|MP|^2/(mean|M|^2 mean|P|^2) = {rec['ratio']:.4f}", flush=True)
json.dump(dict(U=U, r=r, rows=nr, res=out), open(f"slot_U{U}_r{r}.json", "w"), indent=1)
