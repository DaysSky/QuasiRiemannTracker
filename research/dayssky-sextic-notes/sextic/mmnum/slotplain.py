# Plain fourth moment with a prime slot: Sum_u |S_u(N)^2 P_u(Z)|^2 vs (mean|S^2|^2)(mean|P|^2), rows 6pf.
# Paper cost: 2m + 6*kappa*z <= 1; natural cost-1: 2m + z <= 1.
import sys, math, json, numpy as np
from mmdata import Ideals, rows, run_kernel, W
U = int(float(sys.argv[1])); m = float(sys.argv[2]); zs = [float(x) for x in sys.argv[3].split(',')]
N = U**m; Zmax = U**max(zs)
I = Ideals(int(max(N, Zmax)) + 2); rx, ry = rows(U); nr = len(rx)
sys.argv = ['g', str(int(Zmax) + 2)]
src = open('gbias.py').read().split("# squarefree n via twisted")[0]
ns = {}; exec(src, ns); PR = ns['PR']; gt_p = ns['gt_p']
pidx = {pr['g']: j for j, pr in enumerate(PR)}
isprime = np.zeros(len(I.N), bool); gt_col = np.zeros(len(I.N), complex)
for i in range(1, len(I.N)):
    if I.cof[i] == 0:
        isprime[i] = True; j = pidx.get(I.g[i])
        if j is not None: gt_col[i] = gt_p[j]
rng = np.random.default_rng(5)
specs = []; key = []
lo, hi = I.rng(N/2, N); specs.append((lo, hi, W(I.N[lo:hi]/N))); key.append(('S', None, None))
lo, hi = I.rng(N/2, N); specs.append((lo, hi, I.mu[lo:hi]*W(I.N[lo:hi]/N))); key.append(('M', None, None))
for z in zs:
    Z = U**z; lo, hi = I.rng(Z/2, Z); base = isprime[lo:hi]*W(I.N[lo:hi]/Z)
    for lab, a in (('one', np.ones(hi-lo)), ('gconj', np.conj(gt_col[lo:hi])), ('rand', np.exp(2j*np.pi*rng.random(hi-lo)))):
        specs.append((lo, hi, (base*a).real)); specs.append((lo, hi, (base*a).imag)); key.append(('P', z, lab))
V = run_kernel(I, rx, ry, specs, f"sp{U}")
S = V[:, 0]; S2 = S*S; s4 = np.mean(np.abs(S2)**2); M = V[:, 1]; m2 = np.mean(np.abs(M)**2)
c = 2
print(f"U={U:.0e} m={m} rows={nr} mean|S|^2={np.mean(np.abs(S)**2):.3f} mean|S|^4={s4:.3f} max|S|^4/mean={np.max(np.abs(S2)**2)/s4:.1f}")
for k in key[2:]:
    P = V[:, c] + 1j*V[:, c+1]; c += 2
    p2 = np.mean(np.abs(P)**2)
    rs = np.mean(np.abs(S2*P)**2)/(s4*p2); rm = np.mean(np.abs(M*P)**2)/(m2*p2)
    print(f"  z={k[1]:.3f} a={k[2]:5s} 2m+z={2*m+k[1]:.3f} 2m+4.5z={2*m+4.5*k[1]:.3f}: |S^2P|^2 ratio={rs:.4f}   (|MP|^2 ratio={rm:.4f})", flush=True)
