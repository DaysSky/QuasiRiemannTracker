# Restricted Gram structure for large-value row sets of the Mobius polynomial M_u(D), D = U^r.
# M_B[u, n] = chi_n(u) * W(N(n)/D) for squarefree n in (D/2, D];  Gram K = M_B M_B^*.
import sys, math, json, numpy as np
from mmdata import Ideals, rows, run_kernel, W
U = int(float(sys.argv[1])); r = float(sys.argv[2]); sizes = [int(s) for s in sys.argv[3].split(',')]
seed = 12345; rng = np.random.default_rng(seed)
D = U**r; I = Ideals(int(D)+2); rx, ry = rows(U)
lo, hi = I.rng(D/2, D)
sq = I.mu[lo:hi] != 0
wts = W(I.N[lo:hi]/D)*sq
# 1) Mobius polynomial values for all rows
V = run_kernel(I, rx, ry, [(lo, hi, I.mu[lo:hi]*W(I.N[lo:hi]/D))], f"g{U}")[:,0]/math.sqrt(D)
a2 = np.abs(V)**2
order = np.argsort(-a2)
idx_n = np.arange(lo, hi)[sq]; w = wts[sq]
Ncol = len(idx_n); Weff = float((w**2).sum())
res = dict(U=U, r=r, D=D, rows=len(rx), ncols=Ncol, Weff=Weff, mean_a2=float(a2.mean()), out=[])
def gram_stats(sel, label):
    # matrix of chi_n(u) for u in sel, n in idx_n (one spec per column)
    specs = [(int(i), int(i)+1, np.array([wi])) for i, wi in zip(idx_n, w)]
    M = run_kernel(I, rx[sel], ry[sel], specs, f"gb{U}{label}")   # |B| x Ncol
    K = M @ M.conj().T
    d = np.real(np.diag(K)); Koff = K - np.diag(np.diag(K))
    s = np.linalg.svd(M, compute_uv=False)
    S2 = float(np.sum(np.abs(Koff)**2))
    S3 = complex(np.trace(Koff @ Koff @ Koff))
    B = len(sel)
    return dict(label=label, B=B, s1sq=float(s[0]**2), s2sq=float(s[1]**2), s5sq=float(s[min(4,len(s)-1)]**2),
                diag_mean=float(d.mean()), S2=S2, S2_over_B2W=S2/(B*B*Weff),
                S3_re=S3.real, S3_im=S3.imag, S3_over_rand=abs(S3)/((B*Weff)**1.5),
                coh_pred=(Weff**2)*B**3/U,
                minA2=float(a2[sel].min()) if label.startswith('top') else None)
for b in sizes:
    for label, sel in (('top', order[:b]), ('rand', rng.choice(len(rx), b, replace=False))):
        st = gram_stats(sel, f"{label}{b}")
        res['out'].append(st); print(json.dumps(st), flush=True)
json.dump(res, open(f"gramB_U{U}_r{r}.json", "w"), indent=1)
