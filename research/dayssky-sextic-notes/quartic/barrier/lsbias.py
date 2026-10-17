# Analyse the quartic large-sieve runs at rows N(k)<=Q=D^{2/3}, columns D/2<N(c)<=D (Y=D=Q^{3/2}).
# Coefficients: conj(gt) -> row k = sum_c gtilde(conj chi_c, -k) W  (generic index -k: theta bias),
#               rand1   -> random +-1 baseline.
# Reports  S/((Q+Y)||b||^2)  (operator-norm proxy relative to the conjectured optimal constant)
# and the inferred mean-square residue  E/Y^{1/2},  E = S_gen(conj gt)/S_gen(rand) - 1   (~ avg |rho(k)|^2).
import numpy as np, json, glob, math
rows = []
for f in sorted(glob.glob('gout_D*_eta0.6667_n0_s0.npz')):
    d = np.load(f, allow_pickle=True); r = json.loads(str(d['res'])); names = list(d['names'])
    D, H = r['D'], r['H']; w2 = r['sumW2']
    ig, ir = names.index('conj(gt)'), names.index('rand1')
    S = d['S']; Sg = d['S_gen']
    norm = (H + D) * w2
    E = Sg[ig] / Sg[ir] - 1
    rows.append((H, D, S[ig] / norm, S[ir] / norm, E, E / math.sqrt(D)))
rows.sort()
print("   Q      Y=D    S_gt/((Q+Y)|b|^2)  S_rand/((Q+Y)|b|^2)   excess E   E/Y^0.5")
for H, D, a, b, E, e in rows:
    print(f"{H:6d} {D:8d}   {a:10.4f}          {b:10.4f}       {E:8.3f}   {e:.5f}")
if len(rows) >= 2:
    Q = np.log([r[0] for r in rows]); e = np.log([r[5] for r in rows]); x = np.log([r[2] - r[3] for r in rows])
    s1 = np.polyfit(Q, e, 1)[0]; s2 = np.polyfit(Q, x, 1)[0]
    print(f"slope of log(avg|rho|^2) vs log Q : {s1:.3f}   => theta ~ {-s1/2:.3f}  (theta<1/8 <=> non-orthogonality at Y=Q^1.5)")
    print(f"slope of log(excess/((Q+Y)|b|^2)) vs log Q : {s2:.3f}  (prediction 1/4-2theta)")
