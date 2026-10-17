import json, glob, math, sys
import numpy as np
files = sorted(glob.glob("res_U*_s12345.json"), key=lambda f: int(f.split('_U')[1].split('_')[0]))
R = [json.load(open(f)) for f in files]
Us = [r['U'] for r in R]
print("U values:", Us, " rows:", [r['rows'] for r in R])
keys = []
for c in R[-1]['configs']:
    k = (c['r'], c['m'], c.get('m2'), c['var'])
    if k not in keys: keys.append(k)
def get(r, k):
    for c in r['configs']:
        if (c['r'], c['m'], c.get('m2'), c['var']) == k: return c
fitU = lambda: [i for i, U in enumerate(Us) if U >= 1e4]
print(f"{'r':>6} {'m':>6} {'var':>5} | " + " ".join(f"{'ratio@'+format(U,'.0e'):>11}" for U in Us)
      + " | z@maxU | diag@U(min..max) | lam_MM(U>=1e4) lam_ratio")
for k in keys:
    r_, m_, m2, v = k
    cs = [get(r, k) for r in R]
    if any(c is None for c in cs): continue
    if cs[-1]['diag'] == 0 or any(not np.isfinite(c['ratio']) for c in cs if c['diag'] > 0): pass
    rat = [c['ratio'] if c['diag'] > 0 else float('nan') for c in cs]
    idx = [i for i in fitU() if cs[i]['MM'] > 0]
    lamMM = np.polyfit(np.log([Us[i] for i in idx]), np.log([cs[i]['MM'] for i in idx]), 1)[0] if len(idx) > 1 else float('nan')
    lamR = np.polyfit(np.log([Us[i] for i in idx]), np.log([cs[i]['ratio'] for i in idx]), 1)[0] if len(idx) > 1 else float('nan')
    lab = f"{'plain' if r_ is None else format(r_, '.3f'):>6} {m_:6.3f} {v if m2 is None else 'x'+str(m2):>5}"
    print(lab + " | " + " ".join(f"{x:11.4f}" for x in rat) + f" | {cs[-1].get('z', float('nan')):+5.2f} | "
          f"{cs[0]['diag']:.4f}..{cs[-1]['diag']:.4f} | {lamMM:+.4f} {lamR:+.5f}")
# aggregate deviation statistics
print("\nRMS(ratio-1) and mean(z) over non-degenerate configs (nS>=5 or m=0), by variant and U:")
for v in ('mu', 'rnd', 'rmul'):
    line = []
    for r in R:
        cs = [c for c in r['configs'] if c['var'] == v and c['diag'] > 0 and (c['m'] == 0 or c.get('nS', 0) >= 5)]
        dev = np.array([c['ratio'] - 1 for c in cs]); z = np.array([c['z'] for c in cs])
        line.append(f"U={r['U']:.0e}: rms={np.sqrt(np.mean(dev**2)):.4f} meanz={z.mean():+.2f} n={len(cs)}")
    print(f"  {v:5s} " + " | ".join(line))
print("\nby total length r+m (mu only), mean ratio over U>=1e4 & non-degenerate:")
for tot in (0.9, 1.0, 1.1, 1.25):
    vals = [c['ratio'] for r in R if r['U'] >= 1e4 for c in r['configs']
            if c['var'] == 'mu' and c.get('tot') == tot and c['diag'] > 0 and (c['m'] == 0 or c.get('nS', 0) >= 5)]
    print(f"  r+m={tot}: mean ratio={np.mean(vals):.4f} sd={np.std(vals):.4f} n={len(vals)}")
print("\nmax_u |M_u S_u|^2 / U at largest U (trivial bound ~ DN/U = U^{r+m-1}):")
for c in R[-1]['configs']:
    if c['var'] == 'mu' and c['diag'] > 0:
        print(f"  r={c['r']:.3f} m={c['m']:.3f}: max/U={c['maxrow_over_U']:.2e}  share of S-supported rows={c['bounded_rows_share']:.2e}")
