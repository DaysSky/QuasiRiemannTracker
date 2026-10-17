# Radical-record pipeline (detector-affine + timbalud B) with inverse cost CI, slot supply cap, geometry checks,
# the floor bin (theta > 51/100 + h/3 + 51 ell/100), and the affine HIGH capacity C_N at d = h used only when the
# detector-affine window conditions (3.3) hold there; otherwise C_O at d = h.
import sys
sys.path.insert(0, '.')
CONST = sys.argv[1]; CI, ALPHA, CPW = map(float, sys.argv[2:5]); hlo, hhi = 0.74, 0.82
import fastRs, thcore_s as T
fastRs.CI = CI; fastRs.ALPHA = ALPHA; fastRs.CPW = T.CPW = CPW
WIDE = CPW < 3.9586
if WIDE: B, LAM = 1.0, CPW; fastRs.WIN.update(ZLO=0.001, ZHI=0.30, NLO=0.0, NHI=0.5)
elif CONST == 'gao': B, LAM = 0.987047, 3.9586
else: B, LAM = 14200/14361 + 7/10000 - 1e-5, 19312/4787
fastRs.B = T.B = B; fastRs.Lam = T.Lam = LAM
def window_ok(th, h, ell, y):
    if WIDE: return True
    dm, dp, t0 = h - 1e-3, h + 1e-10, 1e-3
    n_lo = 1 - y/dm; n_hi = 1 - y/dp + t0
    zA = lambda n: (B - 2*n)/LAM
    return 0.41 < n_lo and n_hi < 264/625 and zA(n_hi) > 17/500 and zA(n_lo) < 1/25
def parts(th, h):
    ell = 11/3 - 4*th; y = h - 2*ell; x = 1 + ell - h
    if not (0 < ell < x < y < 2*x and y < 0.5 < h < 2*y and ell < 1/5 and h > 4*ell and 4*ell - y > 0): return (-1,)
    if not th - (0.51 + h/3 + 0.51*ell) > 1e-4: return (-1,)
    fastRs.LSUP = ell/h
    mode = 'new' if window_ok(th, h, ell, y) else 'old'
    return (th - 1e-5 - T.Gmax(th, h, ell, y, h, mode), th - 1e-5 - T.Gmax(th, h, ell, y, h - 1e-3, 'old'),
            th - 1e-3 - T.Gmax(th, h, ell, y, 0.5, 'old'))
from scipy.optimize import minimize_scalar
import numpy as np
lo, hi = 0.855, 0.8752; hb = None
for _ in range(17):
    mid = (lo + hi)/2
    best = max((min(parts(mid, h)), h) for h in np.linspace(hlo, hhi, 81))
    res = minimize_scalar(lambda h: -min(parts(mid, h)), bounds=(max(hlo, best[1]-0.002), min(hhi, best[1]+0.002)), method='bounded', options={'xatol':1e-8})
    v, hv = max((best[0], best[1]), (-res.fun, res.x))
    if v >= 0: hi = mid; hb = hv
    else: lo = mid
print(f"[window-aware {CONST}] CI={CI} ALPHA={ALPHA:.4f} CPW={CPW:.3f}: theta_min {hi:.6f} (h {hb})", flush=True)
