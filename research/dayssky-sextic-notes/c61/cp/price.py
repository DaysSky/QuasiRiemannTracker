# Price of a plain slot weight CPW (in place of 6*kappa) on the record pipeline with Gao's constants.
# Mode 'gao': keep Gao's affine branch (B, Lam) in its window, change only the old-branch weight.
# Mode 'all': slot weight CPW everywhere (Lam = CPW, B = 1, wide window).
import sys
sys.path.insert(0, '.')
CPW = float(sys.argv[1]); mode = sys.argv[2]
import fastRw, thoptw_core as T
fastRw.CPW = CPW; T.CPW = CPW
if mode == 'gao':
    fastRw.B = T.B = 0.987047; fastRw.Lam = T.Lam = 3.9586
else:
    fastRw.B = T.B = 1.0; fastRw.Lam = T.Lam = CPW
    fastRw.WIN.update(ZLO=0.001, ZHI=0.30, NLO=0.0, NHI=0.5)
def parts(th, h):
    ell = 11/3 - 4*th; y = h - 2*ell
    return (th - 1e-6 - T.Gmax(th, h, ell, y, h, 'new'), th - 1e-5 - T.Gmax(th, h, ell, y, h - 1e-3, 'old'),
            th - 1e-3 - T.Gmax(th, h, ell, y, 0.5, 'old'))
from scipy.optimize import minimize_scalar
lo, hi = 0.860, 0.8752; hb = None
for _ in range(18):
    mid = (lo + hi)/2
    res = minimize_scalar(lambda h: -min(parts(mid, h)), bounds=(0.76, 0.82), method='bounded', options={'xatol':1e-8})
    if -res.fun >= 0: hi = mid; hb = res.x
    else: lo = mid
print(f"slot weight {CPW:.3f} mode={mode}: theta_min {hi:.6f} (h {hb})", flush=True)
