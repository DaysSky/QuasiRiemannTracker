# Joint pricing on the record pipeline (Gao's B, Lam): inverse cost CI, LM slope ALPHA, plain weight CPW, LOW saving D.
import sys
sys.path.insert(0, '.')
CI, ALPHA, CPW, D = map(float, sys.argv[1:5])
import fastRc, thcore_c as T
fastRc.CI = CI; fastRc.ALPHA = ALPHA; fastRc.CPW = T.CPW = CPW
fastRc.B = T.B = 0.987047; fastRc.Lam = T.Lam = 3.9586
def parts(th, h):
    ell = 11/3 - 4*th - 4*D; y = h - 2*ell
    return (th - 1e-6 - T.Gmax(th, h, ell, y, h, 'new'), th - 1e-5 - T.Gmax(th, h, ell, y, h - 1e-3, 'old'),
            th - 1e-3 - T.Gmax(th, h, ell, y, 0.5, 'old'))
from scipy.optimize import minimize_scalar
lo, hi = 0.855, 0.8752; hb = None
for _ in range(17):
    mid = (lo + hi)/2
    res = minimize_scalar(lambda h: -min(parts(mid, h)), bounds=(0.74, 0.82), method='bounded', options={'xatol':1e-8})
    if -res.fun >= 0: hi = mid; hb = res.x
    else: lo = mid
print(f"CI={CI} ALPHA={ALPHA:.4f} CPW={CPW:.3f} D={D}: theta_min {hi:.6f} (h {hb:.4f})", flush=True)
