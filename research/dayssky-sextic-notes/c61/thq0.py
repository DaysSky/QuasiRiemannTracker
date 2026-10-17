import sys, numpy as np, fastRw
from fractions import Fraction as F
sys.path.insert(0,'../q0scan')
from fam import BL
q0 = F(sys.argv[1]); B0q, Lq, Ka, r0 = BL(q0)
Brad = float(B0q) + 7/10000 - 1e-5; Lam = float(Lq)
fastRw.B = Brad; fastRw.Lam = Lam
import thoptw_core as T
T.B = Brad; T.Lam = Lam
def parts(th, h):
    ell = 11/3 - 4*th; y = h - 2*ell
    return (th - 1e-6 - T.Gmax(th, h, ell, y, h, 'new'), th - 1e-5 - T.Gmax(th, h, ell, y, h - 1e-3, 'old'),
            th - 1e-3 - T.Gmax(th, h, ell, y, 0.5, 'old'))
from scipy.optimize import minimize_scalar
lo, hi = 0.87476, 0.87490; hb = None
for _ in range(15):
    mid = (lo + hi)/2
    res = minimize_scalar(lambda h: -min(parts(mid, h)), bounds=(0.785, 0.80), method='bounded', options={'xatol':1e-8})
    if -res.fun >= 0: hi = mid; hb = res.x
    else: lo = mid
print(f"q0={float(q0):.3f} B={Brad:.6f} Lam={Lam:.5f} theta_min {hi:.8f} h {hb}", flush=True)
