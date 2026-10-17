import numpy as np, fastRw
Brad = 14200/14361 + 7/10000 - 1e-5
fastRw.B = Brad
import thoptw_core as T
T.B = Brad
def parts(th, h):
    ell = 11/3 - 4*th; y = h - 2*ell
    a = th - 1e-6 - T.Gmax(th, h, ell, y, h, 'new')
    b = th - 1e-5 - T.Gmax(th, h, ell, y, h - 1e-3, 'old')
    c = th - 1e-3 - T.Gmax(th, h, ell, y, 0.5, 'old')
    return a, b, c
lo, hi = 0.87478, 0.87490; hb = None
from scipy.optimize import minimize_scalar
for _ in range(16):
    mid = (lo + hi)/2
    res = minimize_scalar(lambda h: -min(parts(mid, h)), bounds=(0.785, 0.80), method='bounded', options={'xatol':1e-8})
    if -res.fun >= 0: hi = mid; hb = res.x
    else: lo = mid
print(f"theta_min {hi:.8f} h {hb:.6f}  parts at optimum:", [f"{v:+.2e}" for v in parts(hi, hb)], flush=True)
