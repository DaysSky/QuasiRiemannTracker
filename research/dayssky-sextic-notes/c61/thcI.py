import sys, fastRwc as F
F.CI = float(sys.argv[1]); lo0 = float(sys.argv[2])
import os
F.B = float(os.environ.get('GB', 14200/14361 + 7/10000 - 1e-5))
if 'GL' in os.environ: F.Lam = float(os.environ['GL'])
sys.modules['fastRw'] = F
import thoptw_core as T
T.B = F.B; T.R = F.R
if 'GL' in os.environ: T.Lam = float(os.environ['GL'])
from scipy.optimize import minimize_scalar
def parts(th, h):
    ell = 11/3 - 4*th; y = h - 2*ell; x = 1 + ell - h
    if not (0 < ell < x < y < 2*x and y < 0.5 < h < 2*y and 1/6 < ell < 1/5 and h > 4*ell and 4*ell - y > 0): return (-1,)
    return (th - 1e-6 - T.Gmax(th, h, ell, y, h, 'new'), th - 1e-5 - T.Gmax(th, h, ell, y, h - 1e-3, 'old'),
            th - 1e-3 - T.Gmax(th, h, ell, y, 0.5, 'old'))
lo, hi = lo0, 0.87490; hb = None
for _ in range(14):
    mid = (lo + hi)/2
    res = minimize_scalar(lambda h: -min(parts(mid, h)), bounds=(0.74, 0.80), method='bounded', options={'xatol':1e-7})
    if -res.fun >= 0: hi = mid; hb = res.x
    else: lo = mid
print(f"c_I={F.CI} theta_min {hi:.7f} h {hb}", flush=True)
