import numpy as np, sys
from precopt import Gmax
from scipy.optimize import minimize_scalar
mode = sys.argv[1]
def margin(th, h):
    ell = 11/3 - 4*th; y = h - 2*ell; x = 1 + ell - h
    if not (0 < ell < x < y < 2*x and y < 0.5 < h < 2*y and 1/6 < ell < 1/5 and h > 4*ell and 4*ell - y > 0): return -1
    g1, _ = Gmax(th, h, ell, y, mode, h, 'new')
    g2, _ = Gmax(th, h, ell, y, mode, h - 1e-3, 'old')
    g3, _ = Gmax(th, h, ell, y, mode, 0.5, 'old')
    return min(th - 1e-6 - g1, th - 1e-5 - g2, th - 1e-3 - g3)
def best_h(th):
    res = minimize_scalar(lambda h: -margin(th, h), bounds=(0.78, 0.81), method='bounded', options={'xatol':1e-7})
    return -res.fun, res.x
lo, hi = 0.8745, 0.874861
m, h = best_h(hi); print("at record theta:", m, h, flush=True)
for _ in range(20):
    mid = (lo+hi)/2
    m, h = best_h(mid)
    if m >= 0: hi = mid; hb = h
    else: lo = mid
print(mode, "theta_min", hi, "h", hb, "ell", 11/3-4*hi, flush=True)
