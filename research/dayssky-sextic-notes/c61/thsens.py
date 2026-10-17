import numpy as np, sys, fastRw
import thoptw_core as T
from scipy.optimize import minimize_scalar
def run(Bv, Lv):
    fastRw.B = Bv; T.B = Bv; fastRw.Lam = Lv; T.Lam = Lv
    def parts(th, h):
        ell = 11/3 - 4*th; y = h - 2*ell
        return min(th - 1e-6 - T.Gmax(th, h, ell, y, h, 'new'), th - 1e-5 - T.Gmax(th, h, ell, y, h - 1e-3, 'old'),
                   th - 1e-3 - T.Gmax(th, h, ell, y, 0.5, 'old'))
    lo, hi = 0.8746, 0.87490
    for _ in range(17):
        mid = (lo + hi)/2
        res = minimize_scalar(lambda h: -parts(mid, h), bounds=(0.785, 0.80), method='bounded', options={'xatol':1e-8})
        if -res.fun >= 0: hi = mid
        else: lo = mid
    return hi
B0 = 14200/14361; L0 = 19312/4787
for dB, dL in [(7e-4, 0), (17e-4, 0), (27e-4, 0), (7e-4, -0.05), (7e-4, -0.2), (7e-4, -0.5)]:
    print(f"B = B0 + {dB:.4f}, Lambda = {L0+dL:.4f}: theta = {run(B0 + dB - 1e-5, L0 + dL):.8f}", flush=True)
