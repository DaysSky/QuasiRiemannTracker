import numpy as np, sys
from c61model import Gmax, low
from scipy.optimize import minimize
def margin(v, th, cap=True, capcost=None):
    h, ell, y = v
    x = 1 + ell - h
    M = x + y
    L = low(ell, x, y)
    if capcost is not None and M + ell > 1: L += capcost*(M + ell - 1)
    G, _ = Gmax(th, h, ell, y, nd=41, np_=41)
    pen = 0
    cons = [x - ell, y - x, 2*x - y - ell, 4*ell - y, 0.2 - ell, h - 4*ell]
    if cap and capcost is None: cons.append(1 - M - ell)
    pen = sum(max(0, -c) for c in cons)*10
    return max(L - th, G - th) + pen
def feasible(th, **kw):
    best = 1e9
    for h0 in np.linspace(0.76, 0.83, 8):
        for l0 in np.linspace(0.16, 0.19, 7):
            for yoff in np.linspace(-0.04, 0.06, 6):
                y0 = h0 - 2*l0 + yoff
                m = margin((h0, l0, y0), th, **kw)
                if m < best: best, arg = m, (h0, l0, y0)
    res = minimize(lambda v: margin(v, th, **kw), arg, method='Nelder-Mead', options={'xatol':1e-7,'fatol':1e-9,'maxiter':400})
    return res.fun, res.x
def bisect(**kw):
    lo, hi = 0.85, 0.876
    for _ in range(18):
        mid = (lo+hi)/2
        f, v = feasible(mid, **kw)
        if f <= 0: hi = mid; vb = v
        else: lo = mid
    return hi, vb
if __name__ == "__main__":
    mode = sys.argv[1]
    if mode == "cap": print("with cap", bisect(cap=True))
    elif mode == "nocap": print("no cap", bisect(cap=False))
    else: print("capcost", mode, bisect(cap=True, capcost=float(mode)))
