import numpy as np, sys
from fastR import R, Lam, B
from scipy.optimize import minimize, minimize_scalar
mode = sys.argv[1]
def G(d, dl, p, th, h, ell, y, Cmode):
    kap = 2*th - 1
    Rv = R(dl, p, kap, mode)
    Cd = (2*y - (2-B)*d)/Lam if Cmode == 'new' else min(d, 2*y - d)/(6*kap)
    return (1+dl)/2 + 13*h/75 - y/2 - ell/2 + ell*dl*p + d*(0.25 + 0.75*Rv - 17/50) - (dl*p/2)*Cd
def Gmax(th, h, ell, y, d, Cmode):
    kap = 2*th - 1; best = (-9, None)
    for dl in np.linspace(0, kap, 61):
        for p in np.linspace(0, 0.5, 21):
            g = G(d, dl, p, th, h, ell, y, Cmode)
            if g > best[0]: best = (g, (dl, p))
    dl0, p0 = best[1]
    f = lambda v: -G(d, min(max(v[0],0),kap), min(max(v[1],0),0.5), th, h, ell, y, Cmode)
    res = minimize(f, [dl0, p0], method='Nelder-Mead', options={'xatol':1e-10,'fatol':1e-13,'maxiter':600})
    return max(best[0], -res.fun)
def margin(th, h):
    ell = 11/3 - 4*th; y = h - 2*ell; x = 1 + ell - h
    if not (0 < ell < x < y < 2*x and y < 0.5 < h < 2*y and 1/6 < ell < 1/5 and h > 4*ell and 4*ell - y > 0): return -1
    g1 = Gmax(th, h, ell, y, h, 'new'); g2 = Gmax(th, h, ell, y, h - 1e-3, 'old'); g3 = Gmax(th, h, ell, y, 0.5, 'old')
    return min(th - 1e-6 - g1, th - 1e-5 - g2, th - 1e-3 - g3)
def best_h(th):
    res = minimize_scalar(lambda h: -margin(th, h), bounds=(0.785, 0.80), method='bounded', options={'xatol':1e-8})
    return -res.fun, res.x
m, h = best_h(0.874861); print(mode, "margin at record theta:", m, "h", h, flush=True)
lo, hi = 0.87475, 0.874861; hb = h
for _ in range(16):
    mid = (lo+hi)/2
    m, h = best_h(mid)
    if m >= 0: hi = mid; hb = h
    else: lo = mid
print(mode, "theta_min", hi, "h", hb, "ell", 11/3-4*hi, flush=True)
