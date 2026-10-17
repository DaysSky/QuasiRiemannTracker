# Cycle61-framework theta as a function of the plain-slot cost coefficient c (z_P = (1-2m)/(c*kappa); paper: c = 6).
import numpy as np
from scipy.optimize import brentq
def short(t, dl, p, kap, c):
    q = dl*p; cI, eI = 1 - q, -(dl - q)
    lo, hi = max(t - 0.5, 0.0), min(t, 1.0)
    a = 1 - 2*dl*t - 2*q*(1 - 2*t)/(c*kap); b = 2*dl - 4*q/(c*kap)
    cands = [lo, hi]
    if abs(eI - b) > 1e-15:
        rc = (a - cI)/(eI - b)
        if lo <= rc <= hi: cands.append(rc)
    return max(min(cI + eI*r, a + b*r) for r in cands)
def R(dl, p, kap, c):
    lng = lambda t: 1 - dl + (5/6 - dl)*(t - 1)
    g = lambda t: short(t, dl, p, kap, c) - lng(t)
    if g(1.0) <= 0: return lng(1.0)
    if g(1.5) > 0: return short(1.5, dl, p, kap, c)
    return lng(brentq(g, 1.0, 1.5, xtol=1e-13))
def Gmax(th, h, c, nd=41, np_=21):
    ell = 11/3 - 4*th; y = h - 2*ell; kap = 2*th - 1
    best = -9
    for d in (0.5, h):
        for dl in np.linspace(0.0, kap, nd):
            for p in np.linspace(0, 0.5, np_):
                Rv = R(dl, p, kap, c)
                Cd = min(d, 2*y - d)/(c*kap)
                G = (1+dl)/2 + 13*h/75 - y/2 - ell/2 + ell*dl*p + d*(0.25 + 0.75*Rv - 17/50) - (dl*p/2)*Cd
                best = max(best, G)
    return best
def feasible(th, c):
    for h in np.linspace(0.70, 0.86, 33):
        if Gmax(th, h, c) < th: return h
    return None
def theta(c):
    lo, hi = 0.86, 0.8760
    if feasible(hi, c) is None: return None
    for _ in range(18):
        mid = (lo+hi)/2
        if feasible(mid, c) is not None: hi = mid
        else: lo = mid
    return hi
for c in (6, 5, 4.5, 4, 3, 2, 1):
    print(f"c={c}: theta ~ {theta(c):.6f}", flush=True)
