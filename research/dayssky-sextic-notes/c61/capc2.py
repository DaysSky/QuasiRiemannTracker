# Cycle61-framework theta with generalized slot costs:
#  plain:   n1+n2+cP*kappa*z <= 1  -> zP(m) = (1-2m)/(cP*kappa)   (paper cP = 6)
#  inverse: r + cI*z <= 1          -> zM(r) = (1-r)/cI            (paper cI = 2)
# capacities capped by the supply ell/d.  Numerical piecewise max-min over r.
import numpy as np, sys
from scipy.optimize import brentq
def short(t, dl, p, kap, cP, cI, sup):
    q = dl*p
    lo, hi = max(t - 0.5, 0.0), min(t, 1.0)
    rs = np.linspace(lo, hi, 401)
    AI = 1 - dl*rs - 2*q*np.minimum((1-rs)/cI, sup)
    m = t - rs
    AP = 1 - 2*dl*m - 2*q*np.minimum(np.maximum(1-2*m, 0)/(cP*kap), sup)
    return np.max(np.minimum(AI, AP))
def R(dl, p, kap, cP, cI, sup):
    lng = lambda t: 1 - dl + (5/6 - dl)*(t - 1)
    ts = np.linspace(1.0, 1.5, 101)
    return min(max(short(t, dl, p, kap, cP, cI, sup), lng(t)) for t in ts)
def Gmax(th, h, cP, cI, nd=31, np_=11):
    ell = 11/3 - 4*th; y = h - 2*ell; kap = 2*th - 1
    best = -9
    for d in (0.5, h):
        sup = ell/d
        for dl in np.linspace(0.0, kap, nd):
            for p in np.linspace(0, 0.5, np_):
                Rv = R(dl, p, kap, cP, cI, sup)
                Cd = min(d, 2*y - d)/(cP*kap)
                G = (1+dl)/2 + 13*h/75 - y/2 - ell/2 + ell*dl*p + d*(0.25 + 0.75*Rv - 17/50) - (dl*p/2)*Cd
                best = max(best, G)
    return best
def feasible(th, cP, cI):
    for h in np.linspace(0.70, 0.86, 17):
        if Gmax(th, h, cP, cI) < th: return h
    return None
def theta(cP, cI):
    lo, hi = 0.855, 0.8760
    if feasible(hi, cP, cI) is None: return None
    for _ in range(14):
        mid = (lo+hi)/2
        if feasible(mid, cP, cI) is not None: hi = mid
        else: lo = mid
    return hi
for cP, cI in [(6,2),(4/0.75,2),(6,1),(4/3,2),(4/3,1),(1,1)]:
    print(f"cP={cP:.3f} cI={cI}: theta ~ {theta(cP,cI):.5f}", flush=True)
