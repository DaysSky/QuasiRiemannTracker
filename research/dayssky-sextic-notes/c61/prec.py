# Precise row-count solver: R = min_t max(short(t), long(t)), short(t) = max_r min(AI(r), AP(t-r)).
import numpy as np
from scipy.optimize import brentq, minimize_scalar
Lam = 19312/4787; B = 14200/14361 - 1e-5
def zP(m, kap, mode):
    old = max((1 - 2*m)/(6*kap), 0.0)
    if mode == 'old': return old
    new = (B - 2*m)/Lam
    if 0.41 <= m <= 0.44 and 0.034 <= new <= 0.04 and old >= 0.034: return max(old, new)  # old moment must cover short selected I
    return old
def short(t, dl, p, kap, mode):
    q = dl*p
    lo, hi = max(t - 0.5, 0.0), min(t, 1.0)
    AI = lambda r: 1 - q - (dl - q)*r
    AP = lambda r: 1 - 2*dl*(t - r) - 2*q*zP(t - r, kap, mode)
    f = lambda r: AI(r) - AP(r)
    rs = np.linspace(lo, hi, 2001); vals = [min(AI(r), AP(r)) for r in rs]
    i = int(np.argmax(vals)); best = vals[i]
    # refine with golden search near the max
    a, b = rs[max(i-1,0)], rs[min(i+1,len(rs)-1)]
    res = minimize_scalar(lambda r: -min(AI(r), AP(r)), bounds=(a, b), method='bounded', options={'xatol':1e-12})
    return max(best, -res.fun)
def R(dl, p, kap, mode):
    lng = lambda t: 1 - dl + (5/6 - dl)*(t - 1)
    g = lambda t: short(t, dl, p, kap, mode) - lng(t)
    if g(1.0) <= 0: return lng(1.0) if short(1.0,dl,p,kap,mode) <= lng(1.0) else short(1.0,dl,p,kap,mode)
    if g(1.5) > 0: return short(1.5, dl, p, kap, mode)
    t = brentq(g, 1.0, 1.5, xtol=1e-13)
    return lng(t)
def Rformula(dl, p, kap):
    V = 5/6 - dl
    Ck = 2*(1-p)*(1-p/(3*kap))/(3-p-2*p/(3*kap))
    return 1 - dl + dl*Ck*V/(2*(V + dl*Ck))
if __name__ == "__main__":
    kap = 2*0.874861 - 1
    for dl, p in [(0.475,0.5),(0.47,0.5),(0.48,0.5),(0.45,0.5)]:
        print(dl, p, Rformula(dl,p,kap), R(dl,p,kap,'old'), R(dl,p,kap,'new'))
