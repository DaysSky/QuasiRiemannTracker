# Plug Theorem-2.2 capacity (with / without amplitude classes) into the Cycle61-framework theta model.
import numpy as np
from ampclass import feasible
from capc3 import ND, NP, NT, NR
def capcurve(amp):
    ms = np.linspace(0.36, 0.48, 25); caps = []
    for m in ms:
        best = 0
        for q0 in np.linspace(0.55, 0.74, 20):
            lo, hi = 0.0, 0.2
            for _ in range(30):
                mid = (lo+hi)/2
                if feasible(2*m, mid, q0, amp): lo = mid
                else: hi = mid
            best = max(best, lo)
        caps.append(best)
    return ms, np.array(caps)
def Rgrid(kap, ms, caps, sup):
    dl = np.linspace(0, kap, ND)[:,None,None,None]; p = np.linspace(0, .5, NP)[None,:,None,None]
    t = np.linspace(1, 1.5, NT)[None,None,:,None]; a = np.linspace(0, 1, NR)[None,None,None,:]
    lo = np.maximum(t - .5, 0); hi = np.minimum(t, 1.0); r = lo + (hi - lo)*a
    q = dl*p; m = t - r
    AI = 1 - dl*r - 2*q*np.minimum((1 - r)/2, sup)
    zold = np.maximum(1 - 2*m, 0)/(6*kap)
    znew = np.interp(m, ms, caps, left=0, right=0)
    AP = 1 - 2*dl*m - 2*q*np.minimum(np.maximum(zold, znew), sup)
    short = np.minimum(AI, AP).max(axis=3)
    lng = (1 - dl + (5/6 - dl)*(t - 1))[..., 0]
    return dl[:,:,0,0], p[:,:,0,0], np.maximum(short, lng).min(axis=2)
def theta(ms, caps):
    def Gmax(th, h):
        ell = 11/3 - 4*th; y = h - 2*ell; kap = 2*th - 1; best = -9
        for d in (0.5, h):
            dl, p, R = Rgrid(kap, ms, caps, ell/d)
            Cd = min(d, 2*y - d)/(6*kap)
            G = (1+dl)/2 + 13*h/75 - y/2 - ell/2 + ell*dl*p + d*(0.25 + 0.75*R - 17/50) - (dl*p/2)*Cd
            best = max(best, G.max())
        return best
    def feas(th):
        for h in np.linspace(0.75, 0.84, 46):
            if Gmax(th, h) < th: return h
        return None
    lo, hi = 0.872, 0.8765
    for _ in range(18):
        mid = (lo+hi)/2
        if feas(mid) is not None: hi = mid
        else: lo = mid
    return hi
for amp in (False, True):
    ms, caps = capcurve(amp)
    print("amplitude classes" if amp else "Theorem 2.2      ", "theta ~", round(theta(ms, caps), 7), flush=True)
ms = np.linspace(0.36,0.48,25); print("old 6kappa only  theta ~", round(theta(ms, np.zeros(25)), 7))
