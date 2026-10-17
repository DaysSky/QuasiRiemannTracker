# Vectorized Cycle61-framework theta with generalized slot costs (plain cP*kappa, inverse cI).
# plain:   n1+n2+cP*kappa*z <= 1  -> zP(m) = (1-2m)/(cP*kappa)   (paper cP = 6)
# inverse: r + cI*z <= 1          -> zM(r) = (1-r)/cI            (paper cI = 2)
import numpy as np, sys
ND, NP, NT, NR = 41, 11, 51, 121
def Rgrid(kap, cP, cI, sup):
    dl = np.linspace(0, kap, ND)[:,None,None,None]; p = np.linspace(0, .5, NP)[None,:,None,None]
    t = np.linspace(1, 1.5, NT)[None,None,:,None]; a = np.linspace(0, 1, NR)[None,None,None,:]
    lo = np.maximum(t - .5, 0); hi = np.minimum(t, 1.0); r = lo + (hi - lo)*a
    q = dl*p; m = t - r
    AI = 1 - dl*r - 2*q*np.minimum((1 - r)/cI, sup)
    AP = 1 - 2*dl*m - 2*q*np.minimum(np.maximum(1 - 2*m, 0)/(cP*kap), sup)
    short = np.minimum(AI, AP).max(axis=3)
    lng = (1 - dl + (5/6 - dl)*(t - 1))[..., 0]
    R = np.maximum(short, lng).min(axis=2)            # [ND, NP]
    return dl[:,:,0,0], p[:,:,0,0], R
def Gmax(th, h, cP, cI):
    ell = 11/3 - 4*th; y = h - 2*ell; kap = 2*th - 1; best = -9
    for d in (0.5, h):
        dl, p, R = Rgrid(kap, cP, cI, ell/d)
        Cd = min(d, 2*y - d)/(cP*kap)
        G = (1+dl)/2 + 13*h/75 - y/2 - ell/2 + ell*dl*p + d*(0.25 + 0.75*R - 17/50) - (dl*p/2)*Cd
        best = max(best, G.max())
    return best
def feasible(th, cP, cI):
    for h in np.linspace(0.70, 0.86, 33):
        if Gmax(th, h, cP, cI) < th: return h
    return None
def theta(cP, cI):
    lo, hi = 0.855, 0.8765
    if feasible(hi, cP, cI) is None: return None
    for _ in range(16):
        mid = (lo+hi)/2
        if feasible(mid, cP, cI) is not None: hi = mid
        else: lo = mid
    return hi
if __name__ == "__main__":
    for cP, cI in [(6,2),(4,2),(3,2),(4/3,2),(6,1.5),(6,1),(4,1.5),(3,1),(4/3,1)]:
        print(f"cP={cP:.3f} cI={cI}: theta ~ {theta(cP,cI):.5f}", flush=True)
