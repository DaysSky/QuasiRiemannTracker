# capc3 + long slope alpha' (LM) : long(t) = 1 - dl + (alpha' - dl)(t - 1) floored at ... (alpha'=5/6 paper)
import numpy as np, sys
ND, NP, NT, NR = 41, 11, 51, 121
def Rgrid(kap, cP, cI, sup, al):
    dl = np.linspace(0, kap, ND)[:,None,None,None]; p = np.linspace(0, .5, NP)[None,:,None,None]
    t = np.linspace(1, 1.5, NT)[None,None,:,None]; a = np.linspace(0, 1, NR)[None,None,None,:]
    lo = np.maximum(t - .5, 0); hi = np.minimum(t, 1.0); r = lo + (hi - lo)*a
    q = dl*p; m = t - r
    AI = 1 - dl*r - 2*q*np.minimum((1 - r)/cI, sup)
    AP = 1 - 2*dl*m - 2*q*np.minimum(np.maximum(1 - 2*m, 0)/(cP*kap), sup)
    short = np.minimum(AI, AP).max(axis=3)
    lng = (1 - dl + (al - dl)*(t - 1))[..., 0]
    R = np.maximum(short, lng).min(axis=2)
    return dl[:,:,0,0], p[:,:,0,0], R
def Gmax(th, h, cP, cI, al):
    ell = 11/3 - 4*th; y = h - 2*ell; kap = 2*th - 1; best = -9
    for d in (0.5, h):
        dl, p, R = Rgrid(kap, cP, cI, ell/d, al)
        Cd = min(d, 2*y - d)/(cP*kap)
        G = (1+dl)/2 + 13*h/75 - y/2 - ell/2 + ell*dl*p + d*(0.25 + 0.75*R - 17/50) - (dl*p/2)*Cd
        best = max(best, G.max())
    return best
def feasible(th, cP, cI, al):
    for h in np.linspace(0.70, 0.86, 33):
        if Gmax(th, h, cP, cI, al) < th: return h
    return None
def theta(cP, cI, al):
    lo, hi = 0.8600, 0.8765
    if feasible(hi, cP, cI, al) is None: return None
    for _ in range(16):
        mid = (lo+hi)/2
        if feasible(mid, cP, cI, al) is not None: hi = mid
        else: lo = mid
    return hi
if __name__ == "__main__":
    kap = 0.7497
    for (cPk, cI) in [(4.498, 2), (2, 2), (4.498, 1), (2, 1), (1, 1)]:
        row = []
        for al in (5/6, 0.7, 0.5, 0.3, 0.0):
            row.append(theta(cPk/kap, cI, al))
        print(f"cP*kappa={cPk:.3f} cI={cI}: " + "  ".join(f"a'={al:.2f}:{v:.5f}" for al, v in zip((5/6,0.7,0.5,0.3,0.0), row)), flush=True)
