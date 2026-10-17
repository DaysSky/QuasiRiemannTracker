# Sensitivity of theta (Cycle61 framework, capc3 grid) to the LM slope alpha' in long(t)
# and to uniform improvements of the Mobius (AI) and plain (AP) counts.
import numpy as np, sys
ND, NP, NT, NR = 41, 11, 51, 121
def Rgrid(kap, sup, al, eI, eP, cP=6, cI=2):
    dl = np.linspace(0, kap, ND)[:,None,None,None]; p = np.linspace(0, .5, NP)[None,:,None,None]
    t = np.linspace(1, 1.5, NT)[None,None,:,None]; a = np.linspace(0, 1, NR)[None,None,None,:]
    lo = np.maximum(t - .5, 0); hi = np.minimum(t, 1.0); r = lo + (hi - lo)*a
    q = dl*p; m = t - r
    AI = 1 - dl*r - 2*q*np.minimum((1 - r)/cI, sup) - eI*dl
    AP = 1 - 2*dl*m - 2*q*np.minimum(np.maximum(1 - 2*m, 0)/(cP*kap), sup) - eP*dl
    short = np.minimum(AI, AP).max(axis=3)
    lng = (1 - dl + (al - dl)*(t - 1))[..., 0]
    lng = np.maximum(lng, (1-dl)[...,0])
    R = np.maximum(short, lng).min(axis=2)
    return dl[:,:,0,0], p[:,:,0,0], R
def Gmax(th, h, al, eI, eP):
    ell = 11/3 - 4*th; y = h - 2*ell; kap = 2*th - 1; best = -9
    for d in (0.5, h):
        dl, p, R = Rgrid(kap, ell/d, al, eI, eP)
        Cd = min(d, 2*y - d)/(6*kap)
        G = (1+dl)/2 + 13*h/75 - y/2 - ell/2 + ell*dl*p + d*(0.25 + 0.75*R - 17/50) - (dl*p/2)*Cd
        best = max(best, G.max())
    return best
def feasible(th, al, eI, eP):
    for h in np.linspace(0.66, 0.86, 41):
        if Gmax(th, h, al, eI, eP) < th: return h
    return None
def theta(al=5/6, eI=0, eP=0):
    lo, hi = 0.80, 0.8765
    if feasible(hi, al, eI, eP) is None: return None
    for _ in range(16):
        mid = (lo+hi)/2
        if feasible(mid, al, eI, eP) is not None: hi = mid
        else: lo = mid
    return hi
if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else 'al'
    if mode == 'al':
        for al in (5/6, 0.8, 0.75, 0.7, 0.6, 0.5, 0.4, 0.2, 0.0):
            print(f"alpha'={al:.3f}: theta ~ {theta(al=al):.5f}", flush=True)
    elif mode == 'eI':
        for e in (0.01, 0.02, 0.05, 0.1, 0.2):
            print(f"eI={e}: theta ~ {theta(eI=e):.5f}", flush=True)
    elif mode == 'eP':
        for e in (0.01, 0.02, 0.05, 0.1, 0.2):
            print(f"eP={e}: theta ~ {theta(eP=e):.5f}", flush=True)
