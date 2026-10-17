# theta if the sextic large sieve were sharp ((U+D)^{1+eps}) for squarefree coefficients up to D = U^{a0}.
# Effects: (i) marked Mobius: slot length z = max((1-r)/2 [reflection, cost 2], a0 - r [direct sieve, cost 1])
#          (ii) plain pair: z_P = max((1-2m)/(6 kappa), a0 - 2m)
#          (iii) pure slot moments: count 1 - 2 q a0
import numpy as np
ND, NP, NT, NR = 41, 11, 51, 121
def Rgrid(kap, sup, a0):
    dl = np.linspace(0, kap, ND)[:,None,None,None]; p = np.linspace(0, .5, NP)[None,:,None,None]
    t = np.linspace(1, 1.5, NT)[None,None,:,None]; a = np.linspace(0, 1, NR)[None,None,None,:]
    lo = np.maximum(t - .5, 0); hi = np.minimum(t, 1.0); r = lo + (hi - lo)*a
    q = dl*p; m = t - r
    zI = np.maximum((1 - r)/2, a0 - r); zP = np.maximum(np.maximum(1 - 2*m, 0)/(6*kap), a0 - 2*m)
    AI = 1 - dl*r - 2*q*np.minimum(zI, sup)
    AP = 1 - 2*dl*m - 2*q*np.minimum(zP, sup)
    short = np.minimum(AI, AP).max(axis=3)
    lng = (1 - dl + (5/6 - dl)*(t - 1))[..., 0]
    R = np.maximum(short, lng).min(axis=2)
    return dl[:,:,0,0], p[:,:,0,0], np.minimum(R, 1 - 2*(dl*p)[...,0,0]*a0)
def Gmax(th, h, a0):
    ell = 11/3 - 4*th; y = h - 2*ell; kap = 2*th - 1; best = -9
    for d in (0.5, h):
        dl, p, R = Rgrid(kap, ell/d, a0)
        Cd = min(d, 2*y - d)/(6*kap)
        G = (1+dl)/2 + 13*h/75 - y/2 - ell/2 + ell*dl*p + d*(0.25 + 0.75*R - 17/50) - (dl*p/2)*Cd
        best = max(best, G.max())
    return best
def theta(a0):
    def feas(th):
        return any(Gmax(th, h, a0) < th for h in np.linspace(0.62, 0.86, 49))
    lo, hi = 0.80, 0.8765
    for _ in range(16):
        mid=(lo+hi)/2
        if feas(mid): hi=mid
        else: lo=mid
    return hi
for a0 in (0.5, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0):
    print(f"a0={a0}: theta ~ {theta(a0):.5f}", flush=True)
