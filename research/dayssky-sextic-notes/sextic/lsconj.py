# theta under a sharp sextic large sieve (U+N)^{1+eps} for general coefficients up to N = U^{a0}:
#  (a) pure slot moments: count <= 1 - 2*q*a0  (q = delta*p)
#  (b) slot costs cI, cP reduced
import numpy as np, sys
import lmsens as L
ND, NP, NT, NR = 41, 11, 51, 121
def Rgrid(kap, sup, a0, cP, cI):
    dl = np.linspace(0, kap, ND)[:,None,None,None]; p = np.linspace(0, .5, NP)[None,:,None,None]
    t = np.linspace(1, 1.5, NT)[None,None,:,None]; a = np.linspace(0, 1, NR)[None,None,None,:]
    lo = np.maximum(t - .5, 0); hi = np.minimum(t, 1.0); r = lo + (hi - lo)*a
    q = dl*p; m = t - r
    AI = 1 - dl*r - 2*q*np.minimum((1 - r)/cI, sup)
    AP = 1 - 2*dl*m - 2*q*np.minimum(np.maximum(1 - 2*m, 0)/(cP*kap), sup)
    short = np.minimum(AI, AP).max(axis=3)
    lng = (1 - dl + (5/6 - dl)*(t - 1))[..., 0]
    R = np.maximum(short, lng).min(axis=2)
    Rslot = 1 - 2*(dl*p)[...,0,0]*a0
    return dl[:,:,0,0], p[:,:,0,0], np.minimum(R, Rslot)
def Gmax(th, h, a0, cP, cI):
    ell = 11/3 - 4*th; y = h - 2*ell; kap = 2*th - 1; best = -9
    for d in (0.5, h):
        dl, p, R = Rgrid(kap, ell/d, a0, cP, cI)
        Cd = min(d, 2*y - d)/(cP*kap)
        G = (1+dl)/2 + 13*h/75 - y/2 - ell/2 + ell*dl*p + d*(0.25 + 0.75*R - 17/50) - (dl*p/2)*Cd
        best = max(best, G.max())
    return best
def theta(a0, cP=6, cI=2):
    def feas(th):
        for h in np.linspace(0.62, 0.86, 49):
            if Gmax(th, h, a0, cP, cI) < th: return True
        return False
    lo, hi = 0.80, 0.8765
    if not feas(hi): return None
    for _ in range(16):
        mid=(lo+hi)/2
        if feas(mid): hi=mid
        else: lo=mid
    return hi
if __name__=="__main__":
    for a0,cP,cI in [(0,6,2),(0.5,6,2),(0.75,6,2),(0.9,6,2),(1,6,2),(1,6,1),(1,4/3,1),(1,1,1)]:
        print(f"a0={a0} cP={cP:.3f} cI={cI}: theta ~ {theta(a0,cP,cI):.5f}", flush=True)
