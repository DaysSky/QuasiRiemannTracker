# Effect of widening Cycle61 Theorem 2.2's validity window (plain capacity z = (B - 2m)/Lambda) on theta.
import numpy as np, sys
ND, NP, NT, NR = 41, 11, 51, 161
B_ = 14200/14361 - 1e-5; Lam = 19312/4787
def Rgrid(kap, sup, mlo, mhi, zlo, zhi):
    dl = np.linspace(0, kap, ND)[:,None,None,None]; p = np.linspace(0, .5, NP)[None,:,None,None]
    t = np.linspace(1, 1.5, NT)[None,None,:,None]; a = np.linspace(0, 1, NR)[None,None,None,:]
    lo = np.maximum(t - .5, 0); hi = np.minimum(t, 1.0); r = lo + (hi - lo)*a
    q = dl*p; m = t - r
    zold = np.maximum(1 - 2*m, 0)/(6*kap)
    znew = (B_ - 2*m)/Lam
    ok = (m >= mlo) & (m <= mhi) & (znew >= zlo) & (znew <= zhi)
    zP = np.where(ok, np.maximum(zold, znew), zold)
    AI = 1 - dl*r - 2*q*np.minimum((1 - r)/2, sup)
    AP = 1 - 2*dl*m - 2*q*np.minimum(zP, sup)
    short = np.minimum(AI, AP).max(axis=3)
    lng = (1 - dl + (5/6 - dl)*(t - 1))[..., 0]
    R = np.maximum(short, lng).min(axis=2)
    return dl[:,:,0,0], p[:,:,0,0], R
def Gmax(th, h, win):
    ell = 11/3 - 4*th; y = h - 2*ell; kap = 2*th - 1; best = -9
    for d in (0.5, h):
        dl, p, R = Rgrid(kap, ell/d, *win)
        Cd = min(d, 2*y - d)/(6*kap) if d == 0.5 else max(min(d, 2*y - d)/(6*kap), (2*y - (2 - B_)*d)/Lam)
        G = (1+dl)/2 + 13*h/75 - y/2 - ell/2 + ell*dl*p + d*(0.25 + 0.75*R - 17/50) - (dl*p/2)*Cd
        best = max(best, G.max())
    return best
def theta(win):
    def feas(th): return any(Gmax(th, h, win) < th for h in np.linspace(0.70, 0.86, 81))
    lo, hi = 0.86, 0.8765
    for _ in range(18):
        mid = (lo+hi)/2
        if feas(mid): hi = mid
        else: lo = mid
    return hi
wins = {'none': (1, 0, 1, 0), 'cycle61': ((B_-.04*Lam)/2, min(.44, (B_-.034*Lam)/2), .034, .04),
        'n[.40,.44] z[.03,.045]': (.40, .44, .03, .045), 'n[.38,.45] z[.024,.06]': (.38, .45, .024, .06),
        'n[.33,.5] z[0,.1]': (.33, .5, 0, .1)}
for k, w in wins.items():
    print(f"{k:28s} theta ~ {theta(w):.6f}", flush=True)
