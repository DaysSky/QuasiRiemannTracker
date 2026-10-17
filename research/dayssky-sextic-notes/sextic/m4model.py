import numpy as np
NR, NT, ND, NQ = 61, 61, 80, 9
def high_margin(ell, x, r4=0.0, lam4=0.0, alpha=5/6, plainc=6.0, invcap=0.5, tmax=1.5):
    h = 1 - x + ell; C0 = -x/3 - 1/12 + 13*ell/12
    theta = 11/12 - ell/4; kappa = 2*theta-1; P = plainc*kappa
    supply = ell/h
    ts = np.linspace(1.0, tmax, NT); s = np.linspace(0, 1, NR)
    T, Sg = np.meshgrid(ts, s, indexing='ij'); r = T - Sg/2; m = Sg/2
    e = np.maximum(1.0, 1 - alpha + alpha*r)
    zM = np.where(r < 1, np.minimum(invcap*(1-r), supply), 0.0)
    zP = np.minimum((1-2*m)/P, supply)
    deltas = np.linspace(1/50, kappa, ND); xs = np.linspace(0, 0.5, NQ)
    D = deltas[:,None,None,None]; X = xs[None,:,None,None]; Qv = X*D
    c1 = e[None,None] - D*r[None,None]
    c2 = np.where(r[None,None] < 1, 1 - D*r[None,None] - 2*Qv*zM[None,None], np.inf)
    c3 = 1 - 2*D*m[None,None] - 2*Qv*zP[None,None]
    c4 = np.where(r[None,None] <= r4+1e-12, 1 + lam4 - 2*D*r[None,None], np.inf)   # Mobius 4th moment witness
    val = np.minimum(np.minimum(np.minimum(c1, c2), c3), c4)
    R = val.max(axis=3).min(axis=2)
    E = C0 + (0.5+ell)*D[:,:,0,0] + Qv[:,:,0,0]*ell - h*(1-R)
    return E.max()
def best_theta(**kw):
    def feas(ell):
        for x in np.linspace(0.2, (1-ell)/2, 50):
            if high_margin(ell, x, **kw) < 0: return x
        return None
    lo, hi = 0.10, 0.2
    if feas(lo) is None: return None
    for _ in range(26):
        mid=(lo+hi)/2
        if feas(mid) is not None: lo=mid
        else: hi=mid
    return 11/12-lo/4, lo, feas(lo)
if __name__=="__main__":
    import sys
    print("baseline", best_theta())
    for r4 in (0.5, 0.6, 2/3, 0.7, 0.75, 0.8, 0.9, 1.0):
        print(f"M4 r<= {r4:.3f} lam=0   ", best_theta(r4=r4), flush=True)
    for lam in (0.01, 0.03, 0.05):
        print(f"M4 r<= 1     lam={lam}", best_theta(r4=1.0, lam4=lam), flush=True)
