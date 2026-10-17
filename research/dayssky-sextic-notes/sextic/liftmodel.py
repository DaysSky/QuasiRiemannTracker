import numpy as np
from p6model import NR, NT, ND, NQ
def high_margin(ell, x, ks=(), mob4lift=False, alpha=5/6, plainc=6.0, invcap=0.5, tmax=1.5, liftslope=5/6):
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
    rr = r[None,None]; mm = m[None,None]
    c1 = e[None,None] - D*rr
    c2 = np.where(rr < 1, 1 - D*rr - 2*Qv*zM[None,None], np.inf)
    c3 = 1 - 2*D*mm - 2*Qv*zP[None,None]
    val = np.minimum(np.minimum(c1, c2), c3)
    for k in ks:   # plain 2k-th moment: Lindelof at critical k*m<=1, sixth-power lift beyond
        val = np.minimum(val, 1 + liftslope*np.maximum(0, k*mm - 1) - k*D*mm)
    if mob4lift:   # Mobius 4th: critical 2r<=1, lift beyond
        val = np.minimum(val, np.where(rr <= 1, 1 + liftslope*np.maximum(0, 2*rr - 1) - 2*D*rr, np.inf))
    R = val.max(axis=3).min(axis=2)
    E = C0 + (0.5+ell)*D[:,:,0,0] + Qv[:,:,0,0]*ell - h*(1-R)
    return E.max()
def best_theta(**kw):
    def feas(ell):
        for x in np.linspace(0.25, (1-ell)/2, 30):
            if high_margin(ell, x, **kw) < 0: return x
        return None
    lo, hi = 0.12, 0.2
    if feas(lo) is None: return None
    for _ in range(20):
        mid=(lo+hi)/2
        if feas(mid) is not None: lo=mid
        else: hi=mid
    return round(11/12-lo/4, 6), round(lo,5), feas(lo)
if __name__=="__main__":
    print("baseline", best_theta(), flush=True)
    print("plain 6th crit+lift", best_theta(ks=(3,)), flush=True)
    print("plain 6th,8th crit+lift", best_theta(ks=(3,4)), flush=True)
    print("plain 2k k=3..8 crit+lift", best_theta(ks=(3,4,5,6,7,8)), flush=True)
    print("Mobius 4th crit+lift", best_theta(mob4lift=True), flush=True)
    print("all", best_theta(ks=(3,4,5,6), mob4lift=True), flush=True)
