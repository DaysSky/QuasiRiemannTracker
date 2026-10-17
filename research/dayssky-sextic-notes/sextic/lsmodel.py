import numpy as np
from p6model import NR, NT, ND, NQ
def high_margin(ell, x, p6ls=False, k=3, m6=0.5, m4ls=False, alpha=5/6, plainc=6.0, invcap=0.5, tmax=1.5):
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
    if p6ls:   # plain 2k-th moment in large-sieve shape: sum |S_m|^{2k} << U^{max(1,k m)}
        val = np.minimum(val, np.where(mm <= m6+1e-12, np.maximum(1, k*mm) - k*D*mm, np.inf))
    if m4ls:   # Mobius 4th moment in large-sieve shape
        val = np.minimum(val, np.where(rr <= 1, np.maximum(1, 2*rr) - 2*D*rr, np.inf))
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
    return round(11/12-lo/4, 6), round(lo,5)
if __name__=="__main__":
    print("baseline", best_theta(), flush=True)
    print("plain 6th LS-shape m<=1/2", best_theta(p6ls=True), flush=True)
    print("plain 8th LS-shape m<=1/2", best_theta(p6ls=True, k=4), flush=True)
    print("plain 6th+8th LS", best_theta(p6ls=True, k=3), flush=True)
    print("Mobius 4th LS-shape", best_theta(m4ls=True), flush=True)
    print("both", best_theta(p6ls=True, m4ls=True), flush=True)
