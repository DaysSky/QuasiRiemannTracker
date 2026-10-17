import numpy as np
NR, NT, ND, NQ = 41, 41, 80, 9
def high_margin(ell, x, kappa=None, alpha=5/6, mmax=0.5, loss=0.0):
    h = 1 - x + ell; C0 = -x/3 - 1/12 + 13*ell/12
    theta = 11/12 - ell/4
    if kappa is None: kappa = 2*theta-1
    P = 6*kappa; supply = ell/h
    ts = np.linspace(1.0, 1.5, NT); s = np.linspace(0, 1, NR)
    T, Sg = np.meshgrid(ts, s, indexing='ij'); r = T - Sg/2; m = Sg/2
    e = np.maximum(1.0, 1 - alpha + alpha*r)
    zM = np.where(r < 1, np.minimum(0.5*(1-r), supply), 0.0)
    zP = np.minimum((1-2*m)/P, supply)
    deltas = np.linspace(1/50, kappa, ND); xs = np.linspace(0, 0.5, NQ)
    D = deltas[:,None,None,None]; Qv = xs[None,:,None,None]*D
    c = [e[None,None] - D*r[None,None],
         np.where(r[None,None] < 1, 1 - D*r[None,None] - 2*Qv*zM[None,None], np.inf),
         1 - 2*D*m[None,None] - 2*Qv*zP[None,None],
         np.where(((r+m) <= 1+1e-9)[None,None] & (m <= mmax+1e-9)[None,None], 1 + loss - D*(r+m)[None,None], np.inf)]
    c = np.broadcast_arrays(*c)
    R = np.minimum.reduce(c).max(axis=3).min(axis=2)
    E = C0 + (0.5+ell)*D[:,:,0,0] + Qv[:,:,0,0]*ell - h*(1-R)
    return E.max()
def best_theta(label, **kw):
    def feas(ell):
        for x in np.linspace(0.2, (1-ell)/2, 30):
            if high_margin(ell, x, **kw) < 0: return x
        return None
    lo, hi = 0.10, 0.2001
    for _ in range(22):
        mid=(lo+hi)/2
        if feas(mid) is not None: lo=mid
        else: hi=mid
    print(f"{label:46s} ell_max={lo:.5f} theta={11/12-lo/4:.6f}", flush=True)
# barrier-free range of the all-row (raw) MM is m <= 1/3 (principal rows a^6); evaluate the model there
import sys
for m0 in [float(a) for a in sys.argv[1:]] or (0.25, 1/3, 0.35, 0.375, 0.4, 0.5):
    best_theta(f"MM only for m <= {m0:.4f}", mmax=m0)
