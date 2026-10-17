import numpy as np
import fastmodel as fm
NR, NT, ND, NQ = 41, 41, 80, 9
def high_margin(ell, x, kappa=None, alpha=5/6, P=None, invcap=0.5, mixed=None, mixedcap=0.0):
    """mixed: None, or total-length limit L for a mixed moment sum|M_r S_m|^2 << U  when r+m<=L,
       giving count 1 - delta(r+m) - 2 q min(mixedcap*(L-r-m), supply)."""
    h = 1 - x + ell; C0 = -x/3 - 1/12 + 13*ell/12
    theta = 11/12 - ell/4
    if kappa is None: kappa = 2*theta-1
    if P is None: P = 6*kappa
    supply = ell/h
    ts = np.linspace(1.0, 1.5, NT); s = np.linspace(0, 1, NR)
    T, Sg = np.meshgrid(ts, s, indexing='ij'); r = T - Sg/2; m = Sg/2
    e = np.maximum(1.0, 1 - alpha + alpha*r)
    zM = np.where(r < 1, np.minimum(invcap*(1-r), supply), 0.0)
    zP = np.minimum((1-2*m)/P, supply)
    deltas = np.linspace(1/50, kappa, ND); xs = np.linspace(0, 0.5, NQ)
    D = deltas[:,None,None,None]; Qv = xs[None,:,None,None]*D
    c = [e[None,None] - D*r[None,None],
         np.where(r[None,None] < 1, 1 - D*r[None,None] - 2*Qv*zM[None,None], np.inf),
         1 - 2*D*m[None,None] - 2*Qv*zP[None,None]]
    if mixed is not None:
        tot = (r+m)[None,None]
        zX = np.minimum(mixedcap*np.maximum(mixed - tot, 0), supply)
        c.append(np.where(tot <= mixed, 1 - D*tot - 2*Qv*zX, np.inf))
    c = np.broadcast_arrays(*c)
    val = np.minimum.reduce(c)
    R = val.max(axis=3).min(axis=2)
    E = C0 + (0.5+ell)*D[:,:,0,0] + Qv[:,:,0,0]*ell - h*(1-R)
    return E.max()
def best_theta(label, **kw):
    def feas(ell):
        for x in np.linspace(0.2, (1-ell)/2, 40):
            if high_margin(ell, x, **kw) < 0: return x
        return None
    lo, hi = 0.10, 0.2
    for _ in range(26):
        mid=(lo+hi)/2
        if feas(mid) is not None: lo=mid
        else: hi=mid
    print(f"{label:58s} ell_max={lo:.5f} theta={11/12-lo/4:.6f}", flush=True)
best_theta("baseline (OpenAI inputs)")
best_theta("mixed moment |M_r S_m|^2, r+m<=1, no slots", mixed=1.0)
best_theta("mixed moment r+m<=1, slot capacity (1-r-m)/2", mixed=1.0, mixedcap=0.5)
best_theta("mixed moment r+m<=0.9, no slots", mixed=0.9)
best_theta("mixed moment r+m<=0.8, no slots", mixed=0.8)
best_theta("mixed moment r+m<=0.75, no slots", mixed=0.75)
