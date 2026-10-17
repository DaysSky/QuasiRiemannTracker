import numpy as np
NR, NT, ND, NQ = 41, 41, 120, 9
def high_margin(ell, x, dlo, dhi, kappa=None, alpha=5/6):
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
         1 - 2*D*m[None,None] - 2*Qv*zP[None,None]]
    c = np.broadcast_arrays(*c)
    R = np.minimum.reduce(c).max(axis=3).min(axis=2)
    dmask = (deltas >= dlo) & (deltas <= dhi)
    R = np.where(dmask[:,None], np.minimum(R, 1 - deltas[:,None]), R)   # DH available on [dlo,dhi]
    E = C0 + (0.5+ell)*D[:,:,0,0] + Qv[:,:,0,0]*ell - h*(1-R)
    return E.max()
def feasible(ell, dlo, dhi):
    for x in np.linspace(0.25, (1-ell)/2, 30):
        if high_margin(ell, x, dlo, dhi) < 0: return True
    return False
targets = [0.8740, 0.8730, 0.8720, 0.8710, 0.8700, 0.8690, 0.8680]
print("theta target | minimal DH window [dlo,dhi] needed (search over a grid)")
for th in targets:
    ell = 4*(11/12 - th)
    best=None
    for dlo in np.arange(0.0, 0.75, 0.025):
        for dhi in np.arange(dlo+0.05, 0.80, 0.025):
            if feasible(ell, dlo, dhi):
                w = dhi - dlo
                if best is None or w < best[0]: best = (w, dlo, dhi)
                break
    print(f"  {th:.4f} (ell={ell:.4f}): ", "none needed" if feasible(ell, 1, 0) else (f"[{best[1]:.3f},{best[2]:.3f}] width {best[0]:.3f}" if best else "not reachable even with DH on a window"), flush=True)
