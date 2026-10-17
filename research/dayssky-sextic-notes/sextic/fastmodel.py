import numpy as np
# vectorised version of model.high_margin
NR, NT, ND, NQ = 41, 41, 80, 9
def R_table(deltas, qs, supply, kappa, alpha=5/6, P=None, plainc=6.0, invcap=0.5, e_override=None):
    """returns R[d_i, q_j] = min_t max_{(r,m)} min(counts)"""
    if P is None: P = plainc*kappa
    ts = np.linspace(1.0, 1.5, NT)
    s = np.linspace(0, 1, NR)                       # r = t - s/2, m = s/2  (r+m=t, m in [0,1/2])
    T, Sg = np.meshgrid(ts, s, indexing='ij')       # [NT, NR]
    r = T - Sg/2; m = Sg/2
    e = np.maximum(1.0, 1 - alpha + alpha*r)
    zM = np.where(r < 1, np.minimum(invcap*(1-r), supply), 0.0)
    zP = np.minimum((1-2*m)/P, supply)
    D = deltas[:, None, None, None]; Q = qs[None, :, None, None]   # q given as absolute values per delta? handle below
    return r, m, e, zM, zP
def high_margin(ell, x, kappa=None, alpha=5/6, P=None, plainc=6.0, invcap=0.5):
    h = 1 - x + ell; C0 = -x/3 - 1/12 + 13*ell/12
    theta = 11/12 - ell/4
    if kappa is None: kappa = 2*theta-1
    if P is None: P = plainc*kappa
    supply = ell/h
    ts = np.linspace(1.0, 1.5, NT); s = np.linspace(0, 1, NR)
    T, Sg = np.meshgrid(ts, s, indexing='ij'); r = T - Sg/2; m = Sg/2
    e = np.maximum(1.0, 1 - alpha + alpha*r)
    zM = np.where(r < 1, np.minimum(invcap*(1-r), supply), 0.0)
    zP = np.minimum((1-2*m)/P, supply)
    deltas = np.linspace(1/50, kappa, ND)
    xs = np.linspace(0, 0.5, NQ)                    # q = xs*delta
    D = deltas[:,None,None,None]; X = xs[None,:,None,None]; Qv = X*D
    c1 = e[None,None] - D*r[None,None]
    c2 = np.where(r[None,None] < 1, 1 - D*r[None,None] - 2*Qv*zM[None,None], np.inf)
    c3 = 1 - 2*D*m[None,None] - 2*Qv*zP[None,None]
    val = np.minimum(np.minimum(c1, c2), c3)        # [ND,NQ,NT,NR]
    R = val.max(axis=3).min(axis=2)                 # adversary over r, we choose t
    E = C0 + (0.5+ell)*D[:,:,0,0] + Qv[:,:,0,0]*ell - h*(1-R)
    i = np.unravel_index(np.argmax(E), E.shape)
    return E.max(), (deltas[i[0]], xs[i[1]], R[i])
def best_theta(label, **kw):
    def feas(ell):
        for x in np.linspace(0.2, (1-ell)/2, 40):
            if high_margin(ell, x, **kw)[0] < 0: return x
        return None
    lo, hi = 0.10, 0.20
    if feas(lo) is None: print(label, "infeasible at", lo); return None
    for _ in range(30):
        mid=(lo+hi)/2
        if feas(mid) is not None: lo=mid
        else: hi=mid
    x=feas(lo)
    print(f"{label:52s} ell_max={lo:.6f} theta={11/12-lo/4:.6f} x={x:.4f}")
    return lo
if __name__=="__main__":
    print(high_margin(1/6, 17/48, kappa=3/4))
    best_theta("OpenAI inputs (P=6kappa)")
    best_theta("Cycle61-type plain P=22/5", P=22/5)
    for P in (4.0, 3.5, 3.0, 2.5, 2.0, 1.0):
        best_theta(f"plain P={P}", P=P)
    for a in (0.8, 0.75, 2/3, 0.5):
        best_theta(f"inverse amplification alpha={a:.3f}", alpha=a)
    for c in (0.6, 0.75, 1.0):
        best_theta(f"inverse capacity (1-r)*{c}", invcap=c)
    best_theta("all: alpha=1/2, P=1, invcap=1", alpha=0.5, P=1.0, invcap=1.0)
