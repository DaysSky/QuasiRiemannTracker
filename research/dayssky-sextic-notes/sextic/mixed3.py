import numpy as np
# adversary chooses r in [t-1/2, t], m in [max(0,t-r), 1/2]  (so r+m in [t, t+1/2])
NT, NRr, NM, ND, NQ = 31, 26, 26, 80, 9
def high_margin(ell, x, kappa=None, alpha=5/6, mm=None, mmloss=0.0):
    """mm=None: no mixed moment; mm='crit': sum|M_r S_m|^2 << U^{1+loss} only for r+m<=1;
       mm='all': << U^{max(1,r+m)?}  -> we model 'all' as << U^{1+loss} for every r<=1, m<=1/2 (Lindelof-on-average type)."""
    h = 1 - x + ell; C0 = -x/3 - 1/12 + 13*ell/12
    theta = 11/12 - ell/4
    if kappa is None: kappa = 2*theta-1
    P = 6*kappa; supply = ell/h
    ts = np.linspace(1.0, 1.5, NT)
    a = np.linspace(0, 1, NRr); b = np.linspace(0, 1, NM)
    T, A, B = np.meshgrid(ts, a, b, indexing='ij')
    r = T - 0.5 + 0.5*A                      # r in [t-1/2, t]
    mlo = np.maximum(0.0, T - r); m = mlo + (0.5 - mlo)*B
    valid = (m <= 0.5 + 1e-12)
    e = np.maximum(1.0, 1 - alpha + alpha*r)
    zM = np.where(r < 1, np.minimum(0.5*(1-r), supply), 0.0)
    zP = np.minimum((1-2*m)/P, supply)
    deltas = np.linspace(1/50, kappa, ND); xs = np.linspace(0, 0.5, NQ)
    D = deltas[:,None,None,None,None]; Qv = xs[None,:,None,None,None]*D
    R_ = r[None,None]; M_ = m[None,None]
    cands = [e[None,None] - D*R_,
             np.where(R_ < 1, 1 - D*R_ - 2*Qv*zM[None,None], np.inf),
             1 - 2*D*M_ - 2*Qv*zP[None,None]]
    if mm == 'crit':
        cands.append(np.where((R_+M_) <= 1 + 1e-9, 1 + mmloss - D*(R_+M_), np.inf))
    elif mm == 'all':
        cands.append(np.where(R_ <= 1 + 1e-9, 1 + mmloss - D*(R_+M_), np.inf))
    cands = np.broadcast_arrays(*cands)
    val = np.minimum.reduce(cands)
    val = np.where(valid[None,None], val, -np.inf)
    Rv = val.max(axis=(3,4)).min(axis=2)       # adversary over (r,m), we choose t
    E = C0 + (0.5+ell)*D[:,:,0,0,0] + Qv[:,:,0,0,0]*ell - h*(1-Rv)
    return E.max()
def best_theta(label, **kw):
    def feas(ell):
        for x in np.linspace(0.25, (1-ell)/2, 25):
            if high_margin(ell, x, **kw) < 0: return x
        return None
    lo, hi = 0.12, 0.2001
    for _ in range(20):
        mid=(lo+hi)/2
        if feas(mid) is not None: lo=mid
        else: hi=mid
    print(f"{label:60s} ell_max={lo:.5f} theta={11/12-lo/4:.6f}", flush=True)
best_theta("baseline (r+m free)")
best_theta("MM critical only (r+m<=1)", mm='crit')
best_theta("MM' all r<=1, m<=1/2 (Lindelof-on-average product)", mm='all')
best_theta("MM' all, loss U^0.02", mm='all', mmloss=0.02)
