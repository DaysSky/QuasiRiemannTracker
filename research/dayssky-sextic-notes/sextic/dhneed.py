import numpy as np
NR, NT, ND, NQ = 41, 61, 300, 11
def Rtab(ell, x, kappa, alpha=5/6):
    h = 1 - x + ell; P = 6*kappa; supply = ell/h
    ts = np.linspace(1.0, 1.5, NT); s = np.linspace(0, 1, NR)
    T, Sg = np.meshgrid(ts, s, indexing='ij'); r = T - Sg/2; m = Sg/2
    e = np.maximum(1.0, 1 - alpha + alpha*r)
    zM = np.where(r < 1, np.minimum(0.5*(1-r), supply), 0.0)
    zP = np.minimum((1-2*m)/P, supply)
    deltas = np.linspace(1/50, kappa, ND); xs = np.linspace(0, 0.5, NQ)
    D = deltas[:,None,None,None]; Qv = xs[None,:,None,None]*D
    c = np.broadcast_arrays(e[None,None] - D*r[None,None],
         np.where(r[None,None] < 1, 1 - D*r[None,None] - 2*Qv*zM[None,None], np.inf),
         1 - 2*D*m[None,None] - 2*Qv*zP[None,None])
    R = np.minimum.reduce(c).max(axis=3).min(axis=2)
    return deltas, xs, R
def failing(ell, x):
    theta = 11/12 - ell/4; kappa = 2*theta - 1
    h = 1 - x + ell; C0 = -x/3 - 1/12 + 13*ell/12
    deltas, xs, R = Rtab(ell, x, kappa)
    Qv = xs[None,:]*deltas[:,None]
    E = C0 + (0.5+ell)*deltas[:,None] + Qv*ell - h*(1-R)
    E_DH = C0 + (0.5+ell)*deltas[:,None] + Qv*ell - h*deltas[:,None]
    return deltas, (E.max(axis=1) >= 0), (E_DH.max(axis=1) < 0)
print("target theta | best x | delta-range where existing counts fail (DH needed there)")
for th in [0.8745, 0.8740, 0.8730, 0.8720, 0.8710, 0.8700, 0.8690, 0.8680, 0.8670]:
    ell = 4*(11/12 - th); best=None
    for x in np.linspace(max(0.2,(13*ell-1)/4+1e-4), (1-ell)/2, 40):
        d, bad, ok = failing(ell, x)
        if not np.all(ok[bad]): continue
        w = bad.sum()
        if best is None or w < best[0]: best = (w, x, d[bad], d)
    if best is None: print(f"  {th:.4f}: infeasible even with DH"); continue
    w, x, bd, d = best
    if w == 0: print(f"  {th:.4f}: x={x:.4f}  no DH needed"); continue
    segs=[]; start=bd[0]; prev=bd[0]; step=(d[1]-d[0])*1.5
    for v in bd[1:]:
        if v-prev>step: segs.append((start,prev)); start=v
        prev=v
    segs.append((start,prev))
    print(f"  {th:.4f}: x={x:.4f}  DH needed for delta in " + ", ".join(f"[{a:.3f},{b:.3f}]" for a,b in segs) + "  i.e. sigma in " + ", ".join(f"[{(1+a)/2:.3f},{(1+b)/2:.3f}]" for a,b in segs), flush=True)
