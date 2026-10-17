# Generalized 7/8-architecture model with a free theta and a partially removed exceptional Gram term.
# low exponent  = x/2 + (1-alpha)*b/12 + pen,  b = 1 - ell - 2x,  pen = max(0, (2b-(y-ell)) - (1-alpha) b/6)/2 [nonexceptional Gram term]
# principal C(th) = th - 5/6 + x/3 + ell/6 ; need C(th) > low
# rows: E = ell + (1/2+ell) d + q ell - h (1-R) - C(th) < 0  for all bins d in (0,kappa], amplitudes q<=d/2
import numpy as np
NR, NT, ND, NQ = 41, 41, 60, 7
def Rtab(ell, x, kappa, alpha5=5/6, plainc=6.0, invcap=0.5):
    h = 1 - x + ell; P = plainc*kappa; supply = ell/h
    ts = np.linspace(1.0, 1.5, NT); s = np.linspace(0, 1, NR)
    T, Sg = np.meshgrid(ts, s, indexing='ij'); r = T - Sg/2; m = Sg/2
    e = np.maximum(1.0, 1 - alpha5 + alpha5*r)
    zM = np.where(r < 1, np.minimum(invcap*(1-r), supply), 0.0)
    zP = np.minimum((1-2*m)/P, supply)
    deltas = np.linspace(1/50, kappa, ND); xs = np.linspace(0, 0.5, NQ)
    D = deltas[:,None,None,None]; X = xs[None,:,None,None]; Qv = X*D
    c1 = e[None,None] - D*r[None,None]
    c2 = np.where(r[None,None] < 1, 1 - D*r[None,None] - 2*Qv*zM[None,None], np.inf)
    c3 = 1 - 2*D*m[None,None] - 2*Qv*zP[None,None]
    val = np.minimum(np.minimum(c1, c2), c3)
    R = val.max(axis=3).min(axis=2)
    return deltas, xs, R
def feasible(th, ell, x, alpha):
    if not (0 < ell <= 0.25 and 0.05 < x): return False
    b = 1 - ell - 2*x
    if b < 0: return False
    y = x + b
    pen = max(0.0, (2*b - (y - ell)) - (1-alpha)*b/6)/2
    low = x/2 + (1-alpha)*b/12 + pen
    C = th - 5/6 + x/3 + ell/6
    if C <= low + 1e-9: return False
    kappa = 2*th - 1
    h = 1 - x + ell
    d, xs, R = Rtab(ell, x, kappa)
    Qv = xs[None,:]*d[:,None]
    E = ell + (0.5+ell)*d[:,None] + Qv*ell - h*(1-R) - C
    if E.max() >= 0: return False
    # generic rows (delta=0): ell - C < 0 ; completed-row excess needs ell<=1/5 roughly
    if ell - C >= 0: return False
    return True
def best(alpha, grid_ell=np.linspace(0.08,0.2,25), grid_x=np.linspace(0.2,0.45,26)):
    best=(1,None)
    for ell in grid_ell:
        for x in grid_x:
            lo, hi = 0.80, 0.92
            if not feasible(hi, ell, x, alpha): continue
            for _ in range(22):
                mid=(lo+hi)/2
                if feasible(mid, ell, x, alpha): hi=mid
                else: lo=mid
            if hi < best[0]: best=(hi,(round(ell,4),round(x,4)))
    return best
if __name__ == "__main__":
    import sys
    for alpha in [float(a) for a in sys.argv[1:]] or [0.0, 0.5, 1.0]:
        print("alpha", alpha, best(alpha), flush=True)
