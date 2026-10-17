import numpy as np, sys
from model import high_margin
def plain_P(P):  # capacity z_P = (1-2m)/P  -> plainc*kappa = P
    return P
def best_theta(label, nd=40, nq=5, **kw):
    # bisection on ell; for each ell scan x
    def feasible(ell):
        h_best=None
        for x in np.linspace(0.25, (1-ell)/2, 14):
            theta = 11/12 - ell/4; kappa = 2*theta-1
            kk = dict(kw)
            if 'P' in kk:
                P = kk.pop('P'); kk['plainc'] = P/kappa
            w,_ = high_margin(ell, x, kappa=kappa, nd=nd, nq=nq, **kk)
            if w < 0: return True, x
        return False, None
    lo, hi = 0.15, 0.2
    ok,_ = feasible(lo)
    if not ok: print(label, "infeasible at ell=0.15"); return
    for _ in range(18):
        mid=(lo+hi)/2
        ok,_=feasible(mid)
        if ok: lo=mid
        else: hi=mid
    ok,x=feasible(lo)
    print(f"{label:55s} ell_max={lo:.5f}  theta={11/12-lo/4:.6f}  x={x}")
if __name__=="__main__":
    best_theta("OpenAI inputs (P=6kappa)")
    best_theta("Cycle61 plain coefficient P=22/5", P=22/5)
    best_theta("plain P=4", P=4.0)
    best_theta("plain P=3", P=3.0)
    best_theta("plain P=2", P=2.0)
    best_theta("inverse amp alpha=3/4 (P=6kappa)", alpha=0.75)
    best_theta("inverse amp alpha=2/3 (P=6kappa)", alpha=2/3)
    best_theta("alpha=2/3, P=3", alpha=2/3, P=3.0)
