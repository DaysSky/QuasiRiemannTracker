# Cycle61 row count with Moebius 4th moment M4 valid for r <= r4 on every family, plus sixth-power lift beyond r4:
#   sum_{u~U} |M_u(U^r)|^4 << U^{e4(r)},  e4 = 1 (r<=r4),  1/6 + 5r/(6 r4) (r>r4)   [H=U^{r/r4}, P=(H/U)^{1/6}]
import numpy as np, sys
Lam = 19312/4787; B = 14200/14361 - 1e-5
def e4(r, r4): return np.where(r <= r4, 1.0, 1/6 + 5*r/(6*r4))
def Rcount(dl, p, kap, r4, slot=True, nt=101, nr=151):
    q = dl*p; best = 9
    for t in np.linspace(1, 1.5, nt):
        rs = np.linspace(max(t-0.5, 0), min(t, 1), nr); m = t - rs
        AI = 1 - q - (dl - q)*rs
        AP = 1 - q/(3*kap) - 2*(dl - q/(3*kap))*m
        A4 = (e4(rs, r4) - q - (2*dl - q)*rs) if slot else (e4(rs, r4) - 2*dl*rs)
        short = np.minimum(np.minimum(AI, AP), A4).max()
        rl = np.linspace(1, t, 21) if t > 1 else np.array([1.0])
        lng = max(1 - dl + (5/6 - dl)*(t - 1), (np.minimum(1/6 + (5/6 - dl)*rl, e4(rl, r4) - 2*dl*rl)).max() if t>1 else 0)
        best = min(best, max(short, lng))
    return best
def Gmax(th, h, ell, y, r4, slot, nd=21, np_=11):
    kap = 2*th - 1; best = -9
    for dl in np.linspace(0, kap, nd):
        for p in np.linspace(0, 0.5, np_):
            R = Rcount(dl, p, kap, r4, slot)
            for d in (h - 1e-3, h):
                Cd = (2*y - (2-B)*d)/Lam
                G = (1+dl)/2 + 13*h/75 - y/2 - ell/2 + ell*dl*p + d*(0.25 + 0.75*R - 17/50) - (dl*p/2)*Cd
                best = max(best, G)
    return best
def feas(th, r4, slot):
    ell = 4*(11/12 - th); best = 9
    for h in np.linspace(0.76, 0.84, 17):
        y = h - 2*ell; x = 1 + ell - h
        if not (x > ell and y >= x and 4*ell > y and ell < 0.2): continue
        best = min(best, Gmax(th, h, ell, y, r4, slot) - th)
    return best
r4 = float(sys.argv[1]); slot = sys.argv[2] == 'slot'
lo, hi = 0.866, 0.8752
for _ in range(12):
    mid = (lo+hi)/2
    if feas(mid, r4, slot) <= 0: hi = mid
    else: lo = mid
print("r4", r4, "slot" if slot else "noslot", round(hi, 6), flush=True)
