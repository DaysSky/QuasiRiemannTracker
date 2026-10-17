import numpy as np, sys
from c61m4 import Gmax
def feas(th, r4):
    ell = 4*(11/12 - th)
    best = 9
    for h in np.linspace(0.76, 0.84, 17):
        y = h - 2*ell; x = 1 + ell - h
        if not (x > ell and y >= x and 4*ell > y and ell < 0.2): continue
        g = Gmax(th, h, ell, y, r4=r4, nd=21, np_=11) - th
        best = min(best, g)
    return best
def bisect(r4):
    lo, hi = 0.866, 0.8752
    for _ in range(12):
        mid = (lo+hi)/2
        if feas(mid, r4) <= 0: hi = mid
        else: lo = mid
    return hi
r4 = sys.argv[1]
r4 = r4 if r4 in ('none','formula') else float(r4)
print(sys.argv[1], round(bisect(r4), 6), flush=True)
