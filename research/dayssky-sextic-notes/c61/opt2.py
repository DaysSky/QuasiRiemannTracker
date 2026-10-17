import numpy as np, sys
from c61model import Gmax, low
def ok_geom(h, ell, y):
    x = 1 + ell - h
    return x > ell and y >= x and 2*x - y - ell > 0 and 4*ell - y > 0 and ell < 0.2 and h > 4*ell
def feas(th, capcost=None, nh=40):
    # returns best margin over (h, ell, y); ell from low-side equality for given M
    best = (1e9, None)
    for h in np.linspace(0.74, 0.84, nh):
        for ell in np.linspace(0.15, 0.199, 50):
            x = 1 + ell - h
            # y: low side requires 5/6 - ell/6 + M/12 + (1+3ell-2M)_+/8 + capcost*(M+ell-1)_+ <= th
            for y in np.linspace(x, min(2*x-ell, 4*ell)-1e-9, 25):
                if not ok_geom(h, ell, y): continue
                M = x + y
                if capcost is None and M + ell > 1 + 1e-12: continue
                L = low(ell, x, y) + (0 if capcost is None else capcost*max(0, M+ell-1))
                if L > th: continue
                G, a = Gmax(th, h, ell, y, nd=31, np_=31)
                if G - th < best[0]: best = (G - th, (h, ell, y, a))
    return best
def bisect(capcost=None):
    lo, hi = 0.86, 0.8752
    out=None
    for _ in range(14):
        mid = (lo+hi)/2
        f, v = feas(mid, capcost)
        if f <= 0: hi = mid; out = v
        else: lo = mid
    return hi, out
if __name__ == "__main__":
    cc = None if sys.argv[1] == "cap" else float(sys.argv[1])
    print(sys.argv[1], bisect(cc), flush=True)
