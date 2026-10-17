# Gao's (B, Lam) with a wider Theorem 2.2 window; old-branch weight stays 6*kappa.
import sys; sys.path.insert(0, '.')
zlo, zhi, nlo, nhi = map(float, sys.argv[1:5])
import fastRw, thoptw_core as T
fastRw.CPW = T.CPW = 6*0.749672
fastRw.B = T.B = 0.987047; fastRw.Lam = T.Lam = 3.9586
fastRw.WIN.update(ZLO=zlo, ZHI=zhi, NLO=nlo, NHI=nhi)
def parts(th, h):
    ell = 11/3 - 4*th; y = h - 2*ell
    return (th - 1e-6 - T.Gmax(th, h, ell, y, h, 'new'), th - 1e-5 - T.Gmax(th, h, ell, y, h - 1e-3, 'old'),
            th - 1e-3 - T.Gmax(th, h, ell, y, 0.5, 'old'))
from scipy.optimize import minimize_scalar
lo, hi = 0.870, 0.8752; hb = None
for _ in range(18):
    mid = (lo + hi)/2
    res = minimize_scalar(lambda h: -min(parts(mid, h)), bounds=(0.76, 0.82), method='bounded', options={'xatol':1e-8})
    if -res.fun >= 0: hi = mid; hb = res.x
    else: lo = mid
print(f"window z[{zlo},{zhi}] n[{nlo},{nhi}]: theta_min {hi:.7f} (h {hb:.5f})", flush=True)
