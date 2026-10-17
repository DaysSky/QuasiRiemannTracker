# Record pipeline with a partial LM saving: long(t) -> long(t) - eta on t in [s1, s2]  (MODEL)
import sys, numpy as np, fastRw
from scipy.optimize import brentq, minimize_scalar
eta, s1, s2 = map(float, sys.argv[1:4])
import os
Brad = float(os.environ.get('GB', 14200/14361 + 7/10000 - 1e-5))
fastRw.B = Brad
if 'GL' in os.environ: fastRw.Lam = float(os.environ['GL'])
def R(dl, p, kap):
    sh = lambda t: fastRw.short(t, dl, p, kap)
    base = lambda t: 1 - dl + (5/6 - dl)*(t - 1)
    pieces = [(1.0, s1, 0.0), (s1, s2, eta), (s2, 1.5, 0.0)] if eta > 0 else [(1.0, 1.5, 0.0)]
    best = 9.0
    for a, b, e in pieces:
        if b <= a: continue
        g = lambda t: sh(t) - (base(t) - e)
        ga, gb = g(a), g(b)
        if ga > 0 and gb < 0:
            tc = brentq(g, a, b, xtol=1e-13); v = base(tc) - e
        else:
            v = min(max(sh(a), base(a) - e), max(sh(b), base(b) - e))
        best = min(best, v)
    return best
fastRw.R = R
import thoptw_core as T
T.B = Brad; T.R = R
if 'GL' in os.environ: T.Lam = float(os.environ['GL'])
def parts(th, h):
    ell = 11/3 - 4*th; y = h - 2*ell; x = 1 + ell - h
    if not (0 < ell < x < y < 2*x and y < 0.5 < h < 2*y and 1/6 < ell < 1/5 and h > 4*ell and 4*ell - y > 0): return (-1,)
    return (th - 1e-6 - T.Gmax(th, h, ell, y, h, 'new'), th - 1e-5 - T.Gmax(th, h, ell, y, h - 1e-3, 'old'),
            th - 1e-3 - T.Gmax(th, h, ell, y, 0.5, 'old'))
lo, hi = float(sys.argv[4]) if len(sys.argv)>4 else 0.8740, 0.87490; hb = None
for _ in range(14):
    mid = (lo + hi)/2
    res = minimize_scalar(lambda h: -min(parts(mid, h)), bounds=(0.74, 0.80), method='bounded', options={'xatol':1e-7})
    if -res.fun >= 0: hi = mid; hb = res.x
    else: lo = mid
print(f"eta={eta} window=[{s1},{s2}] theta_min {hi:.7f} h {hb}", flush=True)
