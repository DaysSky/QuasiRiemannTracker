import numpy as np, sys
from fastR import R, Lam, B, short
from scipy.optimize import minimize, minimize_scalar, brentq
def G(d, dl, p, th, h, ell, y, Cmode, mode, kap=None):
    kap = 2*th - 1 if kap is None else kap
    Rv = R(dl, p, kap, mode)
    Cd = (2*y - (2-B)*d)/Lam if Cmode == 'new' else min(d, 2*y - d)/(6*kap)
    return (1+dl)/2 + 13*h/75 - y/2 - ell/2 + ell*dl*p + d*(0.25 + 0.75*Rv - 17/50) - (dl*p/2)*Cd
def dense_max(th, h, ell, y, d, Cmode, mode, kap=None, n1=301, n2=101):
    kap0 = 2*th - 1
    best = (-9, None)
    for dl in np.linspace(0, kap0, n1):
        for p in np.linspace(0, 0.5, n2):
            g = G(d, dl, p, th, h, ell, y, Cmode, mode, kap)
            if g > best[0]: best = (g, (dl, p))
    dl0, p0 = best[1]
    f = lambda v: -G(d, min(max(v[0],0),kap0), min(max(v[1],0),0.5), th, h, ell, y, Cmode, mode, kap)
    res = minimize(f, [dl0, p0], method='Nelder-Mead', options={'xatol':1e-11,'fatol':1e-14,'maxiter':2000})
    v = res.x; v = (min(max(v[0],0),kap0), min(max(v[1],0),0.5))
    return max(best[0], -res.fun), v
th = float(sys.argv[1]); h = float(sys.argv[2])
ell = 11/3 - 4*th; x = 1 + ell - h; y = h - 2*ell; M = x + y
print("theta", th, "ell", ell, "h", h, "x", x, "y", y, "M+ell", M+ell)
for (d, Cm, v) in [(h, 'new', 1e-6), (h - 1e-3, 'old', 1e-5), (0.5, 'old', 1e-3)]:
    g, arg = dense_max(th, h, ell, y, d, Cm, 'new')
    print(f"row d={d:.6f} {Cm}: Gmax={g:.10f}  margin(th - v - G)={th - v - g:.3e} at (delta,p)={arg}")
# binding configuration detail
g, (dl, p) = dense_max(th, h, ell, y, h, 'new', 'new', n1=121, n2=41)
kap = 2*th - 1
lng = lambda t: 1 - dl + (5/6 - dl)*(t - 1)
t = brentq(lambda t: short(t, dl, p, kap, 'new') - lng(t), 1.0, 1.5)
print("binding delta,p", dl, p, "t", t, "R", lng(t))
# geometry side conditions (Cycle61 sec 3-5)
dm, dp, t0 = h - 1/1000, h + 1e-10, 1/1000
n_lo, n_hi = 1 - y/dm, 1 - y/dp + t0
zA = lambda n: (B - 2*n)/Lam
print("n range", n_lo, n_hi, "need (0.41, 264/625=0.4224]", n_lo > 0.41 and n_hi < 264/625)
print("zA range", zA(n_hi), zA(n_lo), "need (0.034, 0.037)", zA(n_hi) > 0.034 and zA(n_lo) < 0.037)
print("supply: ell - d+ zmax", ell - dp*zA(n_lo) > 0, "ell/d+ - 4/21", ell/dp - 4/21 > 0, "2ell-1/3", 2*ell - 1/3 > 0, "4ell-y", 4*ell - y > 0)
print("order: ell<x<y<2x", ell < x < y < 2*x, " y<1/2<h<2y", y < 0.5 < h < 2*y, " 1/6<ell<1/5", 1/6 < ell < 1/5, " h>4ell", h > 4*ell)
print("low side: y-2x+ell<0", y - 2*x + ell < 0, " 1+3ell-2M<0", 1 + 3*ell - 2*M < 0)
print("(4.6) reserve", th - (7/8 - 49*h/150 + 7*ell/8 + 33*(1/200)/50), " (4.7) floor reserve", th - (0.51 + h/3 + 0.51*ell))
print("(5.14)", y/20 - (1+h)/1000 - 1/4000, h/600 - 1/1000 - 1/4000, " (5.15)", 13*h/75 - y/2 + 63*(1/200)/50)
b = 2*y/dm - 1
print("b=2y/d-1", b, "<4/25", b < 4/25)
zold = lambda tt: max(0, (b - 2*tt)/4.5)
Sb = lambda tt: tt/2 - 0.25*((b - (1 - B))/Lam - zold(tt))
print("(3.8) S_b(t0)", Sb(t0))
print("(4.5) junction C_new-C_old at d_-", (2*y - (2-B)*dm)/Lam - min(dm, 2*y - dm)/(6*(2*th-1)))
