import numpy as np, sys
from prec import R, Rformula, Lam, B
from scipy.optimize import minimize
def G(d, dl, p, th, h, ell, y, mode, Cmode):
    kap = 2*th - 1
    Rv = Rformula(dl, p, kap) if mode == 'formula' else R(dl, p, kap, mode)
    Cd = (2*y - (2-B)*d)/Lam if Cmode == 'new' else min(d, 2*y - d)/(6*kap)
    return (1+dl)/2 + 13*h/75 - y/2 - ell/2 + ell*dl*p + d*(0.25 + 0.75*Rv - 17/50) - (dl*p/2)*Cd
def Gmax(th, h, ell, y, mode, d, Cmode, coarse=(31,11)):
    kap = 2*th - 1
    best = (-9, None)
    for dl in np.linspace(0, kap, coarse[0]):
        for p in np.linspace(0, 0.5, coarse[1]):
            g = G(d, dl, p, th, h, ell, y, mode, Cmode)
            if g > best[0]: best = (g, (dl, p))
    dl0, p0 = best[1]
    f = lambda v: -G(d, min(max(v[0],0),kap), min(max(v[1],0),0.5), th, h, ell, y, mode, Cmode)
    res = minimize(f, [dl0, p0], method='Nelder-Mead', options={'xatol':1e-9, 'fatol':1e-12, 'maxiter':400})
    return max(best[0], -res.fun), res.x
if __name__ == "__main__":
    th = 874861/1e6; h = 3969737987/5e9; ell = 125417/750000; y = h - 2*ell
    for mode in ['formula', 'old', 'new']:
        g, arg = Gmax(th, h, ell, y, mode, h, 'new')
        print(mode, 'new row Gmax', g, 'slack', th - g, arg)
    for mode in ['old','new']:
        g, arg = Gmax(th, h, ell, y, mode, h - 1e-3, 'old')
        print(mode, 'old row d_- Gmax', g, 'slack', th - g, arg)
