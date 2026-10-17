# Faithful-ish model of the Cycle61 endpoint (paper sections 3-5).
import numpy as np
Lam = 19312/4787; B = 14200/14361 - 1e-5
def Gmax(th, h, ell, y, dgrid=None, nd=61, np_=81):
    kap = 2*th - 1
    dl = np.linspace(0, kap, nd)[:,None]; p = np.linspace(0, 0.5, np_)[None,:]
    V = 5/6 - dl
    Ck = 2*(1-p)*(1-p/(3*kap))/(3-p-2*p/(3*kap))
    R = 1 - dl + dl*Ck*V/(2*(V + dl*Ck))
    best = -1e9; arg=None
    dm = h - 1e-3
    ds = [0.5, dm, h] if dgrid is None else dgrid
    for d in ds:
        Cnew = (2*y - (2-B)*d)/Lam
        Cold = min(d, 2*y-d)/(6*kap)
        Cd = Cnew if d >= dm else Cold
        G = (1+dl)/2 + 13*h/75 - y/2 - ell/2 + ell*dl*p + d*(0.25 + 0.75*R - 17/50) - (dl*p/2)*Cd
        g = G.max()
        if g > best:
            i = np.unravel_index(G.argmax(), G.shape); best = g; arg = (d, dl[i[0],0], p[0,i[1]])
    return best, arg
def low(ell, x, y):
    M = x + y
    return 5/6 - ell/6 + M/12 + max(0, 1 + 3*ell - 2*M)/8
if __name__ == "__main__":
    th = 874861/1e6; h = 3969737987/5e9; ell = 125417/750000
    x = 1 + ell - h; y = h - 2*ell
    print("x,y,M+ell", x, y, x+y+ell)
    print("low", low(ell, x, y), "theta", th)
    print("Gmax", Gmax(th, h, ell, y))
