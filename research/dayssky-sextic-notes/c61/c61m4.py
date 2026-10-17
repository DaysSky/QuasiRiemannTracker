# Cycle61 row count with an optional Moebius fourth moment (r <= r4) or density hypothesis.
import numpy as np, sys
Lam = 19312/4787; B = 14200/14361 - 1e-5
def Rcount(dl, p, kap, r4=None, nt=121, nr=121):
    q = dl*p
    best = 9
    for t in np.linspace(1, 1.5, nt):
        rs = np.linspace(max(t-0.5, 0), min(t, 1), nr)
        m = t - rs
        AI = 1 - q - (dl - q)*rs
        AP = 1 - q/(3*kap) - 2*(dl - q/(3*kap))*m
        short = np.minimum(AI, AP)
        if r4 is not None:
            AI4 = np.where(rs <= r4, 1 - q - (2*dl - q)*rs, 9)   # 4th-moment witness count (slot gain kept)
            short = np.minimum(short, AI4)
        s = short.max()
        lng = 1 - dl + (5/6 - dl)*(t - 1)
        best = min(best, max(s, lng))
    return best
def Rformula(dl, p, kap):
    V = 5/6 - dl
    Ck = 2*(1-p)*(1-p/(3*kap))/(3-p-2*p/(3*kap))
    return 1 - dl + dl*Ck*V/(2*(V + dl*Ck))
def Gmax(th, h, ell, y, r4=None, nd=31, np_=11):
    kap = 2*th - 1
    best = -9
    for dl in np.linspace(0, kap, nd):
        for p in np.linspace(0, 0.5, np_):
            R = Rformula(dl, p, kap) if r4 == 'formula' else Rcount(dl, p, kap, None if r4 == 'none' else r4)
            for d in (h - 1e-3, h):
                Cd = (2*y - (2-B)*d)/Lam
                G = (1+dl)/2 + 13*h/75 - y/2 - ell/2 + ell*dl*p + d*(0.25 + 0.75*R - 17/50) - (dl*p/2)*Cd
                best = max(best, G)
    return best
if __name__ == "__main__":
    kap = 0.749722
    for dl, p in [(0.475, 0.5), (0.3, 0.2), (0.6, 0.4), (0.7, 0.0)]:
        print(dl, p, "formula", round(Rformula(dl,p,kap),5), "count", round(Rcount(dl,p,kap),5), "M4(.75)", round(Rcount(dl,p,kap,0.75),5), "M4(1)", round(Rcount(dl,p,kap,1.0),5))
