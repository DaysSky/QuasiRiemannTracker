# Cycle61 row count, with the plain-witness slot capacity upgraded to Theorem 2.2 where it applies.
import numpy as np, sys
Lam = 19312/4787; B = 14200/14361 - 1e-5
def zP(m, kap, mode):
    old = np.maximum((1 - 2*m)/(6*kap), 0)
    if mode == 'old': return old
    new = (B - 2*m)/Lam
    ok = (m >= 0.41) & (m <= 0.44) & (new >= 0.034) & (new <= 0.04)
    return np.where(ok, np.maximum(old, new), old)
def Rcount(dl, p, kap, mode, nt=241, nr=401):
    q = dl*p; best = 9; arg=None
    for t in np.linspace(1, 1.5, nt):
        rs = np.linspace(max(t-0.5, 0), min(t, 1), nr); m = t - rs
        zi = np.minimum((1 - rs)/2, 1.0)
        AI = 1 - q - (dl - q)*rs            # inverse: 1 - delta r - 2 q zI with zI=(1-r)/2
        AP = 1 - 2*dl*m - 2*q*zP(m, kap, mode)
        sh = np.minimum(AI, AP); i = sh.argmax(); s = sh[i]
        lng = 1 - dl + (5/6 - dl)*(t - 1)
        v = max(s, lng)
        if v < best: best = v; arg = (round(t,4), round(rs[i],4), round(m[i],4))
    return best, arg
if __name__ == "__main__":
    kap = 0.749722
    for dl, p in [(0.475,0.5),(0.45,0.5),(0.5,0.5),(0.475,0.4),(0.4,0.5),(0.55,0.5)]:
        print(dl, p, 'old', Rcount(dl,p,kap,'old'), 'new', Rcount(dl,p,kap,'new'))
