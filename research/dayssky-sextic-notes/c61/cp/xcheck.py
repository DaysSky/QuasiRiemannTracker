# Independent float cross-check of a certified point: grid max of G_d(delta,p) - (theta - v) using fastRs (separate code from certG).
import sys, numpy as np
sys.path.insert(0, '.')
th, h, CI, ALPHA, CPW = map(float, sys.argv[1:6])
import fastRs as R
R.CI = CI; R.ALPHA = ALPHA; R.CPW = CPW
if CPW < 3.9586: R.B, R.Lam = 1.0, CPW; R.WIN.update(ZLO=0.001, ZHI=0.30, NLO=0.0, NHI=0.5)
else: R.B, R.Lam = 14200/14361 + 7/10000 - 1e-5, 19312/4787
ell = 11/3 - 4*th; y = h - 2*ell; kap = 2*th - 1; R.LSUP = ell/h
def G(d, dl, p, Cd):
    return (1+dl)/2 + 13*h/75 - y/2 - ell/2 + ell*dl*p + d*(0.25 + 0.75*R.R(dl, p, kap) - 17/50) - (dl*p/2)*Cd
worst = []
for d, v in ((h, 1e-5), (h - 1e-3, 1e-5), (0.5, 1e-3)):
    Cd = min(d, 2*y - d)/CPW
    best = (-9, None)
    for dl in np.linspace(0, 0.75, 301):
        for p in np.linspace(0, 0.5, 51):
            g = G(d, dl, p, Cd) - (th - v)
            if g > best[0]: best = (g, (round(dl,4), round(p,4)))
    worst.append((round(d,4), best))
print(f"theta={th} h={h} CI={CI} ALPHA={ALPHA} CPW={CPW}: max(G - (theta - v)) per row:", worst)
