# Toolkit closure check: for the adversary configuration (D = {r*, t*}, m*) at each delta (p = 1/2),
# evaluate every available input; the minimum must equal R(delta, 1/2).
import sys; sys.path.insert(0, '../c61')
import numpy as np
from fastR import R, short
th = 0.8748389; kap = 2*th - 1; p = 0.5
B_ = 14200/14361 - 1e-5; Lam = 19312/4787
def zP(m):
    z = max(1 - 2*m, 0)/(6*kap)
    if (B_ - .04*Lam)/2 <= m <= min(.44, (B_ - .034*Lam)/2, (1 - .034*6*kap)/2):
        z = max(z, (B_ - 2*m)/Lam)
    return z
for dl in (0.3, 0.4, 0.466, 0.55, 0.65):
    q = dl*p; Rv = R(dl, p, kap, 'new')
    lng = lambda s: (1 + 5*s)/6 - dl*s
    # adversary configuration: t* where short = long, r* the split
    ts = np.linspace(1, 1.5, 20001)
    tstar = ts[np.argmin([abs(short(t, dl, p, kap, 'new') - lng(t)) for t in ts])]
    rs = np.linspace(max(tstar - .5, 0), min(tstar, 1), 20001)
    AI = 1 - q - (dl - q)*rs; AP = np.array([1 - 2*dl*(tstar - r) - 2*q*zP(tstar - r) for r in rs])
    j = int(np.argmax(np.minimum(AI, AP))); rstar = rs[j]; mstar = tstar - rstar
    inputs = {
      'A_I(r*)': AI[j], 'A_P(m*)': AP[j], 'long(t*) [M_{t*} or block M_r*S_m* via copies]': lng(tstar),
      'M4 copies at r*': (1 + 10*rstar)/6 - 2*dl*rstar if 2*rstar > 1 else 1 - 2*dl*rstar,
      'Mobius x Mobius (r*, t*) via copies': lng(rstar + tstar),
      'MM\' theta (needs r+2m<1)': (1 - dl*(rstar + mstar)) if rstar + 2*mstar < 1 else np.inf,
      'M x M marked (needs 3r<=1)': (1 - 2*dl*rstar) if 3*rstar <= 1 else np.inf,
      'pure slots (BGL a0=1/2)': 1 - q,
      'plain 6th via pointwise (loss .375)': 1 + 0.375 - 3*dl*mstar,
    }
    mn = min(inputs.values())
    print(f"delta={dl}: R={Rv:.4f}  t*={tstar:.3f} r*={rstar:.3f} m*={mstar:.3f}  min over inputs = {mn:.4f}")
    for k, v in inputs.items():
        print(f"     {k:55s} {v:.4f}" if np.isfinite(v) else f"     {k:55s} n/a")
