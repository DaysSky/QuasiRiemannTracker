import numpy as np
from dhneed import Rtab
# For target theta: the maximal admissible NEW (q-independent) row count R_new(delta) such that min(R_existing, R_new) satisfies the high side.
def need_curve(th, x):
    ell = 4*(11/12 - th); kappa = 2*th - 1
    h = 1 - x + ell; C0 = -x/3 - 1/12 + 13*ell/12
    d, xs, R = Rtab(ell, x, kappa)
    Qv = xs[None,:]*d[:,None]
    Rallow = 1 - (C0 + (0.5+ell)*d[:,None] + Qv*ell)/h     # need R(delta,q) < Rallow
    fail = R >= Rallow
    need = np.where(fail, Rallow, np.inf).min(axis=1)     # R_new must be below this where existing fails
    return d, need
for th, xs in [(0.8745, np.linspace(0.33,0.40,15)), (0.874, np.linspace(0.33,0.41,17)), (0.872, np.linspace(0.34,0.42,17)), (0.870, np.linspace(0.36,0.43,15))]:
    best=None
    for x in xs:
        d, need = need_curve(th, x)
        if not np.all(np.isfinite(need) <= True): pass
        sig = (1+d)/2
        A = np.where(np.isfinite(need), 2*need/(1-d), np.inf)   # density exponent: U^{A(1-sigma)}, 1-sigma=(1-d)/2
        worstA = A.min()
        if best is None or worstA > best[0]: best=(worstA, x, d, need, A)
    worstA, x, d, need, A = best
    print(f"theta={th}: best x={x:.4f}; most demanding point needs density exponent A <= {worstA:.4f} (DH = 2)")
    for dd in (0.15,0.2,0.3,0.4,0.5,0.6,0.7):
        i = np.argmin(abs(d-dd))
        s = f"A<={A[i]:.3f}  (R<={need[i]:.4f}, DH R={1-d[i]:.4f})" if np.isfinite(A[i]) else "not needed (existing counts suffice)"
        print(f"    sigma={(1+d[i])/2:.3f}: {s}")
