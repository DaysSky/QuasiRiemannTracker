# Exponent model of OpenAI's Part II high side (sextic), generalised.
# theta = 11/12 - ell/4 ; h = 1 - x + ell ; C0 = -x/3 - 1/12 + 13 ell/12
# E = C0 + (1/2+ell) delta + q ell - h (1 - R(delta,q)) must be < 0 for delta in (0,kappa], q in [0,delta/2]
import numpy as np
from functools import lru_cache

def count_exponent(delta, q, t, supply, kappa, alpha=5/6, plainc=6.0, invamp=None, grid=41):
    """Adversarial witness (r,m) with r<=t, r>=t-1/2, 0<=m<=1/2, r+m>=t; we take min over available counts.
    plainc: coefficient c in plain capacity z_P(m) = (1-2m)/(c*kappa)  (OpenAI: 6 -> 6 kappa = 9/2 at kappa=3/4).
    invamp: function e(r) for no-slot inverse; default max(1,(1+5r)/6) i.e. alpha=5/6."""
    worst = -1e9
    for r in np.linspace(max(t-0.5,0), t, grid):
        m = max(t - r, 0.0)
        if m > 0.5 + 1e-12: continue
        e = max(1.0, 1 - alpha + alpha*r) if invamp is None else invamp(r)
        cands = [e - delta*r]                                  # inverse, no slots (amplified)
        if r < 1:
            zM = min((1-r)/2, supply)
            cands.append(1 - delta*r - 2*q*zM)                 # inverse with slots
        zP = min((1-2*m)/(plainc*kappa), supply)
        cands.append(1 - 2*delta*m - 2*q*zP)                   # plain (two copies) with slots
        val = min(cands)
        worst = max(worst, val)
    return worst

def R_best(delta, q, supply, kappa, **kw):
    ts = np.linspace(1.0, 1.5, 51)
    return min(count_exponent(delta, q, t, supply, kappa, **kw) for t in ts)

def high_margin(ell, x, kappa=None, nd=60, nq=6, **kw):
    h = 1 - x + ell
    C0 = -x/3 - 1/12 + 13*ell/12
    theta = 11/12 - ell/4
    if kappa is None: kappa = 2*theta - 1
    supply = ell/h
    worst = -1e9; arg=None
    for delta in np.linspace(1/50, kappa, nd):
        for q in np.linspace(0, delta/2, nq):
            R = R_best(delta, q, supply, kappa, **kw)
            E = C0 + (0.5+ell)*delta + q*ell - h*(1-R)
            if E > worst: worst, arg = E, (round(delta,4), round(q,4), round(R,4))
    return worst, arg

if __name__ == "__main__":
    # calibration: OpenAI geometry ell=1/6, x=17/48, kappa=3/4 (6 kappa = 9/2)
    w,a = high_margin(1/6, 17/48, kappa=3/4)
    print("OpenAI geometry: max E =", w, "at", a, " (paper: -49/440640 = %.6g)" % (-49/440640))
