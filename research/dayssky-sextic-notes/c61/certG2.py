# Exact-rational certificate for the row inequalities with the upgraded plain capacity.
# G_d(delta,p) = E_d(delta,p) + (3/4) d R(delta,p) <= theta - v_d on [0,3/4]x[0,1/2].
# R(delta,p) <= F(t; delta1,p1) for (delta,p) >= (delta1,p1) componentwise (monotonicity), any t in [1,3/2].
from fractions import Fraction as Fr
import sys, math
from scipy.optimize import brentq
F = Fr
LAM = F(19312, 4787); B = F(14200, 14361) + F(7, 10000) - F(1, 100000)
CI = F(sys.argv[3]); ALPHA = F(sys.argv[4]) if len(sys.argv) > 4 else F(5, 6)
LSUP = None   # inverse slot supply ell/d_plus, set in main
def setup(theta, h):
    ell = F(11, 3) - 4*theta; x = 1 + ell - h; y = h - 2*ell; kap = 2*theta - 1
    return dict(theta=theta, h=h, ell=ell, x=x, y=y, kap=kap)
import os as _os
CPW = F(_os.environ['CPW']) if _os.environ.get('CPW') else None
def cpw(kap):
    return CPW if CPW is not None else 6*kap
def zlines(kap):
    # old(m) = (1-2m)/(6 kap); new(m) = (B-2m)/LAM   as (const, slope) in m
    return (F(1)/cpw(kap), F(-2)/cpw(kap)), (B/LAM, F(-2)/LAM)
def Vrange():
    m_lo = max(F(41,100), (B - F(4,100)*LAM)/2)
    m_hi = min(F(44,100), (B - F(34,1000)*LAM)/2, (1 - F(34,1000)*F(9,2))/2)   # kappa=3/4 coverage of short I
    return m_lo, m_hi
def short_exact(t, dl, p, kap):
    q = dl*p
    lo, hi = max(t - F(1,2), F(0)), min(t, F(1))
    (oc, os), (nc, ns) = zlines(kap)
    m_lo, m_hi = Vrange()
    # breakpoints in r
    m_kink = (oc - nc)/(ns - os) if ns != os else None   # old(m)=new(m)
    pts = {lo, hi}
    for mm in (m_lo, m_hi, m_kink):
        if mm is not None:
            r = t - mm
            if lo < r < hi: pts.add(r)
    rk = 1 - CI*LSUP            # inverse capacity (1-r)/CI hits the supply below rk
    if lo < rk < hi: pts.add(rk)
    m_capo = (1 - cpw(kap)*LSUP)/2   # old plain capacity (1-2m)/cpw hits the supply for m < m_capo
    r_capo = t - m_capo
    if lo < r_capo < hi: pts.add(r_capo)
    pts = sorted(pts)
    if len(pts) == 1:
        pts = [pts[0], pts[0]]
    best = None
    AI = lambda r: 1 - 2*q/CI - (dl - 2*q/CI)*r
    for a, b in zip(pts[:-1], pts[1:]):
        mid = (a + b)/2; mm = t - mid
        inV = (m_lo <= mm <= m_hi)
        if inV:
            zc_o = oc + os*mm; zc_n = nc + ns*mm
            (zc, zs) = (nc, ns) if zc_n >= zc_o else (oc, os)
        else:
            (zc, zs) = (oc, os)
        # AP(r) = 1 - 2 dl (t-r) - 2 q (zc + zs (t - r))
        Ac = 1 - 2*dl*t - 2*q*(zc + zs*t); As = 2*dl + 2*q*zs
        if mm < m_capo and (zc, zs) == (oc, os):   # plain slot length capped by the supply
            Ac, As = 1 - 2*dl*t - 2*q*LSUP, 2*dl
        cands = [a, b]
        if mid < rk: Ic, Is = 1 - 2*q*LSUP, -dl          # capped: A_I = 1 - dl r - 2 q LSUP
        else: Ic, Is = 1 - 2*q/CI, -(dl - 2*q/CI)
        if a == b:
            # degenerate: evaluate with both possible z formulas (upper bound)
            for (zc2, zs2) in ((oc, os), (nc, ns)) if inV else ((oc, os),):
                Ac2 = 1 - 2*dl*t - 2*q*(zc2 + zs2*t); As2 = 2*dl + 2*q*zs2
                if (zc2, zs2) == (oc, os) and (t - a) < m_capo:
                    Ac2, As2 = 1 - 2*dl*t - 2*q*LSUP, 2*dl
                v = min(Ic + Is*a, Ac2 + As2*a)
                if best is None or v > best: best = v
            continue
        if As != Is:
            rc = (Ic - Ac)/(As - Is)
            if a < rc < b: cands.append(rc)
        for r in cands:
            v = min(Ic + Is*r, Ac + As*r)
            if best is None or v > best: best = v
    return best
def F_upper(t, dl, p, kap):
    lng = 1 - dl + (ALPHA - dl)*(t - 1)
    return max(short_exact(t, dl, p, kap), lng)
def t_float(dl, p, kap):
    fl = lambda t: float(short_exact(F(t).limit_denominator(10**9), F(dl), F(p), kap)) - (1 - dl + (float(ALPHA) - dl)*(t - 1))
    try:
        if fl(1.0) <= 0: return F(1)
        if fl(1.5) > 0: return F(3,2)
        return F(brentq(fl, 1.0, 1.5, xtol=1e-12)).limit_denominator(10**9)
    except Exception:
        return F(1)
def E(par, d, Cd, dl, p):
    h, y, ell = par['h'], par['y'], par['ell']
    return (1 + dl)/2 + F(13,75)*h - y/2 - ell/2 + ell*dl*p + d*(F(1,4) - F(17,50)) - dl*p*Cd/2
def certify(par, d, Cd, v, maxdepth=40):
    target = par['theta'] - v; kap = par['kap']
    stack = [(F(0), F(3,4), F(0), F(1,2), 0)]
    nbox = 0; worst = None
    while stack:
        d1, d2, p1, p2, depth = stack.pop()
        nbox += 1
        t = t_float(float(d1), float(p1), kap)
        Rup = F_upper(t, d1, p1, kap)
        Eup = max(E(par, d, Cd, a, b) for a in (d1, d2) for b in (p1, p2))
        U = Eup + F(3,4)*d*Rup
        if U <= target:
            gap = target - U
            if worst is None or gap < worst: worst = gap
            continue
        if depth >= maxdepth:
            return False, nbox, (float(d1), float(d2), float(p1), float(p2), float(U - target))
        dm, pm = (d1 + d2)/2, (p1 + p2)/2
        if (d2 - d1) >= 2*(p2 - p1):
            stack += [(d1, dm, p1, p2, depth+1), (dm, d2, p1, p2, depth+1)]
        else:
            stack += [(d1, d2, p1, pm, depth+1), (d1, d2, pm, p2, depth+1)]
    return True, nbox, float(worst)
if __name__ == "__main__":
    theta = F(sys.argv[1]); h = F(sys.argv[2])
    print('CI', CI, 'ALPHA', ALPHA, 'CPW', CPW, 'B', float(B))
    par = setup(theta, h)
    LSUP = par['ell']/(h + F(1, 10**10))
    print('LSUP', float(LSUP), 'rk', float(1 - CI*LSUP))
    print({k: float(v) for k, v in par.items()})
    y, kap = par['y'], par['kap']
    d_minus = h - F(1, 1000)
    import os
    top = ("old d=h", h, min(h, 2*y - h)/cpw(kap), F(1, 10**5)) if os.environ.get('OLDTOP') == '1' else ("new d=h", h, (2*y - (2 - B)*h)/LAM, F(1, 10**6))
    rows = [top,
            ("old d=d-", d_minus, min(d_minus, 2*y - d_minus)/cpw(kap), F(1, 10**5)),
            ("old d=1/2", F(1,2), min(F(1,2), 2*y - F(1,2))/cpw(kap), F(1, 1000))]
    for name, d, Cd, v in rows:
        ok, nb, info = certify(par, d, Cd, v)
        print(name, "certified" if ok else "FAILED", "boxes", nb, "min slack" if ok else "fail box", info, flush=True)
