# Exact piecewise-linear row count. mode 'old' = Cycle61 (3.10); 'new' = plain-witness capacity from Theorem 2.2 where valid.
from scipy.optimize import brentq
Lam = 19312/4787; B = 14200/14361 - 1e-5
def pieces(t, kap, mode):
    """list of (r_lo, r_hi, a, b) with AP(r) = a + b r on [r_lo, r_hi]"""
    lo, hi = max(t - 0.5, 0.0), min(t, 1.0)
    out = []
    def old_piece(rl, rh):
        # AP = 1 - 2 dl m - 2 q (1-2m)/(6k), m = t - r  -> coefficients filled later (need dl, q)
        out.append((rl, rh, 'old'))
    if mode == 'old':
        old_piece(lo, hi); return out
    # new valid where m in [0.41, mmax], mmax from old capacity >= 0.034 and new z in [.034,.04]
    m_hi = min(0.44, (1 - 0.034*6*kap)/2, (B - 0.034*Lam)/2)
    m_lo = max(0.41, (B - 0.04*Lam)/2)
    r_a, r_b = t - m_hi, t - m_lo     # new region in r
    cuts = sorted(set([lo, hi] + [c for c in (r_a, r_b) if lo < c < hi]))
    for i in range(len(cuts)-1):
        rl, rh = cuts[i], cuts[i+1]; mid = (rl+rh)/2; m = t - mid
        out.append((rl, rh, 'new' if (m_lo <= m <= m_hi and m_lo < m_hi) else 'old'))
    return out
def short(t, dl, p, kap, mode):
    q = dl*p
    cI, eI = 1 - q, -(dl - q)                      # AI = cI + eI r
    best = -9
    for rl, rh, kind in pieces(t, kap, mode):
        if kind == 'old':
            # AP = 1 - 2dl(t-r) - 2q(1-2(t-r))/(6k) = [1 - 2dl t - 2q(1-2t)/(6k)] + [2dl - 4q/(6k)] r
            a = 1 - 2*dl*t - 2*q*(1 - 2*t)/(6*kap); b = 2*dl - 4*q/(6*kap)
            # but use max(old,new)? old region: just old
        else:
            # new capacity z = (B - 2m)/Lam (it is >= old there by construction? take max pointwise below)
            a = 1 - 2*dl*t - 2*q*(B - 2*t)/Lam; b = 2*dl - 4*q/Lam
            a_o = 1 - 2*dl*t - 2*q*(1 - 2*t)/(6*kap); b_o = 2*dl - 4*q/(6*kap)
        cands = [rl, rh]
        if abs(eI - b) > 1e-15:
            rc = (a - cI)/(eI - b)
            if rl <= rc <= rh: cands.append(rc)
        if kind == 'new' and abs(b - b_o) > 1e-15:
            rx = (a_o - a)/(b - b_o)
            if rl <= rx <= rh: cands.append(rx)
        for r in cands:
            AP = a + b*r
            if kind == 'new': AP = min(AP, a_o + b_o*r)   # use the better (smaller) count
            best = max(best, min(cI + eI*r, AP))
    return best
def R(dl, p, kap, mode):
    lng = lambda t: 1 - dl + (5/6 - dl)*(t - 1)
    g = lambda t: short(t, dl, p, kap, mode) - lng(t)
    if g(1.0) <= 0: return max(short(1.0, dl, p, kap, mode), lng(1.0)) if False else lng(1.0)
    if g(1.5) > 0: return short(1.5, dl, p, kap, mode)
    t = brentq(g, 1.0, 1.5, xtol=1e-14)
    return lng(t)
if __name__ == "__main__":
    from prec import Rformula
    kap = 2*0.874861 - 1
    for dl, p in [(0.475,0.5),(0.47,0.5),(0.45,0.5),(0.3,0.2),(0.6,0.4),(0.7,0.0)]:
        print(dl, p, round(Rformula(dl,p,kap),10), round(R(dl,p,kap,'old'),10), round(R(dl,p,kap,'new'),10))
