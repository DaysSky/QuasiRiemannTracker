CPW = 6*0.749672
CI = 2.0
ALPHA = 5/6
# fastR with configurable Theorem 2.2 validity window: z in [ZLO, ZHI], n in [NLO, NHI].
from scipy.optimize import brentq
Lam = 19312/4787; B = 14200/14361 - 1e-5
WIN = dict(ZLO=.034, ZHI=.04, NLO=.41, NHI=.44)
def pieces(t, kap):
    lo, hi = max(t - 0.5, 0.0), min(t, 1.0)
    m_hi = min(WIN['NHI'], (1 - WIN['ZLO']*CPW)/2, (B - WIN['ZLO']*Lam)/2)
    m_lo = max(WIN['NLO'], (B - WIN['ZHI']*Lam)/2)
    r_a, r_b = t - m_hi, t - m_lo
    cuts = sorted(set([lo, hi] + [c for c in (r_a, r_b) if lo < c < hi]))
    out = []
    for i in range(len(cuts)-1):
        rl, rh = cuts[i], cuts[i+1]; m = t - (rl+rh)/2
        out.append((rl, rh, 'new' if (m_lo <= m <= m_hi and m_lo < m_hi) else 'old'))
    return out
def short(t, dl, p, kap):
    q = dl*p; cI, eI = 1 - 2*q/CI, -(dl - 2*q/CI); best = -9
    for rl, rh, kind in pieces(t, kap):
        a_o = 1 - 2*dl*t - 2*q*(1 - 2*t)/(CPW); b_o = 2*dl - 4*q/(CPW)
        if kind == 'old': a, b = a_o, b_o
        else: a = 1 - 2*dl*t - 2*q*(B - 2*t)/Lam; b = 2*dl - 4*q/Lam
        cands = [rl, rh]
        if abs(eI - b) > 1e-15:
            rc = (a - cI)/(eI - b)
            if rl <= rc <= rh: cands.append(rc)
        if abs(eI - b_o) > 1e-15:
            rc = (a_o - cI)/(eI - b_o)
            if rl <= rc <= rh: cands.append(rc)
        if kind == 'new' and abs(b - b_o) > 1e-15:
            rx = (a_o - a)/(b - b_o)
            if rl <= rx <= rh: cands.append(rx)
        for r in cands:
            AP = a + b*r
            if kind == 'new': AP = min(AP, a_o + b_o*r)
            best = max(best, min(cI + eI*r, AP))
    return best
def R(dl, p, kap):
    lng = lambda t: 1 - dl + (ALPHA - dl)*(t - 1)
    g = lambda t: short(t, dl, p, kap) - lng(t)
    if g(1.0) <= 0: return lng(1.0)
    if g(1.5) > 0: return short(1.5, dl, p, kap)
    return lng(brentq(g, 1.0, 1.5, xtol=1e-14))
