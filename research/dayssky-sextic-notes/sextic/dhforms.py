# Which single inputs give the DH count R(delta) <= 1 - delta for the block calculus at t = 1?
# Block at length s in [1, 3/2]: rep (i) pieces (r <= 1, m = s - r <= 1/2); rep (ii) pieces (r' in [1, s], m' = s - r').
# Counts available (exponent of U): Mobius LS 1 - delta*r (r<=1); amplified Mobius e(r) - delta*r (r>1);
# plain pair 1 - 2*delta*m; optional hypotheses below.
import numpy as np
def count(delta, s, hyp):
    best_i = -9
    for r in np.linspace(max(s-0.5, 0), 1.0, 401):
        m = s - r
        c = [1 - delta*r, 1 - 2*delta*m]
        if 'MLV' in hyp:  c.append(r*(1-delta))                     # Mobius large values at length U^r, r<=1
        if 'PLV' in hyp:  c.append(1 - delta if m <= 0.5 else 9)     # plain large values (count U^{1-delta} at any m<=1/2)
        if 'MMp' in hyp:  c.append(1 - delta*(r+m))                  # supercritical mixed moment
        best_i = max(best_i, min(c))
    best_ii = -9
    for rp in np.linspace(1.0, s, 401):
        m = s - rp
        e = hyp.get('e', lambda x: (1+5*x)/6)(rp)
        c = [e - delta*rp, 1 - 2*delta*m]
        best_ii = max(best_ii, min(c))
    return min(best_i, best_ii)
def R(delta, hyp):
    return max(count(delta, s, hyp) for s in np.linspace(1.0, 1.5, 101))
hyps = {'none': {}, 'MLV': {'MLV':1}, 'PLV': {'PLV':1}, "MM'": {'MMp':1},
        'long Mobius e(s)=1+0.45(s-1)': {'e': lambda x: 1 + 0.45*(x-1)},
        'long Mobius e(s)=1': {'e': lambda x: 1.0}}
for name, h in hyps.items():
    vals = [(d, R(d, h)) for d in (0.1, 0.3, 0.467, 0.6, 0.75)]
    print(f"{name:32s}", '  '.join(f"d={d}:R={v:.3f}(DH {1-d:.3f})" for d, v in vals))
