# Jensen/Cauchy dichotomy at the zero (MODEL, exponent algebra).
# Case A: |L| <= U^{2 eta} on the box Re in [sigma-R, sigma], height window  ->  plain witnesses have m <= 2 eta / R,
#         detector at t = 1 gives count_A = 1 - delta + (delta - q) * m_max  (A_I at r = 1 - m_max).
# Case B: |L(s1)| >= U^{2 eta} at Re s1 = sigma1 in [sigma-R, sigma] -> AFE piece of length U^m, m >= 2 eta/(1-sigma1),
#         plain 4th moment: count_B = 1 - 2(4 eta + m (2 sigma1 - 1)); worst m = 1/2 or minimal m, worst sigma1 = sigma - R.
import numpy as np
def counts(delta, p, R, eta):
    q = delta*p; sigma = (1+delta)/2; s1 = sigma - R
    mmax = min(0.5, 2*eta/R)
    cA = 1 - delta + (delta - q)*mmax
    mmin = 2*eta/(1 - s1)
    cB = max(1 - 2*(4*eta + m*(2*s1 - 1)) for m in np.linspace(mmin, 0.5, 51)) if mmin <= 0.5 else 0
    return cA, cB
for delta in (0.3, 0.466, 0.6):
    best = (9, None)
    for R in np.linspace(0.01, delta/2, 40):
        for eta in np.linspace(0.001, 0.2, 400):
            cA, cB = counts(delta, 0.5, R, eta)
            v = max(cA, cB)
            if v < best[0]: best = (v, (R, eta, cA, cB))
    print(f"delta={delta}: best max(count_A, count_B) = {best[0]:.4f} at R,eta,cA,cB = {tuple(round(x,4) for x in best[1])}")
