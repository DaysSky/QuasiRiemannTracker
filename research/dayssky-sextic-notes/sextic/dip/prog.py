# Dip realised by a vertical arithmetic progression of n = c0 log U zeros at real part a (constant W~/g).
# Local part f(xi) = exp[rho_1..rho_n](xi) (divided difference) => |f| = |2 sin(Delta xi/2)/Delta|^{n-1}/(n-1)!.
# In r = xi/log U: deficiency eta(r) = -c0*log|sin(phi(r))|, phi linear, maxima (|sin|=1) at r1 <= r*, r2 >= t*.
# Adversary survives iff min_r [cost(r) + 2 eta(r)] >= R on (r1, r2) (cost = A_I for r<=1, long for r>=1).
import numpy as np
dl, q, R, rs, ts = 0.466, 0.233, 0.5935, 0.745, 1.162
def cost(r):
    return np.where(r <= 1, 1 - q - (dl - q)*r, 1 - dl + (5/6 - dl)*(r - 1))
def worst(c0, e0, e1):
    r1, r2 = rs - e0, ts + e1
    r = np.linspace(r1 + 1e-6, r2 - 1e-6, 20001)
    phi = np.pi/2 + np.pi*(r - r1)/(r2 - r1)
    eta = -c0*np.log(np.abs(np.sin(phi)) + 1e-300)
    # outside [r*, t*] the witness costs are >= R anyway (A_I decreasing, long increasing), but include them
    return (cost(r) + 2*eta).min(), eta[np.argmin(np.abs(r - rs))], eta[np.argmin(np.abs(r - ts))]
for c0 in (0.005, 0.01, 0.02, 0.03, 0.05, 0.08):
    best = None
    for e0 in np.linspace(0, 0.3, 61):
        for e1 in np.linspace(0, 0.3, 61):
            m, a0, a1 = worst(c0, e0, e1)
            # the endpoint deficiencies make the r*, t* witnesses cost R + 2*eta: adversary's true count
            cnt = min(m, R + 2*a0, R + 2*a1)
            if best is None or cnt > best[0]: best = (cnt, e0, e1, a0, a1)
    print(f"c0={c0:.3f}: best adversary count {best[0]:.4f} (R={R}), e0={best[1]:.3f} e1={best[2]:.3f}, eta(r*)={best[3]:.4f} eta(t*)={best[4]:.4f}")
