# Theorem 2.2 high bins with slot-amplitude classes: child with zero parameter delta and slot amplitude p (q = p*delta)
# H(delta,p) = 1 - R(delta,p) - delta*n - 2*p*delta*z - a*(5 - 6R(delta,p)), R = Cycle61 closed-form count.
# Low bins unchanged: a(1-4q0/3) >= (2q0/9) g, g = n + 4.5 z - 1.  Find max capacity z at n = 2m.
import numpy as np
kap = 0.7497
def R(d, p):
    V = 5/6 - d; C = 2*(1-p)*(1-p/(3*kap))/(3-p-2*p/(3*kap))
    return 1 - d + d*C*V/(2*(V + d*C))
def feasible(n, z, q0, amp_classes=True):
    g = n + 4.5*z - 1
    if g <= 0: return True
    a = (2*q0/9)*g/(1 - 4*q0/3)            # minimal amplifier paying the low bins
    ds = np.linspace(q0, 0.75, 40); ps = np.linspace(0, 0.5, 51) if amp_classes else np.array([0.0])
    for d in ds:
        for p in ps:
            Rv = R(d, p)
            slot = 2*p*d*z if amp_classes else d*z
            if 1 - Rv - d*n - slot - a*(5 - 6*Rv) < 0: return False
    return True
def cap(m, amp):
    best = (0, None)
    for q0 in np.linspace(0.55, 0.74, 39):
        lo, hi = 0.0, 0.2
        for _ in range(40):
            mid = (lo+hi)/2
            if feasible(2*m, mid, q0, amp): lo = mid
            else: hi = mid
        if lo > best[0]: best = (lo, q0)
    return best
for m in (0.40, 0.417, 0.43):
    print(f"m={m}: no classes cap={cap(m, False)}, with amplitude classes cap={cap(m, True)}  old6k={(1-2*m)/4.5:.4f}", flush=True)
