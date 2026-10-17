# Conditional Theorem D with a floor. High side (REPORT §3 / OpenAI): E(h) = C0 + (1/2+l)d + q l - h(1-R), C0 = -x/3 - 1/12 + 13 l/12,
# h = 1 - x + l, q <= d/2. Per-bin count R(d,x): DH (R = 1-d) for d >= dDH; existing closed form (alpha=5/6) for df < d < dDH;
# trivial R = 1 at the floor d = df. theta = 11/12 - l/4 maximised over l, x with E < 0 for all bins and q in [0, d/2].
import numpy as np
def Rex(d, p, kap, al=5/6):
    c = 1/(3*kap); D = 3-(1+2*c)*p; P = (2-2*c*p)*(1-p); J = (al-d)*D + d*P
    return 1 - d + (al-d)*d*P/(2*J)
def ok(l, x, dDH, df, kap):
    h = 1 - x + l; C0 = -x/3 - 1/12 + 13*l/12
    ds = np.concatenate([[df], np.linspace(df, kap, 400)[1:]])
    for i, d in enumerate(ds):
        for p in np.linspace(0, 0.5, 11):
            q = d*p
            if i == 0: R = 1.0
            elif d >= dDH: R = min(1 - d, Rex(d, p, kap))
            else: R = Rex(d, p, kap)
            if C0 + (0.5+l)*d + q*l - h*(1-R) >= 0: return False
    return True
def theta(dDH, df):
    lo, hi = 0.15, 0.2
    for _ in range(30):
        l = (lo+hi)/2; th = 11/12 - l/4; kap = 2*th - 1
        feas = any(ok(l, x, dDH, df, kap) for x in np.linspace(0.25, (1-l)/2, 60))
        if feas: lo = l
        else: hi = l
    return 11/12 - lo/4
print("no DH, floor .02 (paper):", round(theta(9, 0.02), 6))
for df in (0.02, 0.01, 0.005, 0.002):
    for dDH in (0.02, 0.01, 0.0):
        if dDH < df and dDH > 0: continue
        print(f"floor sigma={0.5+df/2:.4f} DH from sigma={0.5+max(dDH,df)/2:.4f}: theta={theta(max(dDH,df), df):.6f}")
