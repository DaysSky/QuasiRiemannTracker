# Conditional Theorem D with the intermediate row scales d in [d_min, 1/2] (OpenAI/Liu E_B(d) formula):
# E(d) = a - B + h(17/50 - 1/6) - a*l_y - l/2 + q*l + d*(R_d + dl/2 - 17/50), a=(1+dl)/2, B = 11/12 - l/4,
# l_x = x, l_y = 1 - l - x, h = 1 - x + l. Check d = h (count R_h) and d = 1/2 (count R_half); floor bin with R = 1.
import numpy as np
def E(d, dl, q, l, x, R):
    a=(1+dl)/2; B=11/12-l/4; ly=1-l-x; h=1-x+l
    return a-B+h*(17/50-1/6)-a*ly-l/2+q*l+d*(R+dl/2-17/50)
def ok(l, x, Rh, Rhalf, df):
    th=11/12-l/4; kap=2*th-1; h=1-x+l
    for dl in np.concatenate([[df], np.linspace(df, kap, 300)[1:]]):
        for q in np.linspace(0, dl/2, 6):
            floor = (dl == df)
            if E(h, dl, q, l, x, 1.0 if floor else Rh(dl)) >= 0: return False
            if E(0.5, dl, q, l, x, 1.0 if floor else Rhalf(dl)) >= 0: return False
    return True
def theta(Rh, Rhalf, df):
    lo, hi = 0.15, 0.2
    for _ in range(28):
        l=(lo+hi)/2
        if any(ok(l, x, Rh, Rhalf, df) for x in np.linspace(0.2, (1-l)/2, 80)): lo=l
        else: hi=l
    return 11/12-lo/4
DH = lambda dl: 1-dl
Liu = lambda dl: 76/75-2*dl/3
for df in (0.02, 0.002):
    print(f"floor sigma={0.5+df/2}: DH at d=h, Liu count at d<=1/2: theta={theta(DH, Liu, df):.6f};  DH at all d: theta={theta(DH, DH, df):.6f}")
