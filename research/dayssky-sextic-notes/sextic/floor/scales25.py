# As scales.py, but with the final D' split: DH on d in [2/5, h], the paper's (Liu's) count 76/75 - 2 delta/3 on d <= 2/5.
import numpy as np, importlib.util
spec = importlib.util.spec_from_file_location("sc", "scales.py")
src = open("scales.py").read().split("DH = lambda")[0]
ns = {}; exec(src, ns)
E = ns['E']
def ok(l, x, df):
    th = 11/12 - l/4; kap = 2*th - 1; h = 1 - x + l
    for dl in np.concatenate([[df], np.linspace(df, kap, 300)[1:]]):
        for q in np.linspace(0, dl/2, 6):
            floor = (dl == df)
            R_DH = 1.0 if floor else 1 - dl
            R_paper = max(1.0, 76/75 - 2*dl/3) if floor else 76/75 - 2*dl/3
            if E(h, dl, q, l, x, R_DH) >= 0: return False
            if E(0.4, dl, q, l, x, R_DH) >= 0: return False
            if E(0.4, dl, q, l, x, R_paper) >= 0: return False
    return True
def theta(df):
    lo, hi = 0.15, 0.2
    for _ in range(28):
        l = (lo + hi)/2
        if any(ok(l, x, df) for x in np.linspace(0.2, (1 - l)/2, 80)): lo = l
        else: hi = l
    return 11/12 - lo/4
for df in (0.02, 0.002, 0.0002):
    print(f"floor sigma={0.5+df/2}: DH on d in [2/5,h], paper count on d<=2/5: theta={theta(df):.6f}")
