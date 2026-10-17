# theta when the DH count R <= 1 - delta is granted only on a sigma-window (capc4 grid model, current slot costs, slope 5/6)
import numpy as np, capc4, sys
orig = capc4.Rgrid
def make(s_lo, s_hi):
    d_lo, d_hi = 2*s_lo - 1, 2*s_hi - 1
    def Rgrid(*args):
        dl, p, R = orig(*args)
        mask = (dl >= d_lo - 1e-12) & (dl <= d_hi + 1e-12)
        return dl, p, np.where(mask, np.minimum(R, 1 - dl), R)
    return Rgrid
capc4.Rgrid = orig
base = capc4.theta(4.498/0.7497, 2, 5/6)
print(f"baseline {base:.5f}", flush=True)
for lo, hi in [(0.5, 0.875), (0.6, 0.875), (0.65, 0.875), (0.7, 0.875), (0.72, 0.875), (0.75, 0.875), (0.8, 0.875),
               (0.5, 0.7), (0.5, 0.75), (0.5, 0.8), (0.7, 0.76), (0.68, 0.78), (0.72, 0.75)]:
    capc4.Rgrid = make(lo, hi)
    print(f"DH on sigma in [{lo}, {hi}]: theta = {capc4.theta(4.498/0.7497, 2, 5/6):.5f}", flush=True)
