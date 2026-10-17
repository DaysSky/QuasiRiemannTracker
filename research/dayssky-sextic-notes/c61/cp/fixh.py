# Bisect theta at a fixed h with combo_cert3's window-aware parts (radical constants).
import sys
sys.argv = ['combo_cert3.py', 'rad'] + sys.argv[1:4] + [sys.argv[4]]
CI, ALPHA, CPW, H = map(float, sys.argv[2:6])
sys.path.insert(0, '.')
import fastRs, thcore_s as T
fastRs.CI = CI; fastRs.ALPHA = ALPHA; fastRs.CPW = T.CPW = CPW
B, LAM = 14200/14361 + 7/10000 - 1e-5, 19312/4787
fastRs.B = T.B = B; fastRs.Lam = T.Lam = LAM
def window_ok(th, h, ell, y):
    dm, dp, t0 = h - 1e-3, h + 1e-10, 1e-3
    n_lo = 1 - y/dm; n_hi = 1 - y/dp + t0
    zA = lambda n: (B - 2*n)/LAM
    return 0.41 < n_lo and n_hi < 264/625 and zA(n_hi) > 17/500 and zA(n_lo) < 1/25
def parts(th, h):
    ell = 11/3 - 4*th; y = h - 2*ell
    fastRs.LSUP = ell/h
    mode = 'new' if window_ok(th, h, ell, y) else 'old'
    return mode, (th - 1e-5 - T.Gmax(th, h, ell, y, h, mode), th - 1e-5 - T.Gmax(th, h, ell, y, h - 1e-3, 'old'), th - 1e-3 - T.Gmax(th, h, ell, y, 0.5, 'old'))
lo, hi = 0.870, 0.8752
for _ in range(18):
    mid = (lo + hi)/2
    mode, p = parts(mid, H)
    if min(p) >= 0: hi = mid
    else: lo = mid
print(f"CI={CI} ALPHA={ALPHA} CPW={CPW} h={H}: theta {hi:.6f} mode {parts(hi, H)[0]}", flush=True)
