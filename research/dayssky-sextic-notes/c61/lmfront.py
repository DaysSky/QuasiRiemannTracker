# For each s0, the largest LM slope a' on [1, s0] (5/6 beyond) whose model theta gains >= 0.00008 over the baseline
# (the gap between the model baseline 0.87492 and the record 0.8748389).  Bisection on a'.
import capc4, lmrange, sys
base = 0.87492
for s0 in [float(a) for a in sys.argv[1:]]:
    lo, hi = 0.0, 5/6      # gain(lo) large, gain(hi) = 0
    for _ in range(7):
        mid = (lo + hi)/2
        capc4.Rgrid = lmrange.make_Rgrid(s0, mid)
        th = capc4.theta(4.498/0.7497, 2, 5/6)
        if base - th >= 0.00008: lo = mid
        else: hi = mid
    print(f"s0={s0}: LM slope a' <= {lo:.3f} on [1,{s0}] gains >= 0.00008 (a' = {hi:.3f} does not)", flush=True)
