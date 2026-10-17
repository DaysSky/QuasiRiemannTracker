# Theta residues at prime indices: rho(m) = Res_{s=3/4} sum_c gtilde(conj chi_c, m) N(c)^{-s}  (sqfree c, data file arg1)
# fitted from smoothed sums on independent dyadic blocks (see notheta.py).  Also the alpha-twisted series at m^3 (null test).
import sys, math, numpy as np
fn = sys.argv[1]
sys.argv = [sys.argv[0], fn]
exec(open('notheta.py').read().split('print("data')[0])
G = s.Gt
out = []
for k in range(s.K):
    (ma, mb), Nm, _ = s.tp[k]
    c0, se0, z0, _, _ = fit(G * s.chipow(k, 1))            # untwisted, index m
    c1, se1, z1, _, _ = fit(s.alpha * G * s.chipow(k, 3))  # alpha-twisted, index m^3  (null)
    c3, se3, z3, _, _ = fit(G * s.chipow(k, 3))            # untwisted, index m^3
    out.append((ma, mb, Nm, abs(c0), se0, z0, z1, abs(c3), z3))
A = np.array(out, float)
print("alpha-twisted index m^3 (null, expect E z^2 ~ 1): mean z^2 = %.2f, max z = %.2f, #z>3 = %d of %d"
      % (np.mean(A[:, 6] ** 2), A[:, 6].max(), int(np.sum(A[:, 6] > 3)), len(A)))
print("untwisted index m^3: mean z^2 = %.2f, #z>3 = %d" % (np.mean(A[:, 8] ** 2), int(np.sum(A[:, 8] > 3))))
cls = {}
for ma, mb, Nm, c, se, z, _, _, _ in out:
    key = 'N=5 mod 8' if Nm % 8 == 5 else ('inert' if ma != 0 and mb == 0 else ('b=4 mod 8' if mb % 8 == 4 else 'b=0 mod 8'))
    cls.setdefault(key, []).append((Nm, c, se, z))
for key, v in cls.items():
    v = np.array(v)
    print(f"class {key:10s}: n={len(v):3d}  mean|rho|={v[:,1].mean():.4f}  mean SE={v[:,2].mean():.4f}  #z>3={int(np.sum(v[:,3]>3))}")
    if key != 'N=5 mod 8' and len(v) > 5:
        # bias-corrected |rho|^2 = |c|^2 - SE^2 ; regress log(mean in norm bins) on log N
        bins = np.unique(np.round(np.geomspace(v[:,0].min(), v[:,0].max() + 1, 6)))
        xs, ys = [], []
        for lo, hi in zip(bins[:-1], bins[1:]):
            sel = (v[:,0] >= lo) & (v[:,0] < hi)
            if sel.sum() >= 2:
                xs.append(np.log(np.mean(v[sel,0]))); ys.append(np.log(max(np.mean(v[sel,1]**2 - v[sel,2]**2), 1e-8)))
                print(f"     N in [{lo:6.0f},{hi:6.0f}): n={sel.sum():3d}  mean |rho|^2 (debiased) = {math.exp(ys[-1]):.5f}")
        if len(xs) >= 2:
            sl = np.polyfit(xs, ys, 1)[0]
            print(f"     slope d log|rho|^2 / d log N = {sl:.3f}   => theta_prime ~ {-sl/2:.3f}")
