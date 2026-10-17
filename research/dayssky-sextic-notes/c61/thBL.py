import sys, fastRw
from fractions import Fraction as F
sys.path.insert(0,'../q0scan')
from fam import BL, rho
def Bfull(q0):
    B0q, Lq, Ka, r0 = BL(q0); c = Ka*(5-6*r0)
    dB = F(7,10000)*F(81379,15925)/(q0+c)   # timbalud saving rescaled to this q0
    return float(B0q + dB) - 1e-5, float(Lq)
Bd, Ld, Bn, Ln = map(float, sys.argv[1:5])

fastRw.B = Bd; fastRw.Lam = Ld
import thoptw_core as T
T.B = Bn; T.Lam = Ln
def parts(th, h):
    ell = 11/3 - 4*th; y = h - 2*ell
    return (th - 1e-6 - T.Gmax(th, h, ell, y, h, 'new'), th - 1e-5 - T.Gmax(th, h, ell, y, h - 1e-3, 'old'),
            th - 1e-3 - T.Gmax(th, h, ell, y, 0.5, 'old'))
from scipy.optimize import minimize_scalar
lo, hi = 0.87476, 0.87490; hb = None
for _ in range(15):
    mid = (lo + hi)/2
    res = minimize_scalar(lambda h: -min(parts(mid, h)), bounds=(0.785, 0.80), method='bounded', options={'xatol':1e-8})
    if -res.fun >= 0: hi = mid; hb = res.x
    else: lo = mid
print(f"Bd={Bd:.6f} Ld={Ld:.4f} Bn={Bn:.6f} Ln={Ln:.4f} theta_min {hi:.8f} h {hb:.6f}", flush=True)
