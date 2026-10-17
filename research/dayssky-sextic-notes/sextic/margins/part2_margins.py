# Exact re-check of every Part II margin of the OpenAI paper (Sec. 20) at the Conditional-Theorem-D parameters.
# Geometry: B = sigma0 = 11/12 - l/4, h = 1 - x + l, l_y = 1 - l - x, l_x = x. Hypothesis DH_lambda: R <= 1 + lam - delta
# on scales d in [d_lo, h + zeta]; the paper's own no-prime count R = 76/75 - 2 delta/3 (Prop. 19.2, t = 1) below d_lo.
# Requirement (paper (20.12)): E(d) < 0 for delta <= 2 sigma0 - 1, and E(d) < (1+delta)/2 - sigma0 for larger delta (<= 5/6).
from fractions import Fraction as F
import sys
def check(lam, l, x, d_lo, dmin=F(1, 100), zeta=F(1, 1000), a0=None):
    B = F(11, 12) - l/4; h = 1 - x + l; ly = 1 - l - x
    kap_max = F(5, 6)
    out = {}
    def E(d, dl, q, R):
        a = (1 + dl)/2
        return a - B + h*(F(17, 50) - F(1, 6)) - a*ly - l/2 + q*l + d*(R + dl/2 - F(17, 50))
    def need(dl):  # allowed upper bound for E at bin dl (strictly below)
        return max(F(0), (1 + dl)/2 - B)
    worst = {}
    grid = [F(i, 600) for i in range(1, 501)]   # delta in (0, 5/6]
    for dl in grid:
        q = dl/2
        # DH range d in [d_lo, h + zeta]: E linear in d, check both ends
        for d in (d_lo, h + zeta):
            m = E(d, dl, q, 1 + lam - dl) - need(dl)
            worst['DH d in [d_lo,h+zeta]'] = max(worst.get('DH d in [d_lo,h+zeta]', F(-10)), m)
        # paper's count below d_lo
        for d in (dmin, d_lo):
            m = E(d, dl, q, F(76, 75) - 2*dl/3) - need(dl)
            worst['paper count d in [dmin,d_lo]'] = max(worst.get('paper count d in [dmin,d_lo]', F(-10)), m)
    # floor bin: R = 1, delta0 = 2 a0 - 1, slope positive so d = h + zeta bounds it
    d0 = 2*a0 - 1
    worst['floor'] = E(h + zeta, d0, d0/2, F(1)) - 0
    # small rows (20.6)+(10.22): h(17/50 - 1/6) - ly/2 + 2 dmin
    worst['small rows'] = h*(F(17, 50) - F(1, 6)) - ly/2 + 2*dmin
    # principal contour margins (20.3) must be positive
    worst['-m_w = -ly/20'] = -ly/20
    worst['-m_z = -h/600'] = -h/600
    # prime supply at the top scale: l/(h+zeta) > 7/37 (Sec. 20.5)
    worst['7/37 - supply'] = F(7, 37) - l/(h + zeta)
    # low side (REPORT 3): completed-row excess (5l - 1 + d)/8 - d <= 0 on [0, l]  and Gram b = 1 - l - 2x >= 0
    worst['low-side excess (max over d in [0,l])'] = max((5*l - 1 + d)/8 - d for d in (F(0), l))
    worst['-Gram b'] = -(1 - l - 2*x)
    # Prop. 15.3 proof, length inequalities at d = l (worst): X' >= 1, Y'^{11/6}/Pa >> Z^{1/12}, M' > 0
    b = 1 - l - 2*x
    worst["-(l_x - l)  (X' >= 1)"] = -(x - l)
    worst['-(l_y - l - 11b/6 - 1/12)'] = -(ly - l - 11*b/6 - F(1, 12))
    worst["-(M' = 1 - 3l)"] = -(1 - 3*l)
    return B, h, ly, worst
if __name__ == "__main__":
    for lam in (F(1, 1000), F(1, 100), F(1, 50), F(1, 25)):
        lsup = (1 - 2*lam)/(5 + 6*lam)
        l = lsup - lam/100
        lo = (13*l - 1 + 12*(1 + l)*lam)/(4 + 12*lam); hi = (1 - l)/2
        x = (lo + hi)/2                     # interior: Gram b > 0 strictly
        a0 = F(1, 2) + lam/8
        B, h, ly, w = check(lam, l, x, F(2, 5), zeta=lam/10**6, a0=a0)
        print(f"lambda={float(lam)}: l={float(l):.6f} x={float(x):.6f} B=theta={float(B):.6f} h={float(h):.4f} l_y={float(ly):.4f} floor a0={float(a0)}")
        for k, v in w.items():
            print(f"   {k:42s} {float(v):+.7f}  {'OK' if v < 0 else 'FAIL'}")

def validate_against_paper():
    # the same E reproduces the paper's own numbers at its geometry (l = 1/6, x = 17/48, sigma0 = 7/8)
    l, x = F(1, 6), F(17, 48); ly = 1 - l - x; h = 1 - x + l; B = F(7, 8)
    def E(d, dl, q, R):
        a = (1 + dl)/2
        return a - B + h*(F(17, 50) - F(1, 6)) - a*ly - l/2 + q*l + d*(R + dl/2 - F(17, 50))
    assert all(E(F(1, 2), dl, dl/2, F(76, 75) - 2*dl/3) == F(-529, 2400) + 25*dl/96 for dl in (F(0), F(1, 3), F(5, 6)))  # (20.11)
    assert E(h, F(1, 50), F(1, 100), F(1)) == F(-7, 1200)                                                         # (20.5)
    assert h*(F(17, 50) - F(1, 6)) - ly/2 == F(-79, 800)                                                           # (20.6)
    assert ly/20 == F(23, 960) and h/600 == F(13, 9600)                                                            # (20.3)
    print("formulas reproduce paper (20.3), (20.5), (20.6), (20.11) exactly")
validate_against_paper()
