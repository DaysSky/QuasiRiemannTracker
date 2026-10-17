from fractions import Fraction as F
import sys
from cert import LAM, B, Vrange
theta = F(sys.argv[1]); h = F(sys.argv[2])
ell = F(11,3) - 4*theta; x = 1 + ell - h; y = h - 2*ell; M = x + y; kap = 2*theta - 1
dm, dp, t0, d0 = h - F(1,1000), h + F(1,10**10), F(1,1000), F(1,200)
N = F(264,625)
res = {}
res['geometry ell<x<y<2x'] = ell < x < y < 2*x
res['y<1/2<h<2y'] = y < F(1,2) < h < 2*y
res['1/6<ell<1/5'] = F(1,6) < ell < F(1,5)
res['h>4ell'] = h > 4*ell
res['M+ell<=1'] = M + ell <= 1
res['low: y-2x+ell<0'] = y - 2*x + ell < 0
res['low: 1+3ell-2M<=0'] = 1 + 3*ell - 2*M <= 0
# (3.3): n = 1 - y/d + t on d in [d-, d+], t in [0,t0]
n_lo = 1 - y/dm; n_hi = 1 - y/dp + t0
zA = lambda n: (B - 2*n)/LAM
res['(3.3) 41/100 < n_lo'] = F(41,100) < n_lo
res['(3.3) n_hi < N=264/625'] = n_hi < N
res['(3.3) zA_min > 17/500'] = zA(n_hi) > F(17,500)
res['(3.3) zA_max < 1/25'] = zA(n_lo) < F(1,25)
res['case n<=41/100: 41/100+N+(9/2)zA_max < 1'] = F(41,100) + N + F(9,2)*zA(n_lo) < 1
res['case zI<=17/500: 2N+(9/2)(17/500) < 1'] = 2*N + F(9,2)*F(17,500) < 1
# (3.5) supply
res['(3.5) ell - d+ zmax > 0'] = ell - dp*zA(n_lo) > 0
res['(3.5) ell/d+ - 4/21 > 0'] = ell/dp - F(4,21) > 0
res['(3.5) 2ell - 1/3 > 0'] = 2*ell - F(1,3) > 0
res['(3.5) 4ell - y > 0'] = 4*ell - y > 0
# detector plain capacity <= supply on all d >= 1/2: max plain capacity 0.04 vs L = ell/d >= ell/d+
res['detector capacities <= supply'] = max(F(4,21), F(1,25)) <= ell/dp
# (3.7)-(3.8): b = 2y/d - 1 at d = d- (largest b); S_b(t) = t/2 - 1/4((b-(1-B))/LAM - zold(t))
b = 2*y/dm - 1
zold = lambda t: max(F(0), (b - 2*t)/F(9,2))
Sb = lambda t: t/2 - F(1,4)*((b - (1 - B))/LAM - zold(t))
res['(3.8) S_b(t0) > 0'] = Sb(t0) > 0
# S_b is increasing for t >= t0 (slopes 1/2 - 1/(4*4.5)*... ): check slope pieces positive
res['(3.7) S_b slopes positive'] = (F(1,2) - F(1,4)*F(2,9)) > 0
# (4.5) junction
res['(4.5) C_new(d-) - C_old(d-) > 0'] = (2*y - (2 - B)*dm)/LAM - min(dm, 2*y - dm)/(6*kap) > 0
# (4.6), (4.7)
res['(4.6) reserve > 1/10'] = theta - (F(7,8) - F(49,150)*h + F(7,8)*ell + F(33,50)*d0) > F(1,10)
res['(4.7) floor reserve > 1/100'] = theta - (F(51,100) + h/3 + F(51,100)*ell) > F(1,100)
# (5.14), (5.15)
res['(5.14a) > 0'] = y/20 - (1 + h)/1000 - F(1,4000) > 0
res['(5.14b) > 0'] = h/600 - F(1,1000) - F(1,4000) > 0
res['(5.15) < -1/20'] = F(13,75)*h - y/2 + F(63,50)*d0 < -F(1,20)
# validity window of the upgraded plain capacity
mlo, mhi = Vrange()
res['upgrade window nonempty'] = mlo < mhi
for k, v in res.items(): print(("OK  " if v else "FAIL"), k)
print("values: n", float(n_lo), float(n_hi), "zA", float(zA(n_hi)), float(zA(n_lo)), "b", float(b), "S_b(t0)", float(Sb(t0)), "window", float(mlo), float(mhi))
print("ALL OK" if all(res.values()) else "SOME FAIL")
