# Compare the sextic family's current counts with zeta-analogue exponents at the same sigma.
# Normalisation: count exponent relative to family size; zeta N(sigma,T) <= T^{A(1-sigma)} -> A(1-sigma) = A(1-delta)/2.
import sys; sys.path.insert(0,'../c61')
from fastR import R
th=0.8748389; kap=2*th-1
rho=lambda d: 1-10*d/(3*(5-2*d))
print(" delta  sigma  rho(p=0) R(p=1/2)  Ingham  Huxley   GM     DH")
for d in (0.1,0.2,0.3,0.4,0.44,0.466,0.5,0.55,0.6,0.7):
    s=(1+d)/2
    ing=3*(1-d)/(3-d); hux=1.2*(1-d); gm=(15/13)*(1-d)
    print(f" {d:5.3f}  {s:5.3f}  {rho(d):7.4f}  {R(d,0.5,kap,'new'):7.4f}  {ing:6.4f}  {hux:6.4f}  {gm:6.4f}  {1-d:6.4f}")
