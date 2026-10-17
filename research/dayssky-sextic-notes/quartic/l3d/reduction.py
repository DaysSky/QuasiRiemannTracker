# Exponent bookkeeping (base Z; Q = Z^{2/3}, H = Q^2/Z = Z^{1/3}) for the Cauchy route on the main box.
#  S^2 <= Z (Z S + Off_T + Off_gen);  Off_gen = Z (H/K0)^{1/4} S  (generic large sieve for twists h0 > K0, j^4-multiplicity)
#  Off_T = (Z/Q) sum_{h0<=K0} (H/h0)^{1/4} |X(h0)|^2,  X = X_main + X_dual
#  X_main = Q^{3/4} * |sum_c T(c) rho(h0 c)|  with  |sum_c T rho| <= Z^pi h^beta   (input P4; Cauchy/Rankin-Selberg gives pi = 7/8, beta = 3/8)
#  X_dual <= Z h^{1/2} (h <= Z^{1/15}),  Z^{17/18} h^{4/3} (h > Z^{1/15})   (double reflection + large sieve on reflected rows)
import numpy as np
def best(pi,beta):
    H=1/3; res=(9,None)
    for k0 in np.linspace(0,H,3001):
        offT=-1
        for eta in np.linspace(0,k0,301):
            xm=0.5+pi+beta*eta
            xd=1+eta/2 if eta<=1/15 else 17/18+4*eta/3
            x=max(xm,xd)
            offT=max(offT, 1/3 + eta + (H-eta)/4 + 2*x)
        s=max(2+(H-k0)/4, (1+offT)/2, 2)
        if s<res[0]: res=(s,k0)
    return res
for pi,beta,lab in ((7/8,3/8,"Cauchy (RS mean square of theta)"),(0.8,0.375,""),(0.75,0.25,"one-sided (Q2 trick)"),(0.6,0.375,""),(0.5,0.375,""),(3/8,0.375,"square-root (P4)"),(3/8,0.5,""),(3/8,1.0,"")):
    s,k0=best(pi,beta)
    print("pi=%.3f beta=%.3f -> S << Z^%.4f  (L3 sum << Q^%.4f) K0=Z^%.3f %s"%(pi,beta,s,(s-1)*1.5,k0,lab))
