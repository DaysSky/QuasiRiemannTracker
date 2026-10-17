# S51 variant of reduction.py: both X-terms from the G2 series. X_main = Q^{3/4} Z^pi h^beta, X_dual = (QZ)^{1/2+dpi} h^{Bp}.
# pi = 3/8 (convexity at the residue centre); 5/12 if the long-root (3,1) normaliser L(12X1+2, alpha^12) blocks X1 < -1/12;
# 1/2 if the short-root (1,0) normaliser L(4X1+1, alpha^4) blocks X1 < 0 (S51 §20).
import numpy as np
def best(pi,beta,dpi,Bp):
    H=1/3; res=(9,None)
    for k0 in np.linspace(0,H,1501):
        offT=-1
        for eta in np.linspace(0,k0,201):
            xm=0.5+pi+beta*eta
            xd=5/6+dpi+Bp*eta
            offT=max(offT, 1/3+eta+(H-eta)/4+2*max(xm,xd))
        s=max(2+(H-k0)/4,(1+offT)/2,2)
        if s<res[0]: res=(s,k0)
    return res
for pi in (3/8,5/12,1/2):
    for dpi in (0,1/12,1/6):
        s,k0=best(pi,7/8,dpi,1)
        print("pi=%.4f dual+%.3f -> S << Z^%.4f, (L3) << Q^%.4f, K0=Z^%.3f"%(pi,dpi,s,(s-1)*1.5,k0))
