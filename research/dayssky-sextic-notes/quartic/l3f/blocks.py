# Block-optimal mean-value regions for the double sum sum_{c~C, d~D} T*(c) chibar_d(c) gt(h0;chi_d) C^{-1/2-x} D^{-1/2-y}  (fix (a)).
# Block bounds (log, homogeneous in (tC,tD)):
#  B1 = -x tC - y tD + 1/2 max(tC, tD, 2/3 (tC+tD))                          [large sieve over c for the d-sums]
#  B2 = -y tD - (1/2+x) tC + 1/2 sigma(tC,tD)                                  [mean square over d of the c-sums]
#  sigma: S(D,C) = C^{..}: min of convexity tC+2.5tD; S49 direct if tC<=5/3 tD; reflected if tC>=1.5 tD (C'=D^3/C);
#         negligible if tC > 3 tD (length beyond conductor, psi-hat entire).
import numpy as np
L0=13/8
def s49(tCp,tD):   # log of S49 bound (Y + D Y^{1/3} + D^{1/2} Y^{3/4}) Y  with Y=C'
    return tCp+max(tCp, tD+tCp/3, tD/2+0.75*tCp)
def sigma(tC,tD,lam):
    if tC>3*tD+1e-12: return -1e9
    cands=[tC+2.5*tD]
    if tC<=5/3*tD+1e-12: cands.append(s49(tC,tD) if abs(tC-1.5*tD)>1e-9 else tC+lam*tD)
    if tC>=1.5*tD-1e-12:
        tCp=3*tD-tC; v=s49(tCp,tD) if abs(tCp-1.5*tD)>1e-9 else tCp+lam*tD
        cands.append(tC-tCp+v)      # S(D,C) ~ (C/C') S(D,C')
    return min(cands)
def blk(x,y,tC,tD,lam):
    B1=-x*tC-y*tD+0.5*max(tC,tD,2/3*(tC+tD))
    B2=-y*tD-(0.5+x)*tC+0.5*sigma(tC,tD,lam)
    return min(B1,B2)
def region_d_ok(x,y,lam,n=241):     # directions tC=1, tD=delta in [0,1/2]
    return all(blk(x,y,1.0,dl,lam)<-1e-9 for dl in np.linspace(0,0.5,n))
def region_c_ok(x,y,lam,n=241):     # directions tD=1, tC=rho in [0,1.5]; plus tD alone
    return all(blk(x,y,r,1.0,lam)<-1e-9 for r in np.linspace(0,1.5,n))
def b_d(x,y,lam,n=60):             # h0-exponent: max over tC>=0, 0<=tD<=(1+tC)/2 of blk (tK=1)
    best=-1e9
    for tC in np.linspace(0,30,n):
        for tD in np.linspace(0,(1+tC)/2,n):
            best=max(best,blk(x,y,tC,tD,lam))
    return max(best,0.0)
def boundary_d(y,lam):
    lo,hi=-2.0,2.0
    if not region_d_ok(hi,y,lam): return None
    for _ in range(40):
        mid=(lo+hi)/2
        if region_d_ok(mid,y,lam): hi=mid
        else: lo=mid
    return hi
def boundary_c(x,lam):
    lo,hi=-2.0,4.0
    if not region_c_ok(x,hi,lam): return None
    for _ in range(40):
        mid=(lo+hi)/2
        if region_c_ok(x,mid,lam): hi=mid
        else: lo=mid
    return hi
if __name__=="__main__":
    lam=1.5573
    print("Region_d boundary (y>=0): x_min(y), b")
    for y in (0,0.1,0.25,0.5,0.75,1.0,1.5):
        xb=boundary_d(y,lam); print("  y=%.2f x>%.4f  b=%.3f"%(y,xb,b_d(xb+1e-3,y,lam)) if xb is not None else "  y=%.2f none"%y)
    print("Region_c boundary (x>=0): y_min(x)")
    for x in (0,0.05,0.1,0.2,0.3,0.5):
        print("  x=%.2f y>%s"%(x,boundary_c(x,lam)))
