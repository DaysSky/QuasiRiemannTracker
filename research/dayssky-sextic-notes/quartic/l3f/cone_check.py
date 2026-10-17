# Independent check (fresh code) of the limit mechanism at the Region_d boundary point (1/3+dp, 0), theta=1/2, non-polar.
# Table: L(q)+Delta on [2/3,5/6), L(q)=q+4/3 (=s49 for q>=5/6); reflection below 2/3; 1+q below 1/3.
import numpy as np
def s49(q): return 1+max(1.0,q+1/3,q/2+0.75)
def Stab(q,Delta):
    if q<=1/3: return 1+q
    if q<2/3-1e-12:
        Yp=3*q-1; return 1-Yp+Yp*Stab(q/Yp,Delta)
    return min(s49(q), q+4/3+(Delta if q<5/6 else 0.0))
def sigma(tC,tD,Delta):
    if tD<=0: return 0.0 if tC<=0 else -np.inf
    if tC>3*tD: return -np.inf
    v=tC+2.5*tD
    if tC>0: v=min(v,tC*Stab(tD/tC,Delta))
    return v
def f(x,y,tK,tC,tD,Delta):
    B1=-x*tC-y*tD+0.5*(max(tK+tC,tD,2/3*(tK+tC+tD))-tK)
    B2=-y*tD-(0.5+x)*tC+0.5*sigma(tC,tD,Delta)+0.5*(max(tK,tD,2/3*(tK+tD))-tK-tD)
    return min(B1,B2)
def cone_max(x,Delta,n=1201):
    # slice tK=1 (homogeneity), tC in [0,40], tD in [0,(1+tC)/2]
    best=(-np.inf,None)
    for tC in np.concatenate([np.linspace(0,4,n),np.linspace(4,40,361)]):
        for u in np.linspace(0,1,201):
            tD=0.5*(1+tC)*u; v=f(x,0.0,1.0,tC,tD,Delta)
            if v>best[0]: best=(v,(tC,tD))
    return best
def asym_ok(x,Delta):
    return max(f(x,0.0,0.0,1.0,d,Delta) for d in np.linspace(0,0.5,2001))
for dp in (1e-2,1e-3,1e-4):
    x=1/3+dp
    print("x=1/3+%g: asymptotic max (needs <0): %.3e"%(dp,asym_ok(x,0.0)))
    for Delta in (0.0,0.002,0.01):
        v,arg=cone_max(x,Delta)
        print("   Delta=%.3f: max_{cone, tK=1} min(B1,B2) = %.5f at (tC,tD)=%s  -> b_loc(eta)=%.5f*eta"%(Delta,v,np.round(arg,4),max(v,0)))
