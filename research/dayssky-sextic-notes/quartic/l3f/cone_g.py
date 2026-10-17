# g(rho) = b_loc(1-rho, rho) at (1/3+dp, 0), theta=1/2, table L+Delta on [2/3,5/6); vectorised fine grid.
import numpy as np
def Lq(q,Delta):
    q=np.asarray(q,float); out=np.empty_like(q)
    lo=(3*q-1)<=1e-12; mid=(~lo)&(q<2/3); hi=q>=2/3
    out[lo]=1+q[lo]
    Yp=3*q[mid]-1; qp=q[mid]/Yp
    out[mid]=1-Yp+Yp*(qp+4/3+Delta*((qp>=2/3)&(qp<5/6)))
    out[hi]=np.minimum(1+np.maximum.reduce([np.ones(hi.sum()),q[hi]+1/3,q[hi]/2+0.75]), q[hi]+4/3+Delta*(q[hi]<5/6))
    return out
def f(x,tK,tC,tD,Delta):
    B1=-x*tC+0.5*(np.maximum.reduce([tK+tC,tD,2/3*(tK+tC+tD)])-tK)
    sig=np.full(tC.shape,-np.inf)
    z=(tD<=0)&(tC<=0); sig[z]=0
    act=(tD>0)&(tC<=3*tD)
    v=tC+2.5*tD
    pos=act&(tC>0)
    r=np.where(pos,tD/np.where(tC>0,tC,1),0)
    v=np.where(pos,np.minimum(v,tC*Lq(r,Delta)),v)
    sig[act]=v[act]
    B2=-(0.5+x)*tC+0.5*sig+0.5*(np.maximum.reduce([np.full(tD.shape,tK),tD,2/3*(tK+tD)])-tK-tD)
    return np.minimum(B1,B2)
def g(rho,x,Delta):
    tK=1-rho; tE=rho
    tC=np.concatenate([np.linspace(0,4,2001),np.linspace(4,60,1121)])[:,None]
    u=np.linspace(0,1,1001)[None,:]
    tD=0.5*(tK+tE+tC)*u; TC=np.broadcast_to(tC,tD.shape)
    F=f(x,tK,TC,tD,Delta); assert not np.isnan(F).any(); return max(0.0,float(F.max()))
x=1/3+1e-4
for Delta in (0.0,0.01,0.03,1/12):
    gs=[g(r,x,Delta) for r in np.linspace(0,1,21)]
    rh=np.linspace(0,1,21)
    c1=gs[0]/Delta if Delta>0 else 0.0
    # smallest c2 with g(rho) <= c1*Delta*(1-rho)+c2*rho for all rho
    c2=max((gg-c1*Delta*(1-r))/r for gg,r in zip(gs[1:],rh[1:]))
    print("Delta=%.4f  g(0)=%.5f (c1=%.3f)  g(1)=%.5f  needed c2=%.4f  g=%s"%(Delta,gs[0],c1,gs[-1],c2,np.round(gs[::4],4)))
