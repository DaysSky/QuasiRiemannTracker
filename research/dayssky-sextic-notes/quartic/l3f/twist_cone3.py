# Twisted family (index-conductor F) with its OWN table (no full-rows positivity):
#  q>=2/3 direct: tC*L(q) + kap*tF                         (twisted bookkeeping excess <= c*t/5.5, kap=c/5.5)
#  q<2/3: min( Cauchy floor 2tC + max(0,c*tF - 2.75*v*tC)/5.5  [v=2/3-q],
#              twisted reflection tC-tC'+tC'*L(q')+kap*tF if tC'=3tD+tF-tC in (0,tC), q'=tD/tC'>=2/3 )
#  negligible if tC>3tD+tF; tC=0: sigma<=tD.
import numpy as np
def L(q):
    q=np.asarray(q,float); return np.where(q<=1/3,1+q,np.where(q<2/3,2*q+2/3,q+4/3))
def sig(tC,tD,tF,c):
    kap=c/5.5
    with np.errstate(divide='ignore',invalid='ignore'):
        s=np.where(tC>0,tC,1.0); q=tD/s
        direct=np.where((tC>0)&(q>=2/3),tC*L(q)+kap*tF,np.inf)
        v=2/3-q
        floor=np.where((tC>0)&(q<2/3),2*tC+np.maximum(0,c*tF-2.75*v*tC)/5.5,np.inf)
        tCp=3*tD+tF-tC; tCp=np.where(np.abs(tCp)<1e-9,0.0,tCp); sp=np.where(tCp>0,tCp,1.0); qp=tD/sp
        reflval=np.where(tCp>0,tC-tCp+tCp*L(qp),tC+tD)+kap*tF          # tC'=0: single dual term, sigma'=tD
        refl=np.where((tC>0)&(tCp>=0)&(tCp<tC)&((qp>=2/3)|(tCp<=0)),reflval,np.inf)
    out=np.minimum(np.minimum(direct,floor),refl)
    out=np.where(tC==0,tD,out)
    return np.where(tC>3*tD+tF+1e-9,-np.inf,out)
def bloc(x,tK,tE,tF,c):
    tC=np.concatenate([np.linspace(0,4,2401),np.linspace(4,40,721)])[:,None]
    u=np.linspace(0,1,1201)[None,:]
    tD=0.5*(tK+tE+tC)*u; TC=np.broadcast_to(tC,tD.shape)
    B1=-x*TC+0.5*(np.maximum.reduce([tK+TC,tD,2/3*(tK+TC+tD)])-tK)
    B2=-(0.5+x)*TC+0.5*sig(TC,tD,tF,c)+0.5*(np.maximum.reduce([np.full(tD.shape,tK),tD,2/3*(tK+tD)])-tK-tD)
    F=np.minimum(B1,B2); i=np.unravel_index(np.argmax(F),F.shape)
    return max(0.0,float(F[i])),(float(TC[i]),float(tD[i]))
x=1/3+1e-4
c=7/12
for it in range(10):
    vals=[]
    for tK in (0.0,0.3,1.0):
        for tE in (0.0,0.2):
            for tF in (0.05,0.3,1.0):
                b,arg=bloc(x,tK,tE,tF,c); vals.append(((b-0.5*tE)/tF,tK,tE,tF,arg))
    c3,tK,tE,tF,arg=max(vals)
    cn=7/12+c3
    print("iter %d: c=%.4f -> c3=%.4f (worst at tK=%.1f tE=%.1f tF=%.2f block %s) -> c'=%.4f   [boxes need c<2.75]"%(it,c,c3,tK,tE,tF,np.round(arg,3),cn),flush=True)
    if abs(cn-c)<1e-3: break
    c=cn
