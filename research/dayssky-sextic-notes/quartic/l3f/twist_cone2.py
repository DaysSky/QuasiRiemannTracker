# Self-consistent twist cost. Twisted table hypothesis H_kappa (index-conductor F, tF=log F):
#   direct (q=tD/tC >= 2/3):            sigma <= tC*L(q) + kappa*tF
#   twisted reflection (tC' = 3tD+tF-tC, q'=tD/tC' >= 2/3, 0<tC'<tC):  sigma <= tC - tC' + tC'*L(q') + kappa*tF
#   gap region (neither applies): rows-enlargement by positivity, sigma <= tC*L((tD+tF)/tC)   [crude, width ~tF]
#   negligible: tC > 3tD + tF.   tC = 0: sigma <= tD (single c).
import numpy as np
def L(q):
    q=np.asarray(q,float); return np.where(q<=1/3,1+q,np.where(q<2/3,2*q+2/3,q+4/3))
def sig(tC,tD,tF,kap):
    with np.errstate(divide='ignore',invalid='ignore'):
        safe=np.where(tC>0,tC,1.0)
        q=tD/safe
        direct=np.where((tC>0)&(q>=2/3),tC*L(q)+kap*tF,np.inf)
        tCp=3*tD+tF-tC; sp=np.where(tCp>0,tCp,1.0); qp=tD/sp
        refl=np.where((tCp>0)&(tCp<tC)&(qp>=2/3),tC-tCp+tCp*L(qp)+kap*tF,np.inf)
        crude=np.where(tC>0,tC*L((tD+tF)/safe),tD)
    v=np.minimum(np.minimum(direct,refl),crude)
    v=np.where(tC==0,np.minimum(v,tD),v)
    return np.where(tC>3*tD+tF+1e-12,-np.inf,v)
def bloc(x,tK,tF,kap):
    tC=np.concatenate([np.linspace(0,4,2401),np.linspace(4,40,721)])[:,None]
    u=np.linspace(0,1,1201)[None,:]
    tD=0.5*(tK+tC)*u; TC=np.broadcast_to(tC,tD.shape)
    B1=-x*TC+0.5*(np.maximum.reduce([tK+TC,tD,2/3*(tK+TC+tD)])-tK)
    B2=-(0.5+x)*TC+0.5*sig(TC,tD,tF,kap)+0.5*(np.maximum.reduce([np.full(tD.shape,tK),tD,2/3*(tK+tD)])-tK-tD)
    F=np.minimum(B1,B2); i=np.unravel_index(np.argmax(F),F.shape)
    return max(0.0,float(F[i])),(float(TC[i]),float(tD[i]))
x=1/3+1e-4
kap=0.0
for it in range(8):
    c3=max(bloc(x,tK,tF,kap)[0]/tF for tK in (0.0,0.3,1.0) for tF in (0.05,0.3,1.0))
    c=7/12+c3; kn=c/5.5
    print("iter %d: kappa=%.4f -> c3=%.4f, c=7/12+c3=%.4f (need <= 11/12=0.9167), next kappa=c/5.5=%.4f"%(it,kap,c3,c,kn),flush=True)
    if abs(kn-kap)<1e-4: break
    kap=kn
for tK in (0.0,1.0):
    for tF in (0.05,1.0):
        b,arg=bloc(x,tK,tF,kap); print("   tK=%.0f tF=%.2f: b/tF=%.4f at (tC,tD)=(%.3f,%.3f)"%(tK,tF,b/tF,*arg))
