# Twist cost in the cone lemma: c-node twist of index-conductor F (t_F = log F). Only sigma changes:
#  (i) rows enlargement (positivity): sigma <= tC * L((tD + tF)/tC)
#  (ii) reflection with twisted conductor: tC' = 3tD + tF - tC, sigma <= tC - tC' + tC' * L((tD+tF)/tC')
#  (iii) negligible if tC > 3tD + tF.  L = limit table (Delta = 0).
import numpy as np
def L(q):
    q=np.asarray(q,float)
    return np.where(q<=1/3,1+q,np.where(q<2/3,2*q+2/3,q+4/3))
def sig_tw(tC,tD,tF):
    out=np.full(tC.shape,-np.inf)
    act=tC<=3*tD+tF+1e-12
    with np.errstate(divide='ignore',invalid='ignore'):
        direct=np.where(tC>0,tC*L((tD+tF)/np.where(tC>0,tC,1)),tD+tF)   # tC=0: single c, sigma<=tD (+tF conservative)
        tCp=3*tD+tF-tC
        refl=np.where((tCp>0)&(tCp<tC),tC-tCp+tCp*L((tD+tF)/np.where(tCp>0,tCp,1)),np.inf)
    v=np.minimum(direct,refl)
    out[act]=v[act]
    return out
def f(x,tK,tC,tD,tF):
    B1=-x*tC+0.5*(np.maximum.reduce([tK+tC,tD,2/3*(tK+tC+tD)])-tK)
    B2=-(0.5+x)*tC+0.5*sig_tw(tC,tD,tF)+0.5*(np.maximum.reduce([np.full(tD.shape,tK),tD,2/3*(tK+tD)])-tK-tD)
    return np.minimum(B1,B2)
x=1/3+1e-4
tC=np.concatenate([np.linspace(0,4,2401),np.linspace(4,40,721)])[:,None]
u=np.linspace(0,1,1201)[None,:]
for tK in (1.0,0.0):
    for tF in (0.0,0.02,0.05,0.1,0.3,1.0):
        if tK==0 and tF==0: continue
        tD=0.5*(tK+tC)*u; TC=np.broadcast_to(tC,tD.shape)
        F=f(x,tK,TC,tD,tF); i=np.unravel_index(np.argmax(F),F.shape)
        b=max(0.0,float(F[i]))
        print("tK=%.0f tF=%.2f: b_loc=%.5f  b/tF=%s  at (tC,tD)=(%.3f,%.3f)"%(tK,tF,b,"-" if tF==0 else "%.4f"%(b/tF),TC[i],tD[i]))
