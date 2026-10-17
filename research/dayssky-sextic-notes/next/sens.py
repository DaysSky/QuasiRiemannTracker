import numpy as np
from scipy.optimize import minimize_scalar
def run(L=19312/4787,Bc=1419985639/1436100000,mhi=847/2000,invfac=1.0,region=None,cap=True,extra=None):
    m_lo=(Bc-L/25)/2
    dl=np.linspace(0,.75,751)[:,None]; P=np.linspace(0,.5,401)[None,:]
    def Rold(k,d,p):
        den=3*k*(3-p)-2*p; c=2*(1-p)*(3*k-p)/den; rd=(5/6-d)+d*c
        return 1-d+d*c*(5/6-d)/(2*rd)
    V=5/6-dl; DN=3-P-4*P/L; aN=2-4*P/L; cN=1-2*P*Bc/L; J=V*DN+dl*(1-P)*aN
    t=(V*DN+dl*(1-P)*(aN+cN))/J; RN=1-dl+V*dl*(1-P)*cN/J
    r=(aN*t+2*P*Bc/L-P)/DN; m=t-r; sh=1-dl*(P+(1-P)*r)
    Pold=1-2*dl*mhi-2*dl*P*(1-2*mhi)/(6*.75)
    OK=(r>=23/37)&(r<=1)&(m>=m_lo+1e-4)&(Pold<=sh-1e-5)&(t>=1)&(t<=1.5)
    if region is not None:
        a,b,c,d=region; OK=OK&(dl>=a)&(dl<=b)&(P>=c)&(P<=d)
    def margin(th,h):
        ell=11/3-4*th; y=h-2*ell; k=2*th-1
        R=Rold(k,dl,P); R=np.where(OK,np.minimum(R,RN),R)
        out=[]
        for d,C,v in ((h,(2*y-(2-Bc)*h)/L,1e-6),(h-1e-3,min(h-1e-3,2*y-h+1e-3)/(6*k),1e-5),(.5,min(.5,2*y-.5)/(6*k),1e-3)):
            G=(1+dl)/2+13*h/75-y/2-ell/2+ell*dl*P+d*(1/4+3*R/4-17/50)-dl*P*C/2
            out.append(th-v-G.max())
        return min(out)
    def hcap(th):
        ell=11/3-4*th; nm=(Bc-0.037*L)/2
        return min((4*ell-1.16e-3)/0.84,(2*ell-(1-nm)*1e-3)/nm)
    def best(th):
        hi=hcap(th) if cap else .80
        r=minimize_scalar(lambda h:-margin(th,h),bounds=(.78,hi),method='bounded',options={'xatol':1e-9})
        return -r.fun,r.x
    lo,hi=.8740,.87490
    for _ in range(32):
        mid=(lo+hi)/2; mm,h=best(mid)
        if mm>=0: hi=mid; hb=h
        else: lo=mid
    return hi,hb
base=run(region=(89/200,97/200,99/200,1/2))
print('base region',base)
print('base noregion',run())
for L in (4.0,3.9,3.8):
    print('L',L,run(L=L))
for B in (0.99,0.995,1.0):
    print('B',B,run(Bc=B))
