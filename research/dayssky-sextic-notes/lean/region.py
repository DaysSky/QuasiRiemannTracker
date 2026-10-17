from fractions import Fraction as F
import numpy as np
LAM=F(19312,4787); B=F(1419985639,1436100000)
th=F(437423,500000); h=F(795785,1000000)
ell=F(11,3)-4*th; x=1+ell-h; y=h-2*ell; ks=2*th-1
dM=h-F(1,1000); dP=h+F(1,10**10)
print(float(ell),float(x),float(y),float(ks),float(dM))
def Rold(k,d,p):
    den=3*k*(3-p)-2*p; c=2*(1-p)*(3*k-p)/den; rd=(5/6-d)+d*c
    return 1-d+d*c*(5/6-d)/(2*rd)
L=float(LAM); Bf=float(B)
def Rnew(d,p):
    V=5/6-d; DN=3-p-4*p/L; aN=2-4*p/L; cN=1-2*p*Bf/L
    return 1-d+V*d*(1-p)*cN/(V*DN+d*(1-p)*aN)
def tN(d,p):
    V=5/6-d; DN=3-p-4*p/L; aN=2-4*p/L; cN=1-2*p*Bf/L
    c=aN+cN  # c = 3-4p/L-2pB/L ; c-a = cN
    return (V*DN+d*(1-p)*c)/(V*DN+d*(1-p)*aN)
def rN(p,t):
    DN=3-p-4*p/L
    return ((2-4*p/L)*t+2*p*Bf/L-p)/DN
yf,hf,ef=float(y),float(h),float(ell)
def G(d,C,dl,p,R):
    return (1+dl)/2+13*hf/75-yf/2-ef/2+ef*dl*p+d*(1/4+3*R/4-17/50)-dl*p*C/2
newC=lambda d:(2*yf-(2-Bf)*d)/L
oldC=lambda d,k:min(d,2*yf-d)/(6*k)
k=float(ks); T=float(th)
rows=[('new h',hf,newC(hf),1e-6),('new dM',float(dM),newC(float(dM)),1e-5),('old dM',float(dM),oldC(float(dM),k),1e-5),('old 1/2',.5,oldC(.5,k),1e-3)]
ds=np.linspace(0,.75,1501); ps=np.linspace(0,.5,1001)
for name,d,C,v in rows:
    bad=[]
    for dl in ds:
        for p in ps:
            g=G(d,C,dl,p,Rold(k,dl,p))
            if g>T-v: bad.append((dl,p,g-(T-v)))
    if bad:
        a=np.array(bad); print(name,'fails',len(bad),'delta',a[:,0].min(),a[:,0].max(),'p',a[:,1].min(),a[:,1].max(),'worst',a[:,2].max())
    else: print(name,'ok')
