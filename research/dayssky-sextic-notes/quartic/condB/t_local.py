# verify local structure of the full Kubota series D(m^j,u)=sum_{all primary c} G(m^j,c) N(c)^{-u}
# G(r,c)=sum_{x mod c} conj chi_c(x) psi(r x/c).  Claims (m primary prime, c0 sqfree primary coprime to m):
#  G(m^j, c)=0 for non-squarefree c unless c=c0*m^{j+1};
#  j=1: G(m,c0 m^2)  = chi_{c0}(m) G(1,c0) * N(m) g(chi_m^2)
#  j=2: G(m^2,c0 m^3)= s(c0) G(1,c0) * N(m)^2 g(chi_m),  s=(-1)^{((Nc0-1)/4)((Nm-1)/4)}
#  j=3: G(m^3,c0 m^4)= conj chi_{c0}(m) G(1,c0) * (-N(m)^3)
import sys; sys.path.insert(0,'..')
from gauss4 import *
def gpow(x,k):
    r=(1,0)
    for _ in range(k): r=gmul(r,x)
    return r
bad=0; cnt=0
for m in [(-1,-2),(-1,2),(-3,0),(3,2)]:
    Nm=gnorm(m)
    gq=gauss(m,2,(1,0),[m]); g1=gauss(m,1,(1,0),[m])
    for j in (1,2,3):
        r=gpow(m,j)
        c0s=[(1,0)]+[c for c in [(-1,-2),(-1,2),(3,2),(3,-2),(-3,0),(1,4),(1,-4)] if c!=m and gnorm(c)*Nm**(j+1)<=4000]
        for c0 in c0s:
            for k in range(0,j+3):
                c=gmul(c0,gpow(m,k)); Nc=gnorm(c)
                if Nc>4000: continue
                fac=(gfactor(c0) if c0!=(1,0) else [])+[m]*k
                G=gauss(c,3,r,fac)
                if k!=j+1:
                    pred=0 if k>=2 else None
                    if k>=2:
                        cnt+=1
                        if abs(G)>1e-6: bad+=1; print("nonzero unexpected",m,j,c0,k,G)
                    continue
                fac0=gfactor(c0) if c0!=(1,0) else []
                G1=gauss(c0,3,(1,0),fac0) if c0!=(1,0) else 1
                ch=chi(m,c0,fac0) if c0!=(1,0) else 1
                Nc0=gnorm(c0)
                if j==1: pred=ch*G1*Nm*gq
                if j==2: pred=(-1 if (((Nc0-1)//4)%2 and ((Nm-1)//4)%2) else 1)*G1*Nm**2*g1
                if j==3: pred=ch.conjugate()*G1*(-Nm**3) if c0!=(1,0) else G1*(-Nm**3)
                cnt+=1
                if abs(G-pred)>1e-6*abs(pred)+1e-6: bad+=1; print("MISMATCH",m,j,c0,G,pred)
print("checked",cnt,"bad",bad)
