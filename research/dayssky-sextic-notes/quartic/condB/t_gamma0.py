# Gamma_0(m) picture: index-m coefficient sums at cusp infinity with nebentypus (d/m)_2:
#   S(m;c) = sum_{d mod c, (d,c)=1} (d/m)_2 * conj((d/c)_4) * psi(m d / c),   c = m^k c''
# claim: nonzero only for k=2, and S(m; m^2 c'') = (-1)^{omega(m)} N(m) G(m^3,c''), G(r,c)=sum conj chi_c(x) psi(r x/c)
import sys, math, cmath; sys.path.insert(0,'..')
from gauss4 import *
def S(m,c,facm,facc):
    N=gnorm(c); s=0; nc=gconj(c)
    for x in residues(c):
        q=chi(x,c,facc)
        if q==0: continue
        qm=chi(x,m,facm)
        if qm==0: continue
        mx=gmul(m,x); z=gmul(mx,nc)
        s+=(qm**2)*q.conjugate()*cmath.exp(2j*math.pi*z[1]/N)
    return s
bad=0;cnt=0
for m in [(-1,-2),(3,2),(-3,0)]:
    for cpp in [(1,0),(-1,2),(3,-2),(1,4)]:
        if cpp==m: continue
        for k in (1,2,3):
            c=cpp
            for _ in range(k): c=gmul(c,m)
            if gnorm(c)>3500: continue
            facc=(gfactor(cpp) if cpp!=(1,0) else [])+[m]*k
            v=S(m,c,[m],facc)
            if k!=2:
                cnt+=1; bad+= abs(v)>1e-6; continue
            G=gauss(cpp,3,gmul(m,gmul(m,m)),gfactor(cpp)) if cpp!=(1,0) else 1
            pred=-gnorm(m)*G
            cnt+=1
            if abs(v-pred)>1e-6*abs(pred): bad+=1; print("mismatch",m,cpp,v,pred)
print("checked",cnt,"bad",bad)
