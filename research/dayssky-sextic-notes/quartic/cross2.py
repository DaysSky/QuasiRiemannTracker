import random, math
from gauss4 import *
cs=[]
for a in range(-60,61):
    for b in range(-60,61):
        n=(a,b); N=gnorm(n)
        if 50<N<=3000 and is_primary(n) and N%2==1:
            fac=gfactor(n)
            if mobius(fac)==0: continue
            ps=[gnorm(p) for p in fac]
            if any(p%4==3 for p in ps) or len(set(ps))<len(ps): continue
            cs.append((n,fac,N))
random.seed(7); res={}
for _ in range(400):
    s=random.choice(cs); c=random.choice(cs)
    sv,sf,sN=s; cv,cf,cN=c
    if math.gcd(sN,cN)!=1: continue
    lhs=complex(chi((cN,0),sv,sf)).conjugate()*chi((sN,0),cv,cf)
    cb=primary_assoc(gconj(cv)); rhs=chi(cb,sv,sf)**2
    r=complex(round((lhs/rhs).real,6),round((lhs/rhs).imag,6))
    key=(r, (sv[0]%4,sv[1]%4),(cv[0]%4,cv[1]%4))
    res[key]=res.get(key,0)+1
for k,v in sorted(res.items(),key=str): print(k,v)
