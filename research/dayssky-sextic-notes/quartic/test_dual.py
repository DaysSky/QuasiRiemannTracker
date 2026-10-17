from gauss4 import *
import math, random
random.seed(3); bad=0; cnt=0
ns=[(a,b) for a in range(-25,26) for b in range(-25,26) if 5<gnorm((a,b))<=500 and is_primary((a,b))]
for n in ns:
    fac=gfactor(n); mu=mobius(fac)
    if mu==0: continue
    eps=1 if (n[0]%4,n[1]%4)==(1,0) else -1
    al=complex(*n)/math.sqrt(gnorm(n))
    for _ in range(3):
        h=(random.randint(-30,30),random.randint(-30,30))
        if h==(0,0): continue
        lhs=mu*gauss(n,1,h,fac)
        q=chi(h,n,fac)**2
        rhs=eps*al*q*gauss(n,-1,h,fac)
        cnt+=1
        if abs(lhs-rhs)>1e-6*max(1,abs(lhs)): bad+=1
print("checked",cnt,"mismatches",bad)
