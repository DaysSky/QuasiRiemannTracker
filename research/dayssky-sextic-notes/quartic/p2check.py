# (P2): for coprime squarefree primary c1,c2:  conj chi_{c1}(c2) * chi_{c2}(c1) == (-1)^{((Nc1-1)/4)((Nc2-1)/4)}
# (P3): a(c) * g(conj chi_c) / sqrt(Nc) == mu eps chi_c(-1)   and  conj(a(c)) g(chi_c)/sqrt(Nc) == mu(c) eps(c)
import random, math
from gauss4 import *
random.seed(5)
cs=[(a,b) for a in range(-30,31) for b in range(-30,31) if 1<gnorm((a,b))<=900 and is_primary((a,b))]
cs=[c for c in cs if mobius(gfactor(c))!=0]
bad=0; cnt=0
for _ in range(3000):
    c1,c2=random.sample(cs,2)
    f1,f2=gfactor(c1),gfactor(c2)
    if set(f1)&set(f2): continue
    lhs=chi(c2,c1,f1).conjugate()*chi(c1,c2,f2)
    rhs=(-1)**(((gnorm(c1)-1)//4)*((gnorm(c2)-1)//4))
    cnt+=1; bad+= abs(lhs-rhs)>1e-9
print("P2 pairs",cnt,"mismatches",bad)
bad3=0; n3=0
for c in cs:
    if gnorm(c)>500: continue
    f=gfactor(c); N=gnorm(c); al=complex(*c)/math.sqrt(N)
    eps=1 if (c[0]%4,c[1]%4)==(1,0) else -1
    gm=gauss(c,-1,fac=f); gp=gauss(c,1,fac=f)
    a=al*gm/math.sqrt(N)
    mu=mobius(f); cm1=chi((-1,0),c,f)
    n3+=1
    bad3+= abs(a*gm/math.sqrt(N)-mu*eps*cm1)>1e-8 or abs(a.conjugate()*gp/math.sqrt(N)-mu*eps)>1e-8
print("P3 c checked",n3,"mismatches",bad3)
