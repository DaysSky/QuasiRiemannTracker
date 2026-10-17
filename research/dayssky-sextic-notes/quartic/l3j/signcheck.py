# Check the THEOREM_ZF Step 2 identity on coprime pairs of squarefree primary c1 != c2:
#   T(c1) conj(T(c2)) g(k; conj(chi_c1) chi_c2) = (-1)^{e1 e2 + e2} mu(c1)mu(c2) chi_c1(k) conj(chi_c2(k)) sqrt(N c1 c2),
# where T(c) = alpha(c) gtilde(conj chi_c), e = (N(c)-1)/4 mod 2.  Brute force over residues mod c1 c2.
import sys, math, cmath, random
sys.path.insert(0,'..')
from gauss4 import *
random.seed(7)
def T(c,fac):
    return complex(*c)/math.sqrt(gnorm(c)) * gauss(c,-1,(1,0),fac)/math.sqrt(gnorm(c))
cs=[(a,b) for a in range(-12,13) for b in range(-12,13) if 4<gnorm((a,b))<=90 and is_primary((a,b))]
cs=[c for c in cs if mobius(gfactor(c))!=0]
bad=0; cnt=0; signs={}
for trial in range(60):
    c1,c2=random.sample(cs,2)
    f1,f2=gfactor(c1),gfactor(c2)
    if set(f1)&set(f2): continue
    q=gmul(c1,c2); fq=f1+f2; Nq=gnorm(q); qc=gconj(q)
    R=residues(q)
    vals=[(x, chi(x,c1,f1), chi(x,c2,f2)) for x in R]
    t1,t2=T(c1,f1),T(c2,f2)
    e1=((gnorm(c1)-1)//4)%2; e2=((gnorm(c2)-1)//4)%2
    for _ in range(3):
        k=(random.randint(-20,20),random.randint(-20,20))
        s=0
        for x,a1,a2 in vals:
            if a1==0 or a2==0: continue
            z=gmul(gmul(k,x),qc)
            s+=a1.conjugate()*a2*cmath.exp(2j*math.pi*z[1]/Nq)
        lhs=t1*t2.conjugate()*s
        base=mobius(f1)*mobius(f2)*chi(k,c1,f1)*chi(k,c2,f2).conjugate()*math.sqrt(Nq)
        pred=(-1)**(e1*e2+e2)*base
        cnt+=1
        if abs(base)<1e-9:
            if abs(lhs)>1e-6: bad+=1
            continue
        r=lhs/base
        signs[(e1,e2)]=signs.get((e1,e2),set())|{round(r.real,6)+1j*round(r.imag,6)}
        if abs(lhs-pred)>1e-6*max(1,abs(lhs)): bad+=1
print("checked",cnt,"mismatches",bad)
print("observed ratio lhs/base by (e1,e2):",{k:sorted(v,key=lambda z:(z.real,z.imag)) for k,v in signs.items()})
