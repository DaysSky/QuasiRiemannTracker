# For all l in {0..7}^2: max over nonconstant monomials of e(k) - k1 - k2 (normalised size at (x',y') = (0,0), the corner of the
# quadrant x', y' >= 0 where all region points used by S52_LIMIT_PROOF lie), and for l = 0 the non-first-order convergence margin.
import math, itertools
from fractions import Fraction as Fr
import thm41sym as T
def exps(s,l,r0):
    N,_=T.ppart(s,l,'slice'); res={}
    nv={}; 
    for e,c in N.n.to_dict().items():
        key=(int(e[2]),int(e[3])); nv[key]=nv.get(key,Fr(0))+Fr(int(c))*Fr(r0)**int(e[0])*Fr(1,r0)**int(e[1])
    dv=sum(Fr(int(c))*Fr(r0)**int(e[0])*Fr(1,r0)**int(e[1]) for e,c in N.d.to_dict().items())
    for k,v in nv.items():
        if v!=0:
            w=abs(v/dv); res[k]=(math.log(w.numerator)-math.log(w.denominator))/math.log(r0)/2
    return res
worst=-9; arg=None
for s in (-1,1):
  for l in itertools.product(range(8),range(8)):
    E1=exps(s,l,10007); E2=exps(s,l,10009)
    assert set(E1)==set(E2) and all(abs(E1[k]-E2[k])<2e-3 for k in E1), l
    for k,e in E1.items():
        if k==(0,0): continue
        v=e-k[0]-k[1]
        if v>worst+1e-9: worst,arg=v,(s,l,k,round(e,3))
print("max over s=±1, l in {0..7}^2, nonconstant monomials of e(k)-k1-k2 =",round(worst,4),"at",arg)
