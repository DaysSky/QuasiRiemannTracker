import sys; sys.path.insert(0,'.')
from gauss4 import *
import cmath, math
rows=[]
for a in range(-40,41):
    for b in range(-40,41):
        n=(a,b)
        if gnorm(n)<5 or gnorm(n)>1300 or not is_primary(n): continue
        fac=gfactor(n)
        mu=mobius(fac)
        if mu==0: continue
        N=gnorm(n)
        gm1=gauss(n,-1,fac=fac)/math.sqrt(N)
        g1=gauss(n,1,fac=fac)/math.sqrt(N)
        alpha=complex(*n)/math.sqrt(N)
        eps=mu*gm1/(alpha.conjugate()*g1)
        rows.append((n,len(fac),eps,(a%8,b%8)))
# check eps is +-1 and depends only on n mod 4 (or 8)
from collections import defaultdict
cls=defaultdict(set)
bad=0
for n,k,eps,r in rows:
    e=complex(round(eps.real,6),round(eps.imag,6))
    if abs(abs(eps)-1)>1e-6: bad+=1
    cls[(n[0]%4,n[1]%4)].add(e)
print("count",len(rows),"nonunit",bad)
for k,v in sorted(cls.items()): print(k,v)
