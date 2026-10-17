# verify kernel symbols against gauss4.chi on random (k,c)
import sys; sys.path.insert(0,'..'); sys.path.insert(0,'.')
import numpy as np, random
from pocore import *
from gauss4 import chi, gmul, primary_assoc, gnorm, chi_prime
D=3000
kinds,mods,rs,offs,norms,tab=build_primes(D)
# reconstruct generator for each prime ideal: split: find a+bi with a+b*r = 0 mod p
gens=[]
for j in range(len(kinds)):
    p=int(mods[j])
    if kinds[j]==1: gens.append(primary_assoc((p,0))); continue
    r=int(rs[j]); found=None
    for a in range(1,int(p**.5)+1):
        b2=p-a*a; b=int(round(b2**.5))
        if b*b==b2:
            for (A,B) in ((a,b),(a,-b),(b,a),(b,-a)):
                if (A+B*r)%p==0: found=(A,B)
            if found: break
    gens.append(found)
ptr,idx,cn,om=build_columns(norms,D)
random.seed(1)
cols=random.sample(range(len(cn)),60)
xs=np.array([random.randint(-300,300) for _ in range(40)],np.int64); ys=np.array([random.randint(-300,300) for _ in range(40)],np.int64)
C=len(cn); V=1
bad=0;tot=0
for c in cols:
    cre=np.zeros((C,1)); cim=np.zeros((C,1)); cre[c,0]=1.0
    Fre,Fim,dg=kernel(xs,ys,kinds,mods,rs,offs,tab,ptr,idx,cre,cim,np.ones(C))
    n=(1,0)
    for t in range(ptr[c],ptr[c+1]): n=gmul(n,gens[idx[t]])
    for i in range(len(xs)):
        k=(int(xs[i]),int(ys[i]))
        if k==(0,0): continue
        val=1
        for t in range(ptr[c],ptr[c+1]): val*=chi_prime(k,gens[idx[t]])
        got=Fre[i,0]+1j*Fim[i,0]
        tot+=1
        if abs(got-val)>1e-9: bad+=1
print("checked",tot,"bad",bad, "primes",len(kinds),"cols",C)
