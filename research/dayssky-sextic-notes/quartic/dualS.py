# Check the dual coefficient f(c) = conj(alpha(c)) * gtilde(chi_c)^2  (should be mu(c)*chi_c(-1)*gtilde_2(c))
# and the size of S(H) = sum_c f(c) conj(chi_c(H)) W(Nc/Z) for N(H) in various ranges.
import numpy as np, math, sys, random
from gauss4 import *
random.seed(1)
Z=int(sys.argv[1]) if len(sys.argv)>1 else 3000
cs=[]
R=int(math.isqrt(Z))+2
for a in range(-R,R+1):
    for b in range(-R,R+1):
        n=(a,b); N=gnorm(n)
        if Z/2<N<=Z and is_primary(n) and N%2==1:
            fac=gfactor(n)
            if mobius(fac)==0: continue
            cs.append((n,fac,N))
print("columns",len(cs))
# compute f(c) for a few c directly
for (n,fac,N) in cs[:6]:
    g1=gauss(n,1,(1,0),fac)/math.sqrt(N)
    g2=gauss(n,2,(1,0),fac)/math.sqrt(N)
    al=complex(*n)/math.sqrt(N)
    f=np.conj(al)*g1*g1
    print(n,"mu",mobius(fac),"chi(-1)",chi((-1,0),n,fac),"g2~",np.round(g2,4),"f",np.round(f,4))
