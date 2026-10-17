import sys, math, cmath
sys.path.insert(0,'..')
from gauss4 import gauss, gfactor, is_primary, gnorm, chi, mobius, gmul
import j4test as J
# rebuild prime data like j4test for Z=200 and compare T(c) with brute force alpha(c)*g(conj chi_c)/sqrt(N)
def T_brute(c):
    f=gfactor(c); N=gnorm(c)
    return complex(*c)/math.sqrt(N)*gauss(c,-1,(1,0),f)/math.sqrt(N)
import numpy as np
# replicate prime Gauss sum computation
def prime_T(pi,p):
    a,b=pi; r=(-a*pow(b,-1,p))%p
    g=2
    def isprim(g):
        n=p-1; f=2; fs=set()
        while f*f<=n:
            while n%f==0: fs.add(f); n//=f
            f+=1
        if n>1: fs.add(n)
        return all(pow(g,(p-1)//q,p)!=1 for q in fs)
    while not isprim(g): g+=1
    lg={}; y=1
    for s in range(p-1): lg[y]=s; y=y*g%p
    z=pow(g,(p-1)//4,p); eps=1 if z==r else 3
    gs=sum(cmath.exp(-0.5j*math.pi*((eps*lg[x])%4))*cmath.exp(2j*math.pi*((-x*b)%p)/p) for x in range(1,p))
    return complex(*pi)/math.sqrt(p)*gs/math.sqrt(p)
for pi in [(1,4),(-1,2),(1,-4),(-1,-2),(5,4),(5,-4),(-5,2),(-5,-2),(1,10),(1,-10)]:
    if not is_primary(pi): print("not primary",pi); continue
    p=gnorm(pi)
    print(pi,p, np.round(prime_T(pi,p),6), np.round(T_brute(pi),6))
print("composite check")
import itertools
ps=[((1,4),17),((-1,2),5),((1,-4),17),((-1,-2),5),((5,4),41),((5,-4),41),((-5,2),29),((-5,-2),29)]
bad=0;n=0
for (p1,n1),(p2,n2) in itertools.combinations(ps,2):
    for order in [((p1,n1),(p2,n2))]:
        (a,na),(b,nb)=order
        c=gmul(a,b)
        Ta,Tb=prime_T(a,na),prime_T(b,nb)
        q=chi(a,b,[b])**2
        s=-1 if (na%8==5 and nb%8==5) else 1
        Tc=Ta*Tb*q*s
        tb=T_brute(c); n+=1
        if abs(Tc-tb)>1e-6: bad+=1; print("mismatch",a,b,Tc,tb)
print("pairs",n,"bad",bad)
