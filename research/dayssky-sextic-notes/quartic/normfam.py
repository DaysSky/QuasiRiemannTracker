# (a) check gtilde(chi_c)^2 = mu(c)*alpha(c) and gtilde_2(c)=chi_c(-1) for squarefree primary c
# (b) norm family psi_c(m) = (N(m)/c)_4 ; check psi_c(m)=chi_c(m)*conj(chi_{cbar}(m)) and that
#     the Gauss sum of psi_c modulo c*cbar (psi(z)=e(Im z)) equals N(c)*mu(c)*(ray-type unit)
import math, cmath, random
from gauss4 import *
random.seed(2)
cs=[(a,b) for a in range(-14,15) for b in range(-14,15) if 1<gnorm((a,b))<=200 and is_primary((a,b))]
bad=0;cnt=0
for n in cs:
    fac=gfactor(n)
    if mobius(fac)==0: continue
    N=gnorm(n); g1=gauss(n,1,(1,0),fac)/math.sqrt(N); g2=gauss(n,2,(1,0),fac)/math.sqrt(N)
    al=complex(*n)/math.sqrt(N)
    cnt+=1
    if abs(g1*g1-mobius(fac)*al)>1e-6 or abs(g2-chi((-1,0),n,fac))>1e-6: bad+=1
print("(a) checked",cnt,"bad",bad)
# (b) norm family on c with (c,cbar)=1, small
def gauss_mod(chi_fn, M, h=(1,0)):
    s=0; Mc=gconj(M); NM=gnorm(M)
    for x in residues(M):
        v=chi_fn(x)
        if v==0: continue
        z=gmul(gmul(h,x),Mc)
        s+=v*cmath.exp(2j*math.pi*z[1]/NM)
    return s
out=[]
for n in cs:
    fac=gfactor(n)
    if mobius(fac)==0: continue
    nb=gconj(n)
    # require no rational prime factors and no conjugate pair: gcd(n,nbar)=1 <=> all prime norms distinct primes p≡1 mod 4 appearing once
    ps=[gnorm(p) for p in fac]
    if any(p%4==3 for p in ps) or len(set(ps))<len(ps): continue
    if gnorm(n)>120: continue
    M=gmul(n,nb)  # = N(n) rational
    psi=lambda x: chi((gnorm(x),0),n,fac)
    # check identity psi(m)=chi_n(m)*conj(chi_nbar(m)) on random m
    facb=gfactor(primary_assoc(nb))
    for _ in range(5):
        m=(random.randint(-50,50),random.randint(-50,50))
        lhs=psi(m); rhs=chi(m,n,fac)*complex(chi(m,primary_assoc(nb),facb)).conjugate()
        if abs(lhs-rhs)>1e-9: print("psi identity fails",n,m,lhs,rhs)
    G=gauss_mod(psi,M)
    Gt=G/gnorm(n)
    al=complex(*n)/math.sqrt(gnorm(n))
    out.append((n,mobius(fac),round(abs(G)/gnorm(n),4),complex(round((Gt/(mobius(fac)*al)).real,4),round((Gt/(mobius(fac)*al)).imag,4)), (n[0]%4,n[1]%4)))
for o in out: print("(b) c=%s mu=%d |g|/N=%s  g~/(mu*alpha)=%s  c mod4=%s"%o)
