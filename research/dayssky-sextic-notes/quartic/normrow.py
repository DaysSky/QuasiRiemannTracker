# Norm-family B-row: S_n = sum_c conj(alpha(c)) (n/c)_4 W(Nc/Z), rows = rational n (weight r(n)=#{m:N m=n})
# check mean square ~ (#rows)*sum W^2 (perfect orthogonality) for Q below/above Z^{2/3};
# also check the A/B cross factor conj(chi_s(N c))*psi_c(s) vs quadratic symbol (cbar/s)_2.
import numpy as np, math, sys, random
from gauss4 import *
Z=int(sys.argv[1]); Qs=[int(x) for x in sys.argv[2].split(',')]
R=int(math.isqrt(Z))+2
cs=[]
for a in range(-R,R+1):
    for b in range(-R,R+1):
        n=(a,b); N=gnorm(n)
        if Z/2<N<=Z and is_primary(n) and N%2==1:
            fac=gfactor(n)
            if mobius(fac)==0: continue
            cs.append((n,fac,N))
Qmax=max(Qs)
rn=np.zeros(Qmax+1)
for a in range(-int(Qmax**.5)-1,int(Qmax**.5)+2):
    for b in range(-int(Qmax**.5)-1,int(Qmax**.5)+2):
        N=a*a+b*b
        if 0<N<=Qmax: rn[N]+=1
S=np.zeros(Qmax+1,complex); w2=0
# precompute (n/p)_4 for rational n mod p via chi on residues of n mod each prime factor -> use multiplicativity in c
from functools import lru_cache
@lru_cache(None)
def table(pi):
    p=gnorm(pi)
    # rational residues 0..p-1 if p prime ≡1 mod 4, else (inert q, norm q^2): (n/q)_4=1 for q∤n
    if p%4==1:
        return np.array([chi((r,0),pi,[pi]) for r in range(p)],complex),p
    q=int(round(math.sqrt(p)))
    return np.array([0]+[1]*(q-1),complex),q
ns=np.arange(Qmax+1)
for (c,fac,N) in cs:
    tt=N/Z; W=math.sin(math.pi*(2*tt-1))**2 if 0.5<tt<=1 else 0
    if W==0: continue
    w2+=W*W
    v=np.ones(Qmax+1,complex)
    for pi in fac:
        t,mod=table(pi); v*=t[ns%mod]
    al=complex(*c)/math.sqrt(N)
    S+=W*np.conj(al)*v
print(f"Z={Z} cols={len(cs)} Z^(2/3)={Z**(2/3):.0f}")
for Q in Qs:
    tot=np.sum(rn[1:Q+1]*np.abs(S[1:Q+1])**2); nrows=np.sum(rn[1:Q+1])
    print(f" Q={Q:6d} rows(m)={int(nrows):7d}  ratio={tot/(nrows*w2):.3f}")
# cross factor test
random.seed(5); bad=0; tot=0
prim=[c for c in cs[:40]]
for _ in range(60):
    s=random.choice(prim); c=random.choice(prim)
    if s[0]==c[0]: continue
    sv,sf,sN=s; cv,cf,cN=c
    lhs=complex(chi((cN,0),sv,sf)).conjugate()*chi((sN,0),cv,cf)
    cb=primary_assoc(gconj(cv))
    rhs=chi(cb,sv,sf)**2
    tot+=1
    if abs(lhs-rhs)>1e-9 and abs(lhs+rhs)>1e-9: bad+=1
print("cross factor = +-(cbar/s)_2 :", tot-bad,"of",tot)
