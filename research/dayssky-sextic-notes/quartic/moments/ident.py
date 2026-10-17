# Verify the two "pairing" identities that decide the structure of the w-functional equation:
#  cube index  (psi(m^3,s)=sum_c a_c conj chi_c(m) N c^{-s}, a_c = alpha(c) gt(conj chi_c)):
#        root number of L(w, conj chi_c) ~ gt(conj chi_c)  ->  a_c * gt(conj chi_c) = mu(c)          [Mobius => 1/L]
#        conj(a_c) * gt(chi_c) = mu(c) chi_c(-1)                                                    [Mobius => 1/L]
#  first-power index (psi(m,s)=sum_c alpha(c) gt(conj chi_c, m) ...= sum_c alpha(c) gt(conj chi_c) chi_c(m)):
#        b_c = alpha(c) gt(conj chi_c) paired with L(w, chi_c), root number ~ gt(chi_c):
#        b_c * gt(chi_c) = alpha(c) chi_c(-1)                                                       [Hecke char => L]
# plus the CRT cross factor for g(conj chi_c1 * chi_c2).
import sys, math, cmath
sys.path.insert(0,'..')
from gauss4 import *
rows=[]
for a in range(-30,31):
    for b in range(-30,31):
        n=(a,b); N=gnorm(n)
        if N<5 or N>700 or not is_primary(n): continue
        fac=gfactor(n); mu=mobius(fac)
        if mu==0: continue
        g1=gauss(n,1,fac=fac)/math.sqrt(N); gm=gauss(n,-1,fac=fac)/math.sqrt(N)
        al=complex(*n)/math.sqrt(N); cm1=chi((-1,0),n,fac)
        A=al*gm
        e1=abs(A*gm-mu); e2=abs(A.conjugate()*g1-mu*cm1); e3=abs(A*g1-al*cm1)
        rows.append((n,mu,e1,e2,e3))
print("squarefree primary c, N<=700:",len(rows))
print("max |a_c gt(conj chi_c) - mu(c)|          =",max(r[2] for r in rows))
print("max |conj(a_c) gt(chi_c) - mu(c)chi_c(-1)| =",max(r[3] for r in rows))
print("max |a_c gt(chi_c) - alpha(c)chi_c(-1)|   =",max(r[4] for r in rows))
# CRT cross factor: g(conj chi_c1 chi_c2) = conj chi_c1(c2) chi_c2(c1) g(conj chi_c1) g(chi_c2)
import random
random.seed(1)
pr=[r[0] for r in rows if r[1]==-1 and gnorm(r[0])<120]
bad=0;cnt=0
for c1 in pr:
    for c2 in pr:
        if c1==c2 or gnorm(c1)==gnorm(c2) and False: pass
        if c1==c2: continue
        if gnorm(c1)*gnorm(c2)>1500: continue
        n=gmul(c1,c2); N=gnorm(n); f1=[c1]; f2=[c2]
        s=0
        nc=gconj(n)
        for x in residues(n):
            u=chi(x,c1,f1); v=chi(x,c2,f2)
            if u==0 or v==0: continue
            z=gmul(x,nc)
            s+=u.conjugate()*v*cmath.exp(2j*math.pi*z[1]/N)
        G1=gauss(c1,-1,fac=f1); G2=gauss(c2,1,fac=f2)
        pred=chi(c2,c1,f1).conjugate()*chi(c1,c2,f2)*G1*G2
        cnt+=1
        if abs(s-pred)>1e-6*abs(pred): bad+=1
        # reciprocity sign
print("CRT cross-factor check: pairs",cnt,"failures",bad)
# reciprocity: conj chi_c1(c2) chi_c2(c1) = (-1)^{((N1-1)/4)((N2-1)/4)} for primary primes
bad=0;cnt=0
for c1 in pr:
    for c2 in pr:
        if c1==c2: continue
        x=chi(c2,c1,[c1]).conjugate()*chi(c1,c2,[c2])
        sgn=(-1)**(((gnorm(c1)-1)//4)*((gnorm(c2)-1)//4))
        cnt+=1; bad+= abs(x-sgn)>1e-9
print("quartic reciprocity sign check:",cnt,"pairs, failures",bad)
