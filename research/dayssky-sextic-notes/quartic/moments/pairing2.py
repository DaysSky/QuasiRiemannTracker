# j=2 (square index): root number of L(w, chi_c^2) ~ gt(chi_c^2) (quadratic) -- check it is a unit (+-1), so a_c*gt(chi_c^2)=+-a_c: Gauss sums persist.
# cocycle: a_{c1c2} = s(c1,c2) (c1/c2)_2 a_c1 a_c2 with s=+-1 -> psi(m^3,s) has no Euler product.
import sys, math
sys.path.insert(0,'..')
from gauss4 import *
P=[]
for a in range(-25,26):
    for b in range(-25,26):
        n=(a,b); N=gnorm(n)
        if 5<=N<=400 and is_primary(n) and len(gfactor(n))==1: P.append(n)
vals=set()
for p in P:
    g2=gauss(p,2,fac=[p])/math.sqrt(gnorm(p)); vals.add((round(g2.real,6),round(g2.imag,6)))
print("primes:",len(P)," values of gt(chi_pi^2):",sorted(vals))
def a(c,fac):
    N=gnorm(c); return complex(*c)/math.sqrt(N)*gauss(c,-1,fac=fac)/math.sqrt(N)
cnt=0; bad=0; signs={}
small=[p for p in P if gnorm(p)<=60]
for i,p in enumerate(small):
    for q in small[i+1:]:
        if gnorm(p)*gnorm(q)>1300: continue
        c=gmul(p,q); r=a(c,[p,q])/(a(p,[p])*a(q,[q]))
        q2=chi(p,q,[q])**2      # (p/q)_2
        s=r/q2; cnt+=1
        if abs(abs(s)-1)>1e-8 or min(abs(s-1),abs(s+1))>1e-8: bad+=1
        key=((gnorm(p)%16),(gnorm(q)%16)); signs.setdefault(key,set()).add(round(s.real))
print("cocycle pairs:",cnt,"non-(+-1) ratios:",bad)
print("sign s(c1,c2) by (N c1 mod 16, N c2 mod 16):",{k:sorted(v) for k,v in sorted(signs.items())})
