# Search for the local action matching the true Q(i) slice signs: CG action with RAW Gauss sums and Q_{beta,i} multiplied by
# s^{e(beta,i)}; test (a) the averaged p-part is a polynomial, (b) Theorem 4.1 for the true slice factor f' (signs derived in rank 1).
import sympy as sp, sys, itertools
from ppart import ppart, x1, x2, C, L2
s=int(sys.argv[2]); r=sp.Integer(int(sys.argv[3])); q=r*r; a=sp.Rational(sys.argv[4])
LAMS=[tuple(map(int,z.split(','))) for z in sys.argv[5].split(';')]
graw={0:sp.Integer(-1),1:a,2:1/r,3:sp.Integer(s)/(q*a)}
x=sp.symbols('x')
def B(b,i): return b[i-1]+sum(C[(j,i)]*b[j-1] for j in (1,2) if j>i)
cands={
 'muB'     : lambda b,i,mu: s**((mu*B(b,i))%2),
 'muB+mu'  : lambda b,i,mu: s**((mu*(B(b,i)+1))%2),
 'BD'      : lambda b,i,mu: s**((((-L2[i]*mu)%4)*(((-L2[i]*mu)%4)-1)//2)%2),
}
which=sys.argv[1]
sg=cands[which]
tot=bad=0; nonpoly=0
for lam in LAMS:
    Np,den=ppart(q,graw,lam,sg)
    if den!=1: nonpoly+=1; print(which,"lam",lam,"NOT polynomial"); continue
    coef=Np.as_dict()
    for i in (1,2):
        j=3-i; t=L2[i]
        for bj in sorted(set(e[j-1] for e in coef)):
            l=lam[i-1]+(-C[(j,i)])*bj
            Nk={k:sum(v*x**e[i-1] for e,v in coef.items() if e[j-1]==bj and e[i-1]%4==k) for k in range(4)}
            K=lambda k:(l-2*k)%4
            for k in range(4):
                kp=(l+1-k)%4; Kp=K(kp); f=Nk[k]
                if Kp!=3:
                    S=sp.Integer(s)**((Kp+1)*(kp+(C[(j,i)]*bj if j>i else 0)))
                    f=f-q**Kp*q*graw[(t*(Kp+1))%4]*x**(Kp+1)*S*Nk[kp]
                f=sp.cancel(f/(1-q**3*x**4))
                if f==0: continue
                tot+=1
                if sp.cancel(f/((q*x)**(l-K(k))*f.subs(x,1/(q**2*x))))!=1: bad+=1
print("cand=%s s=%+d r=%s a=%s lams=%s: nonpoly=%d, Thm4.1 classes tested %d, fail %d"%(which,s,r,a,LAMS,nonpoly,tot,bad))
