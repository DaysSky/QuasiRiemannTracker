# SYMBOLIC proof of the group relations for the slice-signed action: q = r^2 and gamma(1) = a are indeterminates,
# gamma(2) = 1/r, gamma(3) = s/(q a), gamma(0) = -1.  Exact arithmetic in the rational function field Q(r,a,x1,x2).
# Relations are checked as identities of class-vectors of rational functions:  1|s_i s_i = 1, 1|(s1 s2)^6 = 1 on x^beta, all 16 classes.
import sys, time
from sympy import QQ
from sympy.polys.fields import field
K,r,a,x1,x2=field("r,a,x1,x2",QQ)
n=4; L2={1:1,2:3}; Cm={(1,1):2,(2,2):2,(1,2):-1,(2,1):-3}; m={1:4,2:4}
pair=lambda b,i: b[0]*Cm[(1,i)]+b[1]*Cm[(2,i)]
def B(b,i): return b[i-1]+sum(Cm[(j,i)]*b[j-1] for j in (1,2) if j>i)
def ev(f,img):
    # evaluate a field element at (r,a,x1,x2) -> img (field elements), as a homomorphism
    def evp(p):
        out=K(0)
        for mon,c in p.terms():
            t=K(c)
            for base,e in zip(img,mon):
                if e: t*=base**e
            out+=t
        return out
    return evp(f.numer)/evp(f.denom)
def check(s,l,mode):
    q=r*r
    gam={0:K(-1),1:a,2:1/r,3:K(s)/(q*a)}
    if mode=='bd': gam={i:gam[i]*(K(s)**(i*(i-1)//2)) for i in gam}
    sg=lambda b,i,mu: (s**((mu*(B(b,i)+1))%2)) if mode=='slice' else 1
    def dot(b,i): b=list(b); b[i-1]+=-pair(b,i)+(l[i-1]+1); return tuple(b)
    cls=lambda b:(b[0]%n,b[1]%n)
    mu=lambda b,i: l[i-1]+1-pair(b,i)
    X=(x1,x2)
    def P(b,i): x=X[i-1]; return (q*x)**(l[i-1]+1-(mu(b,i)%m[i]))*(1-1/q)/(1-(q*x)**m[i]/q)
    def Q(b,i): x=X[i-1]; return -sg(b,i,mu(b,i))*gam[(-L2[i]*mu(b,i))%n]*(q*x)**(l[i-1]+1-m[i])*(1-(q*x)**m[i])/(1-(q*x)**m[i]/q)
    def img(j): return [r,a]+[q**(-Cm[(i,j)])*X[i-1]*X[j-1]**(-Cm[(i,j)]) for i in (1,2)]
    IM={1:img(1),2:img(2)}
    def act(V,i):
        W={}
        for lam,f in V.items():
            g=ev(f,IM[i]); lam2=cls(dot(lam,i))
            W[lam]=W.get(lam,K(0))+P(lam,i)*g
            W[lam2]=W.get(lam2,K(0))+Q(lam2,i)*g
        return W
    bad=0
    for k1 in range(n):
        for k2 in range(n):
            beta=(k1,k2); mono=x1**k1*x2**k2
            for word in ([1,1],[2,2],[1,2]*6):
                V={cls(beta):mono}
                for i in word: V=act(V,i)
                ok=all(V.get(c,K(0))==({cls(beta):mono}.get(c,K(0))) for c in set(V)|{cls(beta)})
                if not ok: bad+=1
    return bad
mode=sys.argv[1]; s=int(sys.argv[2])
for l in [tuple(map(int,z.split(','))) for z in sys.argv[3].split(';')]:
    t0=time.time(); b=check(s,l,mode)
    print("mode=%s s=%+d l=%s: %s (%.0fs)"%(mode,s,l,"ALL RELATIONS HOLD IDENTICALLY" if b==0 else "FAIL %d"%b,time.time()-t0),flush=True)
