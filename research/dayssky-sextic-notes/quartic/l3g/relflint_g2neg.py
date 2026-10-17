# SYMBOLIC proof (exact, FLINT) of s1^2 = s2^2 = (s1 s2)^6 = 1 for the slice-signed action, q=r^2, gamma(1)=a indeterminates.
import sys, time
from ratf import F, sub, R, A, X1, X2
n=4; L2={1:1,2:3}; Cm={(1,1):2,(2,2):2,(1,2):-1,(2,1):-3}; m={1:4,2:4}
pair=lambda b,i: b[0]*Cm[(1,i)]+b[1]*Cm[(2,i)]
def B(b,i): return b[i-1]+sum(Cm[(j,i)]*b[j-1] for j in (1,2) if j>i)
r=F(R); a=F(A); X=(F(X1),F(X2)); q=r*r
def check(s,l,mode):
    gam={0:F(-1),1:a,2:-r.inv(),3:F(s)/(q*a)}
    if mode=='bd': gam={i:gam[i]*F(s**(i*(i-1)//2)) for i in gam}
    sg=lambda b,i,mu: (s**((mu*(B(b,i)+1))%2)) if mode=='slice' else 1
    def dot(b,i): b=list(b); b[i-1]+=-pair(b,i)+(l[i-1]+1); return tuple(b)
    cls=lambda b:(b[0]%n,b[1]%n)
    mu=lambda b,i: l[i-1]+1-pair(b,i)
    cacheP={}; cacheQ={}
    def P(b,i):
        k=(cls(b),i)
        if k not in cacheP:
            x=X[i-1]; cacheP[k]=(q*x)**(l[i-1]+1-(mu(b,i)%m[i]))*(1-q.inv())/(1-(q*x)**m[i]/q)
        return cacheP[k]
    def Q(b,i):
        k=(cls(b),i)
        if k not in cacheQ:
            x=X[i-1]; cacheQ[k]=F(-sg(b,i,mu(b,i)))*gam[(-L2[i]*mu(b,i))%n]*(q*x)**(l[i-1]+1-m[i])*(1-(q*x)**m[i])/(1-(q*x)**m[i]/q)
        return cacheQ[k]
    def act(V,i):
        W={}
        for lam,f in V.items():
            g=sub(f,i); lam2=cls(dot(lam,i))
            W[lam]=W.get(lam,F(0))+P(lam,i)*g
            W[lam2]=W.get(lam2,F(0))+Q(lam2,i)*g
        return W
    bad=0
    for k1 in range(n):
        for k2 in range(n):
            beta=(k1,k2); mono=X[0]**k1*X[1]**k2
            for word in ([1,1],[2,2],[1,2]*6):
                V={cls(beta):mono}
                for i in word: V=act(V,i)
                want={cls(beta):mono}
                if not all((V.get(c,F(0))-want.get(c,F(0))).iszero() for c in set(V)|set(want)): bad+=1
    return bad
mode=sys.argv[1]; s=int(sys.argv[2])
for l in [tuple(map(int,z.split(','))) for z in sys.argv[3].split(';')]:
    t0=time.time(); b=check(s,l,mode)
    print("mode=%s s=%+d l=%s: %s (%.1fs)"%(mode,s,l,"all 48 relations hold IDENTICALLY in Q(r,a,x1,x2)" if b==0 else "FAIL %d of 48"%b,time.time()-t0),flush=True)
