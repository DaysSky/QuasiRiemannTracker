# Group relations s1^2 = s2^2 = (s1 s2)^6 = 1 for the slice-signed action (Q_{beta,i} * s^{mu(B_i(beta)+1)}, raw Gauss sums),
# exact rational evaluation along words, as in l3e/cg_local.py.  Also: BD-signed and raw for comparison.
import random, sys
from fractions import Fraction as Fr
n=4; L2={1:1,2:3}; C={(1,1):2,(2,2):2,(1,2):-1,(2,1):-3}; m={1:4,2:4}
pair=lambda b,i: b[0]*C[(1,i)]+b[1]*C[(2,i)]
def B(b,i): return b[i-1]+sum(C[(j,i)]*b[j-1] for j in (1,2) if j>i)
def run(s,l,seed,mode):
    rnd=random.Random(seed)
    r=Fr(rnd.randint(2,9),rnd.randint(1,3)); q=r*r; a=Fr(rnd.randint(1,9),rnd.randint(1,9))
    gam={0:Fr(-1),1:a,2:1/r,3:Fr(s)/(q*a)}
    if mode=='bd': gam={i:gam[i]*Fr(s)**(i*(i-1)//2) for i in gam}
    def sg(b,i,mu):
        if mode!='slice': return 1
        return s**((mu*(B(b,i)+1))%2)
    def dot(b,i): b=list(b); b[i-1]+=-pair(b,i)+(l[i-1]+1); return tuple(b)
    cls=lambda b:(b[0]%n,b[1]%n)
    mu=lambda b,i: l[i-1]+1-pair(b,i)
    def P(b,i,x): return (q*x)**(l[i-1]+1-(mu(b,i)%m[i]))*(1-1/q)/(1-(q*x)**m[i]/q)
    def Q(b,i,x): return -sg(b,i,mu(b,i))*gam[(-L2[i]*mu(b,i))%n]*(q*x)**(l[i-1]+1-m[i])*(1-(q*x)**m[i])/(1-(q*x)**m[i]/q)
    def sx(j,x): return tuple(q**(-C[(i,j)])*x[i-1]*x[j-1]**(-C[(i,j)]) for i in (1,2))
    def word_eval(beta,word,x0):
        ys=[None]*(len(word)+1); ys[len(word)]=x0
        for j in range(len(word),0,-1): ys[j-1]=sx(word[j-1],ys[j])
        V={cls(beta):ys[0][0]**beta[0]*ys[0][1]**beta[1]}
        for j in range(1,len(word)+1):
            i=word[j-1]; y=ys[j]; W={}
            for lam,v in V.items():
                lam2=cls(dot(lam,i))
                W[lam]=W.get(lam,0)+P(lam,i,y[i-1])*v
                W[lam2]=W.get(lam2,0)+Q(lam2,i,y[i-1])*v
            V=W
        return V
    x0=(Fr(rnd.randint(2,30),rnd.randint(31,60)),Fr(rnd.randint(2,30),rnd.randint(31,60)))
    bad=0
    for k1 in range(n):
        for k2 in range(n):
            beta=(k1,k2)
            for word in ([1,1],[2,2],[1,2]*6):
                V=word_eval(beta,word,x0); want={cls(beta):x0[0]**beta[0]*x0[1]**beta[1]}
                if not all(V.get(c,0)==want.get(c,0) for c in set(V)|set(want)): bad+=1
    return bad
for mode in ('raw','bd','slice'):
    for s in (1,-1):
        res=[run(s,l,seed,mode) for l in [(a,b) for a in range(4) for b in range(4)] for seed in (1,2)]
        print("mode=%-5s s=%+d: %d (l,seed) cases, failures per case: %s"%(mode,s,len(res),"all 0" if not any(res) else res))
