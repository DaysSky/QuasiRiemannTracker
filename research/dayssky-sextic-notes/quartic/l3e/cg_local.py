# Chinta-Gunnells local Weyl action (Def. 3.1, (3.8)-(3.14)) for G2, n=4, checked by exact evaluation.
# gamma(i): gamma(0)=-1, gamma(i)gamma(-i) = s^i/q  (CG: s=+1; Q(i) with N(pi)=5 mod 8: s=-1).
import sys, random
from fractions import Fraction as Fr
n=4
L2={1:1,2:3}                     # |alpha|^2: node 1 short (c), node 2 long (d)
C={(1,1):2,(2,2):2,(1,2):-1,(2,1):-3}   # c(i,j)=2<a_i,a_j>/<a_j,a_j>
m={i:n//__import__('math').gcd(n,L2[i]) for i in (1,2)}
def pair(beta,i):                # <beta, alpha_i^vee>, beta=(k1,k2)
    return beta[0]*C[(1,i)]+beta[1]*C[(2,i)]
def dot(beta,l,i):               # sigma_i . beta = sigma_i beta + (l_i+1) alpha_i
    b=list(beta); b[i-1]+= -pair(beta,i)+(l[i-1]+1); return tuple(b)
def cls(beta): return (beta[0]%n, beta[1]%n)   # Lambda' = 4 Lambda here (all m=4)
def mu(beta,l,i): return l[i-1]+1-pair(beta,i)
def run(s,l,seed=1,fix=None):
    rnd=random.Random(seed)
    r=Fr(rnd.randint(2,9),rnd.randint(1,3)); q=r*r
    a=Fr(rnd.randint(1,9),rnd.randint(1,9))
    gam={0:Fr(-1),1:a,2:1/r,3:Fr(s)/(q*a)}
    if fix=='binom': gam={i:gam[i]*Fr(s)**(i*(i-1)//2) for i in gam}   # (pi,pi)^{binom(i,2)}
    if fix=='naive': gam={i:gam[i]*Fr(s)**i for i in gam}
    def P(beta,i,x):
        mm=m[i]; e=l[i-1]+1-(mu(beta,l,i)%mm)
        return (q*x)**e*(1-1/q)/(1-(q*x)**mm/q)
    def Q(beta,i,x):
        mm=m[i]; e=l[i-1]+1-mm
        return -gam[(-L2[i]*mu(beta,l,i))%n]*(q*x)**e*(1-(q*x)**mm)/(1-(q*x)**mm/q)
    def sx(j,x):
        return tuple(q**(-C[(i,j)])*x[i-1]*x[j-1]**(-C[(i,j)]) for i in (1,2))
    def word_eval(beta,word,x0):
        # points: y_k = x0, y_{j-1} = sigma_{i_j}(y_j); V_j[class] = value of that component of f|s_{i1}..s_{ij} at y_j
        ys=[None]*(len(word)+1); ys[len(word)]=x0
        for j in range(len(word),0,-1): ys[j-1]=sx(word[j-1],ys[j])
        V={cls(beta):ys[0][0]**beta[0]*ys[0][1]**beta[1]}
        for j in range(1,len(word)+1):
            i=word[j-1]; y=ys[j]; W={}
            for lam,v in V.items():
                lam2=cls(dot(lam,l,i))
                W[lam]=W.get(lam,0)+P(lam,i,y[i-1])*v
                W[lam2]=W.get(lam2,0)+Q(lam2,i,y[i-1])*v
            V=W
        return V
    x0=(Fr(rnd.randint(2,30),rnd.randint(31,60)),Fr(rnd.randint(2,30),rnd.randint(31,60)))
    bad=[]
    for k1 in range(n):
        for k2 in range(n):
            beta=(k1,k2)
            for word,name in (([1,1],"s1^2"),([2,2],"s2^2"),([1,2]*6,"(s1s2)^6")):
                V=word_eval(beta,word,x0)
                want={cls(beta):x0[0]**beta[0]*x0[1]**beta[1]}
                ok=all(V.get(c,0)==want.get(c,0) for c in set(V)|set(want))
                if not ok: bad.append((beta,name))
    return bad
for s in (1,-1):
    for fix in (None,'naive','binom'):
        for l in ((0,0),(1,0),(0,1),(2,1)):
            for seed in (1,2):
                bad=run(s,l,seed,fix)
                print("s=%+d fix=%-5s l=%s seed=%d: %s"%(s,fix,l,seed,"all relations hold" if not bad else "FAIL %d"%len(bad)))
