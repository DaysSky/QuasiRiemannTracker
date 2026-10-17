# SYMBOLIC (exact, FLINT) version of thm41b.py: q=r^2, gamma(1)=a indeterminates.
#  (1) N(x;l) = h D is a polynomial for the slice-signed action;
#  (2) Theorem 4.1' : the TRUE slice factor f'_k (raw Gauss sums, self-symbol sign S) satisfies f'_k(x) = (qx)^(l-K(k)) f'_k(1/(q^2x)).
import sys, time
from ratf import F, sub, R, A, X1, X2, ctx, ONE
n=4; L2={1:1,2:3}; Cm={(1,1):2,(2,2):2,(1,2):-1,(2,1):-3}; m={1:4,2:4}
pos=[(1,0),(0,1),(1,1),(2,1),(3,1),(3,2)]
words=[[]]
for k in range(1,7):
    for st in (1,2):
        if k==6 and st==2: continue
        words.append([st if j%2==0 else 3-st for j in range(k)])
pair=lambda b,i: b[0]*Cm[(1,i)]+b[1]*Cm[(2,i)]
def B(b,i): return b[i-1]+sum(Cm[(j,i)]*b[j-1] for j in (1,2) if j>i)
r=F(R); a=F(A); X=(F(X1),F(X2)); q=r*r
def ppart(s,l,mode):
    graw={0:F(-1),1:a,2:r.inv(),3:F(s)/(q*a)}
    gam=graw if mode!='bd' else {i:graw[i]*F(s**(i*(i-1)//2)) for i in graw}
    sg=lambda b,i,mu: (s**((mu*(B(b,i)+1))%2)) if mode=='slice' else 1
    def dot(b,i): b=list(b); b[i-1]+=-pair(b,i)+(l[i-1]+1); return tuple(b)
    cls=lambda b:(b[0]%n,b[1]%n)
    mu=lambda b,i: l[i-1]+1-pair(b,i)
    def P(b,i): x=X[i-1]; return (q*x)**(l[i-1]+1-(mu(b,i)%m[i]))*(1-q.inv())/(1-(q*x)**m[i]/q)
    def Q(b,i): x=X[i-1]; return F(-sg(b,i,mu(b,i)))*gam[(-L2[i]*mu(b,i))%n]*(q*x)**(l[i-1]+1-m[i])*(1-(q*x)**m[i])/(1-(q*x)**m[i]/q)
    def act(V,i):
        W={}
        for lam,f in V.items():
            g=sub(f,i); lam2=cls(dot(lam,i))
            W[lam]=W.get(lam,F(0))+P(lam,i)*g
            W[lam2]=W.get(lam2,F(0))+Q(lam2,i)*g
        return W
    def Delta_w(word):
        # Delta(w x) with w x = sigma_{i1}(sigma_{i2}(...x)): evaluate Delta's factors as monomials, substituting innermost first
        tot=F(1)
        for k1,k2 in pos:
            f=F(1)-(q**(4*(k1+k2)))*X[0]**(4*k1)*X[1]**(4*k2)
            for i in word:   # f(sigma_{i1} y) with y = sigma_{i2}...x : apply subs in order i1, i2, ...
                f=sub(f,i)
            tot=tot*f
        return tot
    h=F(0)
    for w in words:
        V={(0,0):F(1)}
        for i in w: V=act(V,i)
        h=h+sum(V.values(),F(0))/Delta_w(w)
    D=F(1)
    for k1,k2 in pos: D=D*(F(1)-(q**(4*(k1+k2)-1))*X[0]**(4*k1)*X[1]**(4*k2))
    N=h*D
    return N, graw
def main(s,l,mode):
    N,graw=ppart(s,l,mode)
    dg=N.d.degrees()
    if dg[2]!=0 or dg[3]!=0:
        return "N NOT polynomial in x (den degrees %s)"%(dg,), None
    # coefficients of N as polys in (r,a): group by (x1,x2) exponent
    coef={}
    Dra=F(N.d)
    for e,c in N.n.to_dict().items():
        key=(int(e[2]),int(e[3])); coef[key]=coef.get(key,F(0))+F(ctx.from_dict({(e[0],e[1],0,0):c}))
    coef={k:v/Dra for k,v in coef.items()}
    # Theorem 4.1' with x as X1 (reuse variable x1 as the univariate slice variable)
    x=X[0]; tot=bad=0
    for i in (1,2):
        j=3-i; t=L2[i]
        for bj in sorted(set(e[j-1] for e in coef)):
            ll=l[i-1]+(-Cm[(j,i)])*bj
            Nk={k:sum((v*x**e[i-1] for e,v in coef.items() if e[j-1]==bj and e[i-1]%4==k),F(0)) for k in range(4)}
            K=lambda k:(ll-2*k)%4
            for k in range(4):
                kp=(ll+1-k)%4; Kp=K(kp); f=Nk[k]
                if Kp!=3:
                    S=s**(((Kp+1)*(kp+(Cm[(j,i)]*bj if j>i else 0)))%2)
                    f=f-F(S)*q**Kp*q*graw[(t*(Kp+1))%4]*x**(Kp+1)*Nk[kp]
                f=f/(1-q**3*x**4)
                if f.iszero(): continue
                tot+=1
                # f(1/(q^2 x)) = sigma_1 applied to a function of x1 only
                if not (f-(q*x)**(ll-K(k))*sub(f,1)).iszero(): bad+=1
    return "N polynomial; Thm 4.1' classes %d, fail %d"%(tot,bad), coef
if __name__=="__main__":
    mode=sys.argv[1]; s=int(sys.argv[2])
    for l in [tuple(map(int,z.split(','))) for z in sys.argv[3].split(';')]:
        t0=time.time(); msg,_=main(s,l,mode)
        print("mode=%s s=%+d l=%s: %s (%.0fs)"%(mode,s,l,msg,time.time()-t0),flush=True)
