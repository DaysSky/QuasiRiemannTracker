# Numerical test of the quartic dual mean square
#   S(H) = sum_{0<N(h)<=H} | sum_{n sf primary} a(n) chi_n(h) W(N n / X) |^2
# for a(n) = conj(alpha(n))^k * gamma_1(n)  (k=1 is the coefficient from the Mobius identity),
# compared with random unimodular a(n).
import numpy as np, math, sys, time
from gauss4 import primary_assoc, is_primary, gmul, gnorm

def primes_upto(n):
    s=np.ones(n+1,bool); s[:2]=False
    for p in range(2,int(n**.5)+1):
        if s[p]: s[p*p::p]=False
    return np.nonzero(s)[0]

def prim_root(p):
    f=[int(q) for q in primes_upto(int(p**.5)+1) if (p-1)%q==0]
    m=p-1
    for q in list(f):
        while m%q==0: m//=q
    if m>1: f.append(m)
    for g in range(2,p):
        if all(pow(g,(p-1)//q,p)!=1 for q in f): return g

def build(X,H,seed=0):
    t0=time.time()
    # rows h: all nonzero Gaussian integers with N(h)<=H
    R=int(math.isqrt(H))+1
    hs=[(c,d) for c in range(-R,R+1) for d in range(-R,R+1) if 0<c*c+d*d<=H]
    hs.sort(key=lambda z:z[0]**2+z[1]**2)
    hc=np.array([z[0] for z in hs],np.int64); hd=np.array([z[1] for z in hs],np.int64)
    hN=hc*hc+hd*hd
    # primary primes with norm <= X
    P=[]  # (pi, norm, gamma1, chi_exp over rows (int8, -1 for zero))
    for p in primes_upto(X):
        p=int(p)
        if p==2: continue
        if p%4==1:
            # find a,b
            a=next(a for a in range(1,int(p**.5)+1) if int(round((p-a*a)**.5))**2==p-a*a)
            b=int(round((p-a*a)**.5))
            g=prim_root(p)
            ind=np.zeros(p,np.int64); x=1
            for t in range(p-1): ind[x]=t; x=x*g%p
            for pi in (primary_assoc((a,b)),primary_assoc((a,-b))):
                A,B=pi
                r=(-A*pow(B,-1,p))%p           # image of i in F_p
                z=pow(g,(p-1)//4,p)
                e=1 if z==r else (3 if z==(p-r)%p else None)
                assert e is not None
                # chi on integers x: i^(e*ind(x))
                xs=np.arange(1,p)
                k=(e*ind[xs])%4
                gs=np.sum((1j**k)*np.exp(-2j*np.pi*B*xs/p))
                gam=gs/math.sqrt(p)
                # rows
                y=(hc+hd*r)%p
                ex=np.where(y==0,-1,(e*ind[y])%4).astype(np.int8)
                P.append((pi,p,gam,ex))
        else:
            q=p
            if q*q>X: continue
            pi=primary_assoc((q,0))
            # F_{q^2} arithmetic via Gaussian integers mod q
            E=(q*q-1)//4
            def powmod(z,e):
                r=(1,0)
                while e:
                    if e&1: r=((r[0]*z[0]-r[1]*z[1])%q,(r[0]*z[1]+r[1]*z[0])%q)
                    z=((z[0]*z[0]-z[1]*z[1])%q,(2*z[0]*z[1])%q); e>>=1
                return r
            U={(1,0):0,(0,1):1,(q-1,0):2,(0,q-1):3}
            gs=0
            for u in range(q):
                for v in range(q):
                    if u==0 and v==0: continue
                    k=U[powmod((u,v),E)]
                    # psi(x/pi), pi=-q (or q): x/pi = x*conj(pi)/q^2 ; Im
                    num=gmul((u,v),(pi[0],-pi[1]))
                    gs+=(1j**k)*np.exp(2j*np.pi*num[1]/(q*q))
            gam=gs/q
            ex=np.empty(len(hs),np.int8)
            for i,(c,d) in enumerate(hs):
                if c%q==0 and d%q==0: ex[i]=-1
                else: ex[i]=U[powmod((c%q,d%q),E)]
            P.append((pi,q*q,gam,ex))
    t1=time.time()
    return hs,hN,P,t1-t0

def chi_pair(pi1,pi2,P_by):
    pass

def run(X,H,ks=(1,0,2,-1),nrand=2,seed=1):
    hs,hN,P,tb=build(X,H)
    idx={p[0]:j for j,p in enumerate(P)}
    # quartic symbol between primes: chi_pi(pi2) computed via F_p map
    def chi_pp(i,j):
        pi,p,_,_=P[i]; pi2=P[j][0]
        if p%4==1 or True:
            # generic: use gauss4.chi_prime
            from gauss4 import chi_prime
            return chi_prime(pi2,pi)
    # enumerate squarefree primary n with X/2 < N(n) <= X as products of distinct primes
    norms=[p[1] for p in P]
    order=sorted(range(len(P)),key=lambda j:norms[j])
    ns=[]
    def rec(start,prod_norm,lst):
        if prod_norm>X/2 and lst: ns.append(list(lst))
        for t in range(start,len(order)):
            j=order[t]
            if prod_norm*norms[j]>X: break
            lst.append(j); rec(t+1,prod_norm*norms[j],lst); lst.pop()
    rec(0,1,[])
    from gauss4 import chi_prime
    cache={}
    def cpp(i,j):
        if (i,j) not in cache: cache[(i,j)]=chi_prime(P[j][0],P[i][0])
        return cache[(i,j)]
    nH=len(hs)
    results={}
    coefs={k:[] for k in ks}
    rng=np.random.default_rng(seed)
    acc={k:np.zeros(nH,complex) for k in ks}
    accr=[np.zeros(nH,complex) for _ in range(nrand)]
    wsum2=0.0
    for lst in ns:
        # n = product, gamma via twisted multiplicativity
        n=(1,0); gam=1.0+0j; Nn=1
        for t,j in enumerate(lst):
            gam*=P[j][2]
            for s_ in lst[:t]:
                gam*=cpp(s_,j)*cpp(j,s_)
            n=gmul(n,P[j][0]); Nn*=P[j][1]
        alpha=complex(*n)/math.sqrt(Nn)
        tt=Nn/X
        W=math.sin(math.pi*(2*tt-1))**2 if 0.5<tt<=1 else 0.0  # smooth bump on [1/2,1]
        if W==0: continue
        wsum2+=W*W
        # chi_n(h) exponent vector
        ex=np.zeros(nH,np.int16); zero=np.zeros(nH,bool)
        for j in lst:
            e=P[j][3]; zero|=(e<0); ex+=np.where(e<0,0,e)
        v=np.where(zero,0,1j**(ex%4))
        for k in ks:
            acc[k]+=W*(np.conj(alpha)**k)*gam*v
        for r_ in range(nrand):
            ph=np.exp(2j*np.pi*rng.random())
            accr[r_]+=W*ph*v
    return hs,hN,acc,accr,wsum2,len(ns),tb

if __name__=="__main__":
    X=int(sys.argv[1]); H=int(sys.argv[2])
    t=time.time()
    hs,hN,acc,accr,w2,nn,tb=run(X,H)
    print(f"X={X} H={H} rows={len(hs)} columns={nn} sum W^2={w2:.1f} build={tb:.1f}s total={time.time()-t:.1f}s")
    for Hc in [H//16,H//8,H//4,H//2,H]:
        m=hN<=Hc; nh=int(m.sum())
        line=f"  H'={Hc:6d} #h={nh:6d} | "
        for k,a in acc.items():
            S=np.sum(np.abs(a[m])**2)
            line+=f"k={k}: {S/(nh*w2):7.3f}  "
        Sr=np.mean([np.sum(np.abs(a[m])**2) for a in accr])
        line+=f"| rand: {Sr/(nh*w2):7.3f}"
        print(line)
