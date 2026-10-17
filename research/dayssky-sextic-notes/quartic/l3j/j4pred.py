# Predicted contribution of the "trivial-character" dual frequencies k = u*j^4 (u a unit) to S_{G'} - diag,
# for G' = mu*T, rows = all m in Z[i] with V(m)=exp(-pi N(m)/Q), c squarefree primary split-or-inert, N(c) in (Z/2,Z].
# Poisson: sum_m V(m) phi(m) = (Q/N(q)) sum_k ghat(k) exp(-pi Q N(k)/N(q)).
# Coprime pair term at frequency k: s(c1,c2) lambda(c2) chi_c1(k) conj chi_c2(k) sqrt(N c1 c2) W1 W2 (Q/N q) exp(...).
# Compare with the observed (S_G' - diag)/diag from j4test-style direct computation.
import sys, math, numpy as np
sys.path.insert(0,'..')
from gauss4 import chi, gmul, gnorm
import j4test3 as J
def collect_c(Z,inert=True):
    # reuse prime list builder from j4test3 via a tiny run with Q small
    P=[]
    for p in J.primes_upto(Z):
        if p%4==1:
            for a in range(1,int(p**.5)+1):
                b2=p-a*a; b=int(round(b2**.5))
                if b*b==b2: break
            got=set()
            for (x,y) in [(a,b),(a,-b),(b,a),(b,-a)]:
                for u in [(1,0),(0,1),(-1,0),(0,-1)]:
                    z=(x*u[0]-y*u[1],x*u[1]+y*u[0])
                    if J.is_primary(z): got.add(z)
            for z in got: P.append((z,p))
        elif p%4==3 and inert and p*p<=Z: P.append(((-p,0),p*p))
    P.sort(key=lambda t:(t[1],t[0]))
    cs=[]
    def dfs(i,c,N,fac):
        if N>Z//2: cs.append((c,N,tuple(f[0] for f in fac)))
        for j in range(i,len(P)):
            pi,Np=P[j]
            if N*Np>Z: break
            dfs(j+1,gmul(c,pi),N*Np,fac+[P[j]])
    dfs(0,(1,0),1,[])
    return cs
def predict(Z,Q):
    cs=collect_c(Z)
    W=lambda x: math.exp(-1/((x-.5)*(1-x))+16) if .5<x<1 else 0.0
    cs=[(c,N,f,W(N/Z)) for (c,N,f) in cs if W(N/Z)>0]
    n=len(cs); Ns=np.array([t[1] for t in cs],float); Ws=np.array([t[3] for t in cs])
    e=((Ns.astype(np.int64)-1)//4)%2
    facs=[set(t[2]) for t in cs]
    # coprimality matrix
    cop=np.ones((n,n),bool)
    for a in range(n):
        for b in range(n):
            if facs[a]&facs[b]: cop[a,b]=False
    np.fill_diagonal(cop,False)
    sgn=(-1.0)**(np.outer(e,e)+e[None,:])
    # j range: exp(-pi Q N(k)/N(q)) with N(q)<=Z^2 -> need N(k) up to ~ 3 Z^2/Q
    Kmax=4*Z*Z/Q
    tot=0.0
    js=[(a,b) for a in range(0,40) for b in range(-40,40) if 0<a*a+b*b and (a*a+b*b)**4<=Kmax and (a>0 or b>0)]
    # j up to associates: j^4 same for j, ij, -j, -ij -> take one per class: a>0,b>=0 or so
    js=[(a,b) for (a,b) in js if a>0 and b>=0]
    units=[(1,0),(0,1),(-1,0),(0,-1)]
    for (a,b) in js:
        j=(a,b); j2=gmul(j,j); j4=gmul(j2,j2); Nk=gnorm(j4)
        jfac=set()
        # primes dividing j: check chi(j,c)==0 via factor list
        for u in units:
            k=gmul(u,j4)
            ch=np.array([chi(k,t[0],list(t[2])) for t in cs])
            M=np.outer(ch,np.conj(ch))
            wgt=np.exp(-math.pi*Q*Nk/np.outer(Ns,Ns))
            term=sgn*M*np.outer(Ws,Ws)*Q/np.sqrt(np.outer(Ns,Ns))*wgt*cop
            tot+=term.sum().real
    # diag (approx): sum_c W^2 * sum_m V(m)*1_{(m,c)=1} ~ sum_c W^2 * Q * prod(1-1/N(pi))
    diag=0.0
    for (c,N,f,w) in cs:
        dens=1.0
        for pi in f: dens*=(1-1/gnorm(pi))
        diag+=w*w*Q*dens
    return n,tot/diag
for Z in [int(a) for a in sys.argv[1:]]:
    n,r=predict(Z,Z); print(f"Z=Q={Z}: #c={n} predicted trivial-frequency excess / diag = {r:.4f}",flush=True)
