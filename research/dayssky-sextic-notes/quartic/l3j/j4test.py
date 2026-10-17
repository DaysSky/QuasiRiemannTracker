# Test the 4th-power-frequency mechanism behind THEOREM_ZF.
# For coefficient families beta(c) = T(c) and beta(c) = mu(c)T(c) (=: G'), compute
#   S = sum_{m in Z[i]} V(m) |sum_c beta(c) conj(chi_c(m)) W(N c/Z)|^2,  V(m) = exp(-pi N(m)/Q),
# over squarefree primary c with only split prime factors, N(c) in (Z/2, Z], and the diagonal.
# The mu-side theory predicts S_{G'}/diag ~ 1 + const*K^{1/4} (K = Z^2/Q) from k = j^4 frequencies,
# while S_T/diag ~ 1.  Rows: all m (no coset).  Optional class restriction N(c) = 1 mod 8.
import sys, math, numpy as np
sys.path.insert(0,'..')
from gauss4 import gauss, gfactor, is_primary, gnorm, chi
def primes_upto(n):
    s=bytearray([1])*(n+1); s[0]=s[1]=0
    for i in range(2,int(n**.5)+1):
        if s[i]: s[i*i::i]=bytearray(len(s[i*i::i]))
    return [i for i in range(n+1) if s[i]]
def run(Z,Q,cls=None,seed=0):
    ps=[p for p in primes_upto(Z) if p%4==1]
    # split primes pi = a+bi primary, both conjugates
    P=[]
    for p in ps:
        for a in range(1,int(p**.5)+1):
            b2=p-a*a; b=int(round(b2**.5))
            if b*b==b2: break
        for (x,y) in [(a,b),(a,-b),(b,a),(b,-a)]:
            for u in [(1,0),(0,1),(-1,0),(0,-1)]:
                z=(x*u[0]-y*u[1],x*u[1]+y*u[0])
                if is_primary(z) and z not in [q[0] for q in P if q[1]==p]:
                    P.append((z,p))
    P=sorted(set(P),key=lambda t:(t[1],t[0]))
    # rows
    R=int(math.sqrt(6*Q))+1
    A,B=np.meshgrid(np.arange(-R,R+1),np.arange(-R,R+1),indexing='ij')
    A=A.ravel(); B=B.ravel(); Nm=A*A+B*B; keep=(Nm>0)&(Nm<=6*Q)
    A=A[keep]; B=B[keep]; Nm=Nm[keep]; V=np.exp(-math.pi*Nm/Q)
    # per prime: exponent of (m/pi)_4 and quartic char table
    E={}; G={}
    for (pi,p) in P:
        a,b=pi
        # r with r = i mod pi:  a+b r = 0 mod p -> r = -a/b
        r=(-a*pow(b,-1,p))%p
        assert (r*r+1)%p==0
        # primitive root
        g=2
        while any(pow(g,(p-1)//q,p)==1 for q in set(f for f in range(2,p) if (p-1)%f==0 and all(f%d for d in range(2,int(f**.5)+1)))): g+=1
        lg=np.zeros(p,dtype=np.int64); y=1
        for s in range(p-1): lg[y]=s; y=y*g%p
        # quartic symbol (x/pi)_4 = x^{(p-1)/4} mod pi; identify with i^k via r
        z=pow(g,(p-1)//4,p)  # = i^{eps} mod pi with eps = 1 if z==r else 3
        eps=1 if z==r else 3
        assert z in (r,p-r)
        x=(A%p + (B%p)*r)%p
        ex=np.where(x==0,-1,(eps*lg[x])%4)
        E[pi]=ex
        # Gauss sum g(conj chi_pi, 1) = sum_{x mod p} conj(chi(x)) e(Im(x/pi)), x rational reps
        xs=np.arange(1,p); ch=(eps*lg[xs])%4
        ph=np.exp(2j*math.pi*(-xs*b% p)/p)  # Im(x/pi)= Im(x*conj(pi))/p = -x b/p
        gs=np.sum(np.exp(-0.5j*math.pi*ch)*ph)
        G[pi]=gs
    # build squarefree c from primes (split only), N(c) in (Z/2,Z]
    cs=[]
    def dfs(i,c,N,fac):
        if N>Z//2: cs.append((c,N,tuple(fac)))
        for j in range(i,len(P)):
            pi,p=P[j]
            if N*p>Z: break
            if any(f[1]==p for f in fac): continue   # avoid pi and conj(pi)? allowed: conj is different prime; keep both allowed
            nc=(c[0]*pi[0]-c[1]*pi[1], c[0]*pi[1]+c[1]*pi[0])
            dfs(j+1,nc,N*p,fac+[(pi,p)])
    def dfs2(i,c,N,fac):
        if N>Z//2: cs.append((c,N,tuple(fac)))
        for j in range(i,len(P)):
            pi,p=P[j]
            if N*p>Z: break
            nc=(c[0]*pi[0]-c[1]*pi[1], c[0]*pi[1]+c[1]*pi[0])
            dfs2(j+1,nc,N*p,fac+[(pi,p)])
    dfs2(0,(1,0),1,[])
    if cls is not None: cs=[t for t in cs if ((t[1]-1)//4)%2==cls]
    W=lambda x: math.exp(-1/((x-.5)*(1-x))+16) if .5<x<1 else 0.0
    RT=np.zeros(len(A),dtype=complex); RG=np.zeros(len(A),dtype=complex); diag=0.0
    for (c,N,fac) in cs:
        w=W(N/Z)
        if w==0: continue
        # T(c) via twisted multiplicativity from prime values
        mu=(-1)**len(fac)
        # build T recursively: T(ab)=T(a)T(b)(a/b)_2 s(a,b)
        Tc=1+0j; Nc=1; cc=(1,0)
        for (pi,p) in fac:
            Tp=complex(*pi)/math.sqrt(p)*G[pi]/math.sqrt(p)
            if Nc>1:
                # (cc/pi)_2 = chi_pi(cc)^2
                q=chi(cc,pi,[pi])**2
                s=-1 if (Nc%8==5 and p%8==5) else 1
                Tc=Tc*Tp*q*s
            else: Tc=Tp
            cc=(cc[0]*pi[0]-cc[1]*pi[1], cc[0]*pi[1]+cc[1]*pi[0]); Nc*=p
        ex=np.zeros(len(A),dtype=np.int64); zero=np.zeros(len(A),dtype=bool)
        for (pi,p) in fac:
            e=E[pi]; zero|=(e<0); ex+=np.where(e<0,0,e)
        val=np.where(zero,0,np.exp(-0.5j*math.pi*(ex%4)))   # conj chi_c(m)
        RT+=Tc*w*val; RG+=mu*Tc*w*val
        diag+=w*w*np.sum(V*(~zero))
    ST=np.sum(V*abs(RT)**2); SG=np.sum(V*abs(RG)**2)
    return len(cs),ST/diag,SG/diag
if __name__=="__main__":
    # sanity: T(c) vs brute force for a few small c
    for Z in [int(a) for a in sys.argv[1:]]:
        n,rt,rg=run(Z,Z)
        print(f"Z=Q={Z}: #c={n}  S_T/diag={rt:.4f}  S_mu*T/diag={rg:.4f}  K^(1/4)={Z**0.25:.2f}",flush=True)
