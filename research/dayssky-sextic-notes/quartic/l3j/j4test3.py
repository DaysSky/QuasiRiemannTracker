# j4test with inert primes included in c (to match l3/gen_rows.c), primary odd rows, sharp cutoff; compare.
import sys, math, cmath, numpy as np
sys.path.insert(0,'..')
from gauss4 import gauss, gfactor, is_primary, gnorm, chi, gmul
def primes_upto(n):
    s=bytearray([1])*(n+1); s[0]=s[1]=0
    for i in range(2,int(n**.5)+1):
        if s[i]: s[i*i::i]=bytearray(len(s[i*i::i]))
    return [i for i in range(n+1) if s[i]]
def cpow_mod(a,b,e,q):
    ra=np.ones_like(a); rb=np.zeros_like(b); xa=a%q; xb=b%q
    while e:
        if e&1: ra,rb=(ra*xa-rb*xb)%q,(ra*xb+rb*xa)%q
        xa,xb=(xa*xa-xb*xb)%q,(2*xa*xb)%q; e>>=1
    return ra,rb
def run(Z,Q,rows='primary',weight='sharp',inert=True):
    P=[]
    for p in primes_upto(Z):
        if p%4==1:
            for a in range(1,int(p**.5)+1):
                b2=p-a*a; b=int(round(b2**.5))
                if b*b==b2: break
            got=set()
            for (x,y) in [(a,b),(a,-b),(b,a),(b,-a)]:
                for u in [(1,0),(0,1),(-1,0),(0,-1)]:
                    z=(x*u[0]-y*u[1],x*u[1]+y*u[0])
                    if is_primary(z): got.add(z)
            for z in got: P.append((z,p,False))
        elif p%4==3 and inert and p*p<=Z:
            P.append(((-p,0),p*p,True))
    P.sort(key=lambda t:(t[1],t[0]))
    R=int(math.sqrt(6*Q))+1
    A,B=np.meshgrid(np.arange(-R,R+1),np.arange(-R,R+1),indexing='ij')
    A=A.ravel(); B=B.ravel(); Nm=A*A+B*B; keep=(Nm>0)&(Nm<=6*Q)
    A=A[keep]; B=B[keep]; Nm=Nm[keep]
    if rows=='primary':
        sel=((A%4==1)&(B%4==0))|((A%4==3)&(B%4==2)); A=A[sel];B=B[sel];Nm=Nm[sel]
    V=np.exp(-math.pi*Nm/Q) if weight=='gauss' else (Nm<=Q).astype(float)
    E={}; G={}
    for (pi,Np,ine) in P:
        if not ine:
            p=Np; a,b=pi; r=(-a*pow(b,-1,p))%p
            fs=set(); n=p-1; f=2
            while f*f<=n:
                while n%f==0: fs.add(f); n//=f
                f+=1
            if n>1: fs.add(n)
            g=2
            while any(pow(g,(p-1)//qq,p)==1 for qq in fs): g+=1
            lg=np.zeros(p,dtype=np.int64); y=1
            for s in range(p-1): lg[y]=s; y=y*g%p
            z=pow(g,(p-1)//4,p); eps=1 if z==r else 3
            x=(A%p+(B%p)*r)%p
            E[pi]=np.where(x==0,-1,(eps*lg[x])%4)
            xs=np.arange(1,p); ch=(eps*lg[xs])%4
            G[pi]=np.sum(np.exp(-0.5j*math.pi*ch)*np.exp(2j*math.pi*((-xs*b)%p)/p))
        else:
            q=-pi[0]; e=(q*q-1)//4
            ra,rb=cpow_mod(A,B,e,q)
            zero=(A%q==0)&(B%q==0)
            ex=np.full(len(A),-1)
            for k,(ua,ub) in enumerate([(1,0),(0,1),(q-1,0),(0,q-1)]):
                ex=np.where((ra==ua)&(rb==ub)&~zero,k,ex)
            E[pi]=ex
            aa,bb=np.meshgrid(np.arange(q),np.arange(q),indexing='ij'); aa=aa.ravel(); bb=bb.ravel()
            ra2,rb2=cpow_mod(aa,bb,e,q); nz=~((aa==0)&(bb==0))
            ch=np.full(len(aa),-1)
            for k,(ua,ub) in enumerate([(1,0),(0,1),(q-1,0),(0,q-1)]):
                ch=np.where((ra2==ua)&(rb2==ub)&nz,k,ch)
            ok=ch>=0
            # x=aa+bb i ; Im(x/pi) with pi=-q : Im(x/(-q)) = -bb/q
            G[pi]=np.sum(np.exp(-0.5j*math.pi*ch[ok])*np.exp(2j*math.pi*((-bb[ok])%q)/q))
    cs=[]
    def dfs(i,c,N,fac):
        if N>Z//2: cs.append((c,N,tuple(fac)))
        for j in range(i,len(P)):
            pi,Np,ine=P[j]
            if N*Np>Z: break
            dfs(j+1,gmul(c,pi),N*Np,fac+[(pi,Np)])
    dfs(0,(1,0),1,[])
    W=lambda x: math.exp(-1/((x-.5)*(1-x))+16) if .5<x<1 else 0.0
    RT=np.zeros(len(A),dtype=complex); RG=np.zeros(len(A),dtype=complex); diag=0.0
    for (c,N,fac) in cs:
        w=W(N/Z)
        if w==0: continue
        mu=(-1)**len(fac); Tc=None; cc=(1,0); Nc=1
        for (pi,Np) in fac:
            Tp=complex(*pi)/math.sqrt(Np)*G[pi]/math.sqrt(Np)
            if Tc is None: Tc=Tp
            else:
                qq=chi(cc,pi,[pi])**2
                s=-1 if (Nc%8==5 and Np%8==5) else 1
                Tc=Tc*Tp*qq*s
            cc=gmul(cc,pi); Nc*=Np
        ex=np.zeros(len(A),dtype=np.int64); zero=np.zeros(len(A),dtype=bool)
        for (pi,Np) in fac:
            e=E[pi]; zero|=(e<0); ex+=np.where(e<0,0,e)
        val=np.where(zero,0,np.exp(-0.5j*math.pi*(ex%4)))
        RT+=Tc*w*val; RG+=mu*Tc*w*val
        diag+=w*w*np.sum(V*(~zero))
    return len(cs),np.sum(V*abs(RT)**2)/diag,np.sum(V*abs(RG)**2)/diag
if __name__=="__main__":
    # sanity: inert prime T vs brute force
    from tcheck import T_brute
    Z=int(sys.argv[1])
    for inert in [False,True]:
        n,rt,rg=run(Z,Z,inert=inert)
        print(f"Z=Q={Z} inert={inert}: #c={n} S_T/diag={rt:.4f} S_muT/diag={rg:.4f}",flush=True)
