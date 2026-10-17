# Core for the (PO) experiment:
#   S_a(H,D) = sum_{0<N(k)<=H} | sum_{c sqfree odd ideal, D/2<N(c)<=D} a_c chi_c(k) |^2 ,  chi_c(k) = (k/c)_4
# chi_c(k) depends only on the ideal (c); for prime pi | p (p=1 mod 4) with i -> r mod pi:
#   (k/pi)_4 = i^e  where (x+y r)^{(p-1)/4} == r^e mod p  (k=x+iy);  for inert q: (k/q)_4 = i^e with k^{(q^2-1)/4} == i^e in F_q[i].
import numpy as np, math, numba as nb

def primes_upto(n):
    s=np.ones(n+1,bool); s[:2]=False
    for p in range(2,int(n**.5)+1):
        if s[p]: s[p*p::p]=False
    return np.nonzero(s)[0]

def prim_root(p):
    m=p-1; f=[]; d=2
    while d*d<=m:
        if m%d==0:
            f.append(d)
            while m%d==0: m//=d
        d+=1
    if m>1: f.append(m)
    for g in range(2,p):
        if all(pow(g,(p-1)//q,p)!=1 for q in f): return g

@nb.njit(cache=True)
def _split_table(p,g,r,out):
    # out[t] = e with t^{(p-1)/4} = r^e ; out[0]=-1
    z=1
    for _ in range((p-1)//4): z=z*g%p
    # z=g^{(p-1)/4}; find e0 with z=r^e0
    rr=1; e0=-1
    for e in range(4):
        if rr==z: e0=e
        rr=rr*r%p
    x=1
    for j in range(p-1):
        out[x]=(j*e0)%4
        x=x*g%p
    out[0]=-1

@nb.njit(cache=True)
def _inert_table(q,out):
    E=(q*q-1)//4
    for u in range(q):
        for v in range(q):
            if u==0 and v==0:
                out[0]=-1; continue
            # power (u+iv)^E mod q
            ru,rv=1,0; bu,bv=u,v; e=E
            while e>0:
                if e&1:
                    ru,rv=(ru*bu-rv*bv)%q,(ru*bv+rv*bu)%q
                bu,bv=(bu*bu-bv*bv)%q,(2*bu*bv)%q
                e>>=1
            if ru==1 and rv==0: k=0
            elif ru==0 and rv==1: k=1
            elif ru==q-1 and rv==0: k=2
            elif ru==0 and rv==q-1: k=3
            else: k=-9
            out[u*q+v]=k

def build_primes(D):
    """prime ideals of norm<=D (odd). returns arrays: kind(0 split,1 inert), mod, r, off, norm, table"""
    kinds=[];mods=[];rs=[];offs=[];norms=[];tabs=[];off=0
    for p in primes_upto(D):
        p=int(p)
        if p==2: continue
        if p%4==1:
            g=prim_root(p); r=pow(g,(p-1)//4,p)  # r^2=-1
            for rr in (r,p-r):
                t=np.empty(p,np.int8); _split_table(p,g,rr,t)
                kinds.append(0);mods.append(p);rs.append(rr);offs.append(off);norms.append(p);tabs.append(t);off+=p
        else:
            if p*p>D: continue
            t=np.empty(p*p,np.int8); _inert_table(p,t)
            assert t.min()>=-1
            kinds.append(1);mods.append(p);rs.append(0);offs.append(off);norms.append(p*p);tabs.append(t);off+=p*p
    return (np.array(kinds,np.int64),np.array(mods,np.int64),np.array(rs,np.int64),
            np.array(offs,np.int64),np.array(norms,np.int64),np.concatenate(tabs))

def build_columns(norms,D,lo=None):
    """squarefree products of distinct prime ideals with lo<N<=D (default lo=D/2). CSR."""
    if lo is None: lo=D/2
    order=np.argsort(norms,kind='stable'); nl=norms[order]
    cols=[]; 
    def rec(start,prod,lst):
        if prod>lo and lst: cols.append((prod,list(lst)))
        for t in range(start,len(order)):
            if prod*nl[t]>D: break
            lst.append(int(order[t])); rec(t+1,prod*int(nl[t]),lst); lst.pop()
    rec(0,1,[])
    ptr=[0];idx=[];cn=[];om=[]
    for prod,l in cols:
        idx+=l; ptr.append(len(idx)); cn.append(prod); om.append(len(l))
    return np.array(ptr,np.int64),np.array(idx,np.int64),np.array(cn,np.int64),np.array(om,np.int64)

def Wfun(t):
    t=np.asarray(t,float); u=4*t-3; out=np.zeros_like(t)
    m=np.abs(u)<1; out[m]=np.exp(1-1/(1-u[m]**2)); return out

@nb.njit(parallel=True,cache=True)
def kernel(xs,ys,kinds,mods,rs,offs,tab,ptr,idx,coef_re,coef_im,w2):
    nr=xs.shape[0]; P=kinds.shape[0]; C=ptr.shape[0]-1; V=coef_re.shape[1]
    Fre=np.zeros((nr,V)); Fim=np.zeros((nr,V)); diag=np.zeros(nr)
    for i in nb.prange(nr):
        x=xs[i]; y=ys[i]
        e=np.empty(P,np.int8)
        for j in range(P):
            m=mods[j]
            if kinds[j]==0:
                t=(x+y*rs[j])%m
                if t<0: t+=m
                e[j]=tab[offs[j]+t]
            else:
                u=x%m; v=y%m
                if u<0: u+=m
                if v<0: v+=m
                e[j]=tab[offs[j]+u*m+v]
        acc_re=np.zeros(V); acc_im=np.zeros(V); dg=0.0
        for c in range(C):
            s=0; z=False
            for t in range(ptr[c],ptr[c+1]):
                ej=e[idx[t]]
                if ej<0:
                    z=True; break
                s+=ej
            if z: continue
            dg+=w2[c]
            s&=3
            for v in range(V):
                a=coef_re[c,v]; b=coef_im[c,v]
                if s==0:
                    acc_re[v]+=a; acc_im[v]+=b
                elif s==1:
                    acc_re[v]-=b; acc_im[v]+=a
                elif s==2:
                    acc_re[v]-=a; acc_im[v]-=b
                else:
                    acc_re[v]+=b; acc_im[v]-=a
        for v in range(V):
            Fre[i,v]=acc_re[v]; Fim[i,v]=acc_im[v]
        diag[i]=dg
    return Fre,Fim,diag
