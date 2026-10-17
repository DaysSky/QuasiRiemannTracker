# normalized quartic Gauss sums gt(c)=g(chi_c)/sqrt(Nc), chi_c=(./c)_4, psi(z)=e(Im z), for squarefree primary c, via twisted multiplicativity
import numpy as np, math, sys
from pocore import *
def is_primary(A,B): return (A%4,B%4) in ((1,0),(3,2))
def gens_for(kinds,mods,rs):
    G=[]
    for j in range(len(kinds)):
        p=int(mods[j])
        if kinds[j]==1:
            G.append((-p,0) if is_primary(-p,0) else (p,0)); continue
        r=int(rs[j]); found=None
        a=1
        while a*a<p:
            b2=p-a*a; b=math.isqrt(b2)
            if b*b==b2: break
            a+=1
        for (A,B) in ((a,b),(a,-b),(b,a),(b,-a),(-a,b),(-a,-b),(-b,a),(-b,-a)):
            if is_primary(A,B) and (A+B*r)%p==0: found=(A,B)
        assert found is not None
        G.append(found)
    return G
def sym(kind,mod,r,off,tab,A,B):
    """(A+Bi / pi)_4 exponent, -1 if divisible"""
    if kind==0:
        return int(tab[off+((A+B*r)%mod)])
    return int(tab[off+(A%mod)*mod+(B%mod)])
def prime_gauss(kinds,mods,rs,offs,tab,G):
    out=np.zeros(len(kinds),complex)
    for j in range(len(kinds)):
        p=int(mods[j]); A,B=G[j]
        if kinds[j]==0:
            x=np.arange(1,p); e=tab[offs[j]+x].astype(int)
            out[j]=np.sum((1j)**e*np.exp(-2j*np.pi*B*x/p))/math.sqrt(p)
        else:
            q=p; s=0
            u,v=np.meshgrid(np.arange(q),np.arange(q),indexing='ij'); u=u.ravel(); v=v.ravel()
            e=tab[offs[j]+u*q+v].astype(int); m=e>=0
            # x/pi with pi=(A,0), A=+-q : Im(x/A)=v/A
            s=np.sum((1j)**e[m]*np.exp(2j*np.pi*v[m]/A))
            out[j]=s/q
    return out
def column_gauss(ptr,idx,kinds,mods,rs,offs,tab,G,gp):
    C=len(ptr)-1; gt=np.zeros(C,complex); gen=[]
    cache={}
    def sy(i,j):  # (pi_j / pi_i)_4 as complex
        if (i,j) not in cache:
            A,B=G[j]; e=sym(kinds[i],int(mods[i]),int(rs[i]),int(offs[i]),tab,A,B); cache[(i,j)]=(1j)**e
        return cache[(i,j)]
    for c in range(C):
        ps=idx[ptr[c]:ptr[c+1]]; val=1+0j; n=(1,0)
        for a in range(len(ps)):
            val*=gp[ps[a]]
            for b in range(a):
                val*=sy(ps[a],ps[b])*sy(ps[b],ps[a])
            A,B=G[ps[a]]; n=(n[0]*A-n[1]*B,n[0]*B+n[1]*A)
        gt[c]=val; gen.append(n)
    return gt,np.array(gen)
if __name__=="__main__":
    # sanity: compare with gauss4 for small D
    sys.path.insert(0,'..')
    import gauss4
    D=600
    kinds,mods,rs,offs,norms,tab=build_primes(D)
    G=gens_for(kinds,mods,rs); gp=prime_gauss(kinds,mods,rs,offs,tab,G)
    ptr,idx,cn,om=build_columns(norms,D,lo=1)
    gt,gen=column_gauss(ptr,idx,kinds,mods,rs,offs,tab,G,gp)
    bad=0
    for c in range(0,len(cn),7):
        n=tuple(int(t) for t in gen[c]); assert gauss4.is_primary(n)
        ref=gauss4.gauss(n,1)/math.sqrt(gauss4.gnorm(n))
        if abs(ref-gt[c])>1e-8: bad+=1
    print("checked",len(range(0,len(cn),7)),"bad",bad, "mean|gt|",np.abs(gt).mean())
