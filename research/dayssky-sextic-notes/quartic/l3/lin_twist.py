import sys, numpy as np
sys.path.insert(0,'..'); sys.path.insert(0,'.')
import gauss4 as g4
def load(fn):
    with open(fn,'rb') as f:
        NR,J=np.frombuffer(f.read(8),dtype=np.int32)
        rows=np.frombuffer(f.read(8*6*NR),dtype=np.float64).reshape(NR,6)
        win=np.frombuffer(f.read(),dtype=np.float64).reshape(J,NR,3)
    return NR,J,rows,win
NR,J,rows,win=load(sys.argv[1])
N,sq,a,b,gr,gi=rows.T
ok=np.where((sq==1)&(N>1))[0]
g=(gr+1j*gi)/np.sqrt(N)
hs=[(1,0),(0,1),(1,1),(3,0),(2,1),(1,2),(-1,2),(3,2),(5,0),(1,4)]
fac={}
W={}
for h in hs:
    w=[]
    for i in ok:
        n=(int(a[i]),int(b[i]))
        c=g4.chi(h,n)
        w.append(np.conj(g[i])*np.conj(c))
    W[h]=np.array(w)
rng=np.random.default_rng(2)
for j in range(J):
    R=(win[j,:,0]+1j*win[j,:,1])[ok]; S=np.sum(np.abs(R)**2)
    out=["h=%s:%.2f"%(str(h),abs(np.sum(W[h]*R))/np.sqrt(S)) for h in hs]
    print("win",j," ".join(out))
    # cumulative growth for h=1 (rows sorted by norm)
    order=np.argsort(N[ok]); cs=np.cumsum((W[(1,0)]*R)[order]); cS=np.cumsum(np.abs(R[order])**2)
    idx=[int(len(order)*f)-1 for f in (0.1,0.2,0.4,0.6,0.8,1.0)]
    print("   cumulative |X|/sqrt(S) at row fractions .1,.2,.4,.6,.8,1:"," ".join("%.2f"%(abs(cs[k])/np.sqrt(cS[k])) for k in idx))
