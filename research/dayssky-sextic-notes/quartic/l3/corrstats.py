import sys, numpy as np
def load(fn):
    with open(fn,'rb') as f:
        NR,J=np.frombuffer(f.read(8),dtype=np.int32)
        rN=np.frombuffer(f.read(8*NR),dtype=np.float64)
        d=np.frombuffer(f.read(),dtype=np.float64).reshape(J,NR,2)
    return NR,J,rN,d
for fn in sys.argv[1:]:
    NR,J,rN,d=load(fn)
    X=d[:,:,0]/d[:,:,1]           # |R|^2 / rowdiag per window
    ok=np.all(d[:,:,1]>0,axis=0)&(rN>0)   # squarefree m != 1 only
    X=X[:,ok]
    print(fn, "rows",X.shape[1])
    for j in range(J):
        x=X[j]; m=x.mean()
        print("  window %d: S/diag=%.3f  E[x^2]/E[x]^2=%.2f  E[x^3]/E[x]^3=%.2f  max/mean=%.1f"%(j,m,(x**2).mean()/m**2,(x**3).mean()/m**3,x.max()/m))
    for j in range(J):
        for k in range(j+1,J):
            c=np.corrcoef(X[j],X[k])[0,1]; ls=np.corrcoef(np.log(X[j]+1e-12),np.log(X[k]+1e-12))[0,1]
            print("  corr(win%d,win%d) of x: %.3f   of log x: %.3f   (null s.e. ~ %.3f)"%(j,k,c,ls,1/np.sqrt(X.shape[1])))
