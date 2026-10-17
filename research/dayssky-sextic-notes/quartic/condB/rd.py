import numpy as np
def load(fn):
    with open(fn,'rb') as f:
        K,_=np.frombuffer(f.read(8),np.int32); K=int(K)
        n=int(np.frombuffer(f.read(8),np.int64)[0])
        N=np.frombuffer(f.read(8*n),np.int64); a=np.frombuffer(f.read(8*n),np.int64); b=np.frombuffer(f.read(8*n),np.int64)
        Gr=np.frombuffer(f.read(8*n),np.float64); Gi=np.frombuffer(f.read(8*n),np.float64)
        e=np.frombuffer(f.read(n*K),np.uint8).reshape(n,K)
        tp=[]
        for k in range(K):
            v=np.frombuffer(f.read(24),np.int64); g=np.frombuffer(f.read(16),np.float64)
            tp.append(((int(v[0]),int(v[1])),int(v[2]),complex(g[0],g[1])))
    return dict(K=K,n=n,N=N,a=a,b=b,G=Gr+1j*Gi,e=e,tp=tp)
