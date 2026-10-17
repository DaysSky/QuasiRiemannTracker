import numpy as np
d=np.fromfile('pg1e6.bin').reshape(-1,5); a,b,p,Gr,Gi=d.T
g=(Gr+1j*Gi)/np.sqrt(p); th=np.angle(a+1j*b); n=len(g)
res=[]
cls4=((a%4).astype(int),(b%4).astype(int))
for k in range(-16,17):
    v=g*np.exp(1j*k*th)
    res.append((abs(v.sum())/np.sqrt(n),'k=%d all'%k, v.sum()))
    for c in [(1,0),(3,2)]:
        s=(cls4[0]==c[0])&(cls4[1]==c[1]); res.append((abs(v[s].sum())/np.sqrt(s.sum()),'k=%d cls %s'%(k,c), v[s].sum()))
    for m in [3,5]:
        for x in range(m):
            for y in range(m):
                s=((a%m)==x)&((b%m)==y)
                if s.sum()>500: res.append((abs(v[s].sum())/np.sqrt(s.sum()),'k=%d mod%d (%d,%d)'%(k,m,x,y), v[s].sum()))
res.sort(key=lambda r:-r[0])
for r in res[:20]: print("%.2f %s %s"%(r[0],r[1],np.round(r[2],1)))
print("tests",len(res))
# growth of k=0 and k=1 partial sums
for X in [1e5,3e5,1e6]:
    s=p<=X; print(X, s.sum(), "k=0:",np.round(g[s].sum(),1), " k=1 (T):",np.round((g*np.exp(1j*th))[s].sum(),1), " X^{3/4}/logX=%.0f"%(X**0.75/np.log(X)))
