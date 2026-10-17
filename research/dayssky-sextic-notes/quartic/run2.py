import sys, numpy as np, time
from dualms import run
X=int(sys.argv[1]); H=int(sys.argv[2])
t=time.time()
hs,hN,acc,accr,w2,nn,tb=run(X,H,ks=(1,0,2,-1,3),nrand=3)
print(f"X={X} H={H} rows={len(hs)} cols={nn} W2={w2:.1f} time={time.time()-t:.1f}s")
for Hc in [H//64,H//16,H//4,H]:
    m=hN<=Hc; nh=int(m.sum())
    line=f" H'={Hc:6d} #h={nh:6d} |"
    for k,a in acc.items(): line+=f" k={k}:{np.sum(np.abs(a[m])**2)/(nh*w2):6.3f}"
    line+=" | rand:"+",".join(f"{np.sum(np.abs(a[m])**2)/(nh*w2):.3f}" for a in accr)
    print(line)
for k in (0,1):
    a=acc[k]; r=np.abs(a)**2/w2
    top=np.argsort(-r)[:8]
    print(f" k={k} top rows:", [(hs[i], round(float(r[i]),2)) for i in top])
print(" rand top:", sorted((np.abs(accr[0])**2/w2))[-4:])
