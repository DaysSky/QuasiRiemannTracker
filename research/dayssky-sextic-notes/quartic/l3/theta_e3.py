import numpy as np, sys, math
exec(open('theta_fit.py').read().split("cases=")[0])
Zs=np.logspace(5,6,9)
pN=np.round(pr[:,0]**2+pr[:,1]**2).astype(int)
out={(e,r):[] for e in (2,3) for r in (1,5)}
for k in range(K):
    n=pN[k]
    if n%4!=1 or math.isqrt(n)**2==n or n>400: continue
    for e in (2,3):
        v=np.array([K1([k]*e,Z)[0]/Z**0.75 for Z in Zs]).mean(); out[(e,n%8)].append(abs(v))
for (e,r),v in sorted(out.items()): print("e=%d N≡%d mod 8: n=%d median %.4f min %.4f max %.4f"%(e,r,len(v),np.median(v),min(v),max(v)))
