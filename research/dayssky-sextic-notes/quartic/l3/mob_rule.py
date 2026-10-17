import numpy as np, sys, math
exec(open('theta_fit.py').read().split("cases=")[0])
Zs=np.logspace(5,6,9)
pN=np.round(pr[:,0]**2+pr[:,1]**2).astype(int)
res={1:[],5:[]}
for k in range(K):
    n=pN[k]
    if n%4!=1 or math.isqrt(n)**2==n: continue
    v=np.array([K1([k],Z,mob=True)[0]/Z**0.75 for Z in Zs]).mean()
    res[n%8].append(abs(v))
for r in (1,5): print("Moebius-weighted (alpha-twisted family) index pi, N≡%d mod 8: n=%d median %.4f max %.4f"%(r,len(res[r]),np.median(res[r]),max(res[r])))
