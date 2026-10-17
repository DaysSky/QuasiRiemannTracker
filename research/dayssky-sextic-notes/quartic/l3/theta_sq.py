import numpy as np, sys, math
exec(open('theta_fit.py').read().split("cases=")[0])
Zs=np.logspace(5,6,9)
pN=np.round(pr[:,0]**2+pr[:,1]**2).astype(int)
out={1:[],5:[]}
for k in range(K):
    n=pN[k]
    if n%4!=1 or int(math.isqrt(n))**2==n: continue
    v=np.array([K1([k,k],Z)[0]/Z**0.75 for Z in Zs]).mean()
    out[n%8].append((n,abs(v)))
for r in (1,5):
    vals=[a for n,a in out[r]]
    print("index pi^2, N(pi)≡%d mod 8: n=%d median %.4f min %.4f max %.4f"%(r,len(vals),np.median(vals),min(vals),max(vals)))
    print("   ", " ".join("%d:%.3f"%(n,a) for n,a in out[r][::2]))
