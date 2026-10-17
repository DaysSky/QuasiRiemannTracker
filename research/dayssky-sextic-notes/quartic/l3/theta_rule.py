import numpy as np, sys
exec(open('theta_fit.py').read().split("cases=")[0])
Zs=np.logspace(5,6,9)
pN=np.round(pr[:,0]**2+pr[:,1]**2).astype(int)
rows=[]
for k in range(K):
    n=pN[k]
    if n%4!=1: continue   # split primes only (inert have N=q^2)
    v=np.array([K1([k],Z)[0]/Z**0.75 for Z in Zs]); m=v.mean(); sd=v.std()/np.sqrt(len(v))
    rows.append((n,n%8,abs(m),sd,(int(pr[k,0]),int(pr[k,1]))))
rows.sort()
for n,r,am,sd,p in rows: print("N=%4d  N mod 8=%d  |mean K1/Z^{3/4}| = %.4f   (sd of mean ~%.4f)  pi=%s"%(n,r,am,sd,p))
a1=[am for n,r,am,sd,p in rows if r==1]; a5=[am for n,r,am,sd,p in rows if r==5]
print("N≡1 mod 8: n=%d  median %.4f  min %.4f | N≡5 mod 8: n=%d median %.4f max %.4f"%(len(a1),np.median(a1),min(a1),len(a5),np.median(a5),max(a5)))
