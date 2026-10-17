import numpy as np, sys, math
exec(open('theta_fit.py').read().split("cases=")[0])
Zs=np.logspace(5,6,9)
pN=np.round(pr[:,0]**2+pr[:,1]**2).astype(int)
for k in range(K):
    n=pN[k]
    if n%8!=1 or int(math.isqrt(n))**2==n: continue
    v=np.array([K1([k],Z)[0]/Z**0.75 for Z in Zs]).mean()
    a,b=int(pr[k,0]),int(pr[k,1]); th=math.atan2(b,a)
    q4=any(math.isqrt(n-64*y*y)**2==n-64*y*y for y in range(1,int((n/64)**0.5)+1))
    print("p=%4d pi=(%3d,%3d) |v|=%.4f arg v=%+.3f  arg pi=%+.3f  arg v + arg(pi)/2=%+.3f  arg v - arg(pi)/2=%+.3f  2 quartic res: %s"%(n,a,b,abs(v),np.angle(v),th,(np.angle(v)+th/2+np.pi)%(2*np.pi)-np.pi,(np.angle(v)-th/2+np.pi)%(2*np.pi)-np.pi,q4))
