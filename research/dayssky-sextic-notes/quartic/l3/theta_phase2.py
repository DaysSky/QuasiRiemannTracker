import numpy as np, sys, math
exec(open('theta_fit.py').read().split("cases=")[0])
Zs=np.logspace(5,6,9)
pN=np.round(pr[:,0]**2+pr[:,1]**2).astype(int)
A=d['a'].astype(int); Bv=d['b'].astype(int)
for k in range(K):
    n=pN[k]
    if n%8!=1 or int(math.isqrt(n))**2==n: continue
    v=np.array([K1([k],Z)[0]/Z**0.75 for Z in Zs]).mean()
    a,b=int(pr[k,0]),int(pr[k,1])
    j=np.where((A==a)&(Bv==b))[0]
    if len(j)==0: print("no c for",a,b); continue
    gt=g[j[0]]   # gtilde(conj chi_pi)
    al=complex(a,b)/math.sqrt(n)
    out=[]
    for name,z in [("v/g",v/gt),("v*g",v*gt),("v/conj g",v/np.conj(gt)),("v*al^.5/g",v*np.sqrt(al)/gt),("v/(g*al^-.5)",v/(gt*np.sqrt(np.conj(al))))]:
        out.append("%s arg %+.2f"%(name,np.angle(z)))
    q4=any(math.isqrt(n-64*y*y)**2==n-64*y*y for y in range(1,int((n/64)**0.5)+1))
    print("p=%4d (%3d,%3d) |v|=%.4f  "%(n,a,b,abs(v))+" | ".join(out)+"  2q4=%s"%q4)
