import numpy as np, sys, math
exec(open('theta_fit.py').read().split("cases=")[0])
Zs=np.logspace(5,6,9)
pN=np.round(pr[:,0]**2+pr[:,1]**2).astype(int)
split=[k for k in range(K) if pN[k]%4==1 and math.isqrt(pN[k])**2!=pN[k] and pN[k]<=200]
rng=np.random.default_rng(5); mob=[]; blind=[]
for _ in range(60):
    a,b=rng.choice(split,2,replace=True)
    idx=[a,b]
    vm=np.array([K1(idx,Z,mob=True)[0]/Z**0.75 for Z in Zs]).mean(); vb=np.array([K1(idx,Z)[0]/Z**0.75 for Z in Zs]).mean()
    mob.append(abs(vm)); blind.append(abs(vb))
print("composite indices (60 random products of two primes, N<=200 each): Moebius-weighted median %.4f max %.4f | parity-blind median %.4f max %.4f"%(np.median(mob),max(mob),np.median(blind),max(blind)))
