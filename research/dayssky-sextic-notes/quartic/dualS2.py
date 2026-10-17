# Mean square over rows h (N h <= U) of S(h)=sum_{c sf primary, c≡1 mod 4, N c~Z} w(c) chi_c(h)
# for w = mu(c) (the Poisson-dual coefficient), w = conj(alpha)*gtilde (the B-row), w = random phases.
import numpy as np, math, sys, time
from dualms import build
from gauss4 import chi_prime, gmul
Z=int(sys.argv[1]); U=int(sys.argv[2]); cls=sys.argv[3] if len(sys.argv)>3 else 'n1'
hs,hN,P,tb=build(Z,U)
norms=[p[1] for p in P]; order=sorted(range(len(P)),key=lambda j:norms[j])
ns=[]
def rec(start,pn,lst):
    if pn>Z/2 and lst: ns.append(list(lst))
    for t in range(start,len(order)):
        j=order[t]
        if pn*norms[j]>Z: break
        lst.append(j); rec(t+1,pn*norms[j],lst); lst.pop()
rec(0,1,[])
cache={}
def cpp(i,j):
    if (i,j) not in cache: cache[(i,j)]=chi_prime(P[j][0],P[i][0])
    return cache[(i,j)]
nH=len(hs); rng=np.random.default_rng(0)
acc={k:np.zeros(nH,complex) for k in ['mu','B','rnd','mu_conj']}
w2=0
for lst in ns:
    n=(1,0); gam=1.0+0j; Nn=1
    for t,j in enumerate(lst):
        gam*=P[j][2]
        for s_ in lst[:t]: gam*=cpp(s_,j)*cpp(j,s_)
        n=gmul(n,P[j][0]); Nn*=P[j][1]
    c4=1 if (n[0]%4,n[1]%4)==(1,0) else 3
    if cls=='n1' and c4!=1: continue
    if cls=='n3' and c4!=3: continue
    tt=Nn/Z; W=math.sin(math.pi*(2*tt-1))**2 if 0.5<tt<=1 else 0.0
    if W==0: continue
    w2+=W*W
    al=complex(*n)/math.sqrt(Nn); mu=(-1)**len(lst)
    ex=np.zeros(nH,np.int16); zero=np.zeros(nH,bool)
    for j in lst:
        e=P[j][3]; zero|=(e<0); ex+=np.where(e<0,0,e)
    v=np.where(zero,0,1j**(ex%4))
    acc['mu']+=W*mu*v; acc['mu_conj']+=W*mu*np.conj(v)
    acc['B']+=W*np.conj(al)*gam*v
    acc['rnd']+=W*np.exp(2j*np.pi*rng.random())*v
print(f"Z={Z} U={U} cls={cls} cols={len(ns)} W2={w2:.1f} rows={nH}")
for Uc in [U//64,U//16,U//4,U]:
    m=(hN<=Uc); nh=int(m.sum())
    print(f" U'={Uc:7d} U'/Z={Uc/Z:6.2f} #h={nh:7d} | "+"  ".join(f"{k}:{np.sum(np.abs(a[m])**2)/(nh*w2):7.3f}" for k,a in acc.items()))
