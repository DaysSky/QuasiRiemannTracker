# Conductor test for the B-row family: if conj(alpha(c)) g~(chi_c) chi_c(m) were GL2-type in c with
# twist conductor ~N(m)^2, then for N(m)^2 << Z the smooth row sum would be negligible (only residues).
# Compare: GL3-type (conductor N(m)^3) predicts generic size ~ sqrt(sum W^2) unless N(m)^3 << Z.
import numpy as np, math, sys
from dualms import build
from gauss4 import chi_prime, gmul
Z=int(sys.argv[1]); H=int(sys.argv[2])
hs,hN,P,tb=build(Z,H)
norms=[p[1] for p in P]; order=sorted(range(len(P)),key=lambda j:norms[j])
ns=[]
def rec(start,pn,lst):
    if pn>Z/4 and lst: ns.append(list(lst))
    for t in range(start,len(order)):
        j=order[t]
        if pn*norms[j]>Z: break
        lst.append(j); rec(t+1,pn*norms[j],lst); lst.pop()
rec(0,1,[])
cache={}
def cpp(i,j):
    if (i,j) not in cache: cache[(i,j)]=chi_prime(P[j][0],P[i][0])
    return cache[(i,j)]
def bump(t):  # C^infinity bump on (1/4,1)
    u=(t-0.25)/0.75
    if u<=0 or u>=1: return 0.0
    return math.exp(-1/(u*(1-u))+4)
nH=len(hs); acc=np.zeros(nH,complex); accU=np.zeros(nH,complex); w2=0
for lst in ns:
    n=(1,0); gam=1.0+0j; Nn=1
    for t,j in enumerate(lst):
        gam*=P[j][2]
        for s_ in lst[:t]: gam*=cpp(s_,j)*cpp(j,s_)
        n=gmul(n,P[j][0]); Nn*=P[j][1]
    W=bump(Nn/Z)
    if W==0: continue
    w2+=W*W
    ex=np.zeros(nH,np.int16); zero=np.zeros(nH,bool)
    for j in lst:
        e=P[j][3]; zero|=(e<0); ex+=np.where(e<0,0,e)
    v=np.where(zero,0,1j**(ex%4))
    al=complex(*n)/math.sqrt(Nn)
    acc+=W*np.conj(al)*gam*v; accU+=W*gam*v
r=np.abs(acc)**2/w2; rU=np.abs(accU)**2/w2
print(f"Z={Z} rows N(m)<={H}, sumW^2={w2:.1f}")
for lo,hi in [(1,4),(5,10),(11,20),(21,40),(41,80),(81,160),(161,H)]:
    m=(hN>=lo)&(hN<=hi)
    if m.sum()==0: continue
    print(f" N(m) in [{lo},{hi}] (#{int(m.sum())}): mean |B|^2/W2 alpha-twisted={r[m].mean():.3f}  untwisted={rU[m].mean():.3f}  N(m)^2/Z<= {hi*hi/Z:.3f}  N(m)^3/Z<= {hi**3/Z:.2f}")
