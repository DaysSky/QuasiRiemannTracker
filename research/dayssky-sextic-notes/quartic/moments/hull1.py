# First-moment double series Z1(s,w)=sum_m psi(m^3,s) N(m)^{-w} = sum_c a_c L(w,conj chi_c) N(c)^{-s}.
# Cartan data [[2,-1],[-3,2]] (G2), but tau_w lands on sum_c mu(c) L(w',chi_c) N c^{-u} (Mobius) -> only <tau_s> acts.
import numpy as np
from scipy.optimize import linprog
TS=np.array([[-1,0],[3,1]],float); TW=np.array([[1,1],[0,-1]],float); E=np.eye(2)
def P(v,r): return ([np.array(x,float) for x in v],[np.array(x,float) for x in r])
RA=P([(.5,.5)],[E[0],E[1]])
Om=P([(.5,.5),(-.5,2)],[E[0],E[1],(-1,3)])          # y > 1/2 + 3 e(x)
Oc=P([(.5,.5),(1,-.5)],[E[0],E[1],(1,-1)])            # x > 1/2 + e(y)
def Rpm(s0): return P([(s0-.5,.5)],[E[0],E[1]])
def run(base, group):
    V=[];R=[]
    for g in group:
        for (vs,rs) in base: V+=[g@v for v in vs]; R+=[g@r for r in rs]
    n,m=len(V),len(R); c=np.array([v[1] for v in V]+[r[1] for r in R])
    A=np.zeros((2,n+m)); b=np.array([0,1.])
    for j,v in enumerate(V): A[0,j]=v[0]; A[1,j]=1
    for j,r in enumerate(R): A[0,n+j]=r[0]
    res=linprog(c,A_eq=A,b_eq=b,bounds=[(0,None)]*(n+m),method='highs'); return res.fun if res.status==0 else res.status
img=lambda M,Pp:([M@v for v in Pp[0]],[M@r for r in Pp[1]])
for s0 in (1.0,0.75,0.5):
    base=[RA,Om,Oc,img(TW,Oc),img(TW,Rpm(s0))]
    y=run(base,[np.eye(2),TS])
    print(f"first moment, zero-free Re u>{s0}: min y at x=0: {y:.4f} -> error exponent {y+.5:.4f}")
# hypothetical full G2 group (12 elements) for comparison
G=[np.eye(2)]
for _ in range(12):
    G=G+[s@g for g in G for s in (TS,TW)]
    G=list({tuple(np.round(g,9).ravel()):g for g in G}.values())
print("|W(G2)| =",len(G))
y=run([RA,Om,Oc],G); print("hypothetical full G2 action:",y)
