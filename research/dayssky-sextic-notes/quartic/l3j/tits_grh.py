# Bochner hull for the hyperbolic second-moment series with GRH input on the Mobius side (HEURISTIC model).
# T-side centred coords u; sigma_i(u)_j = u_j - c(j,i) u_i.  C_abs = rho + orthant.  GRH: Mobius side v=sigma3(u) has
# v_i > 0 (1/L holomorphic, bounded by K^eps), v3 > 1/2  ->  T-side cone apex p0=(1/2,1/2,-1/2), rays (1,0,0),(0,1,0),(1,1,-1).
import numpy as np, sys
from scipy.optimize import linprog
C={(1,1):2,(2,2):2,(3,3):2,(1,2):0,(2,1):0,(1,3):-1,(2,3):-1,(3,1):-3,(3,2):-3}
S={}
for i in (1,2,3):
    M=np.eye(3)
    for j in (1,2,3): M[j-1,i-1]-=C[(j,i)]
    S[i]=M
L=int(sys.argv[1]) if len(sys.argv)>1 else 10
elems={tuple(np.eye(3).round(6).flatten()):np.eye(3)}; fr=[np.eye(3)]
for _ in range(L):
    nf=[]
    for g in fr:
        for i in (1,2,3):
            h=S[i]@g; k=tuple(h.round(6).flatten())
            if k not in elems: elems[k]=h; nf.append(h)
    fr=nf
G=list(elems.values()); print("elements",len(G))
rho=np.array([.5,.5,.5]); p0=np.array([.5,.5,-.5])
R_abs=[np.array(r,float) for r in ([1,0,0],[0,1,0],[0,0,1])]
R_grh=[np.array(r,float) for r in ([1,0,0],[0,1,0],[1,1,-1])]
def hull_min(use_grh):
    pts=[g@rho for g in G]; rays=[g@r for g in G for r in R_abs]
    if use_grh:
        pts+=[g@p0 for g in G]; rays+=[g@r for g in G for r in R_grh]
    pts=np.array(pts); rays=np.array(rays); rays=rays/np.linalg.norm(rays,axis=1)[:,None]
    def feas(t):
        n,m=len(pts),len(rays)
        A=np.vstack([np.hstack([pts.T,rays.T]),np.hstack([np.ones(n),np.zeros(m)])])
        b=np.array([0,0,t,1.0])
        return linprog(np.zeros(n+m),A_eq=A,b_eq=b,bounds=(0,None),method='highs').status==0
    lo,hi=-10,10
    if not feas(hi): return None
    for _ in range(40):
        mid=(lo+hi)/2
        if feas(mid): hi=mid
        else: lo=mid
    return hi
print("min u3 at (0,0,u3), absolute convergence only:",hull_min(False))
print("min u3 at (0,0,u3), with GRH region on Mobius side:",hull_min(True))
print("(L3) needs u3 <= 1/2 ; Q^{3/2} corresponds to u3 = 1")
