# Correct two-series ("groupoid") Bochner model for the (L^3) second moment.
# Z_T(u): T-side series, centred coords u=(s1-1/2, s2-1/2, w-1/2). FEs: Kubota sigma_i: u_i->-u_i, u3->u3+3u_i (i=1,2).
# Z_mu(v): Mobius side sum_k N(k)^{-1/2-v3} / (L(1/2+v1,xi_k) L(1/2+v2,conj xi_k)). Relation Z_T(u)=gamma*Z_mu(v), v=P(u)=(u1+u3,u2+u3,-u3).
# FE of Z_mu: J: (v1,v2,v3) -> (-v2,-v1, v3-v1-v2)  (both 1/L reflected; root numbers cancel).
# Regions: T abs: u_i>1/2,u3>1/2.  mu abs: v_i>1/2, v3>1/2.  mu GRH: v_i>0, v3>1/2.
# Each region = apex + cone(rays).  Iterate: H_T = conv(W_T (H_T u P^{-1} H_mu)), H_mu = conv(J-orbit(H_mu u P H_T)),
# represented by finite generator sets (points, rays).  Report min u3 with (0,0,u3) in H_T.
import numpy as np, sys
from scipy.optimize import linprog
def sig(i):
    M=np.eye(3); M[i,i]=-1; M[2,i]=3; return M
S1,S2=sig(0),sig(1)
WT=[np.eye(3),S1,S2,S1@S2]
P=np.array([[1,0,1],[0,1,1],[0,0,-1.]]); Pinv=np.linalg.inv(P)
J=np.array([[0,-1,0],[-1,0,0],[-1,-1,1.]])
WM=[np.eye(3),J]
def cone(apex,rays): return ([np.array(apex,float)],[np.array(r,float) for r in rays])
E=[[1,0,0],[0,1,0],[0,0,1]]
T_abs=cone([.5,.5,.5],E)
M_abs=cone([.5,.5,.5],E)
M_grh=cone([0,0,.5],E)
def closure(Tgen,Mgen,iters=6):
    for it in range(iters):
        tp,tr=Tgen; mp,mr=Mgen
        # transfer
        tp2=tp+[Pinv@p for p in mp]; tr2=tr+[Pinv@r for r in mr]
        mp2=mp+[P@p for p in tp]; mr2=mr+[P@r for r in tr]
        # group orbits
        tp3=[g@p for g in WT for p in tp2]; tr3=[g@r for g in WT for r in tr2]
        mp3=[g@p for g in WM for p in mp2]; mr3=[g@r for g in WM for r in mr2]
        def dedup(L):
            out=[];seen=set()
            for x in L:
                k=tuple(np.round(x,9))
                if k not in seen: seen.add(k); out.append(x)
            return out
        Tgen=(dedup(tp3),dedup(tr3)); Mgen=(dedup(mp3),dedup(mr3))
    return Tgen,Mgen
def min_u3(Tgen):
    pts=np.array(Tgen[0]); rays=np.array(Tgen[1]); rays=rays/np.linalg.norm(rays,axis=1)[:,None]
    def feas(t):
        n,m=len(pts),len(rays)
        A=np.vstack([np.hstack([pts.T,rays.T]),np.hstack([np.ones(n),np.zeros(m)])])
        return linprog(np.zeros(n+m),A_eq=A,b_eq=np.array([0,0,t,1.]),bounds=(0,None),method='highs').status==0
    lo,hi=-10,20
    if not feas(hi): return None
    for _ in range(45):
        mid=(lo+hi)/2
        if feas(mid): hi=mid
        else: lo=mid
    return hi
for name,Mstart in [("unconditional",M_abs),("GRH on Mobius side",(M_abs[0]+M_grh[0],M_abs[1]+M_grh[1]))]:
    Tg,Mg=closure(T_abs,Mstart)
    print(name,": generators",len(Tg[0]),len(Tg[1]),"  min u3 =",round(min_u3(Tg),6))
print("Q^{3/2} <-> u3=1 ; (L3) <-> u3=1/2")
