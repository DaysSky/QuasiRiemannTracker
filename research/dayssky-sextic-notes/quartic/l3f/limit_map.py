# Delta_n -> Delta_{n+1}: S50 bookkeeping at every q in [2/3,5/6] using only the explicit two-chamber hull point at the
# Region_d boundary point (x',0), x'=1/3+dp, and the cone bound b_loc(tK,tE) <= c1*Delta*(tK+tE) + c2*tE (+ margin).
import numpy as np
def lXi(a,b): return max(a,b,2/3*(a+b))
def X(q,eta,eps,Delta,dp,c1=3.0,c2=0.5,marg=0.0):
    tK,tE=eta,2*eps
    bloc=c1*Delta*(tK+tE)+c2*tE+marg
    return 0.75*q+5/12+(3*(1/3+dp)-0.25)*(tK+tE)+bloc          # main term; dual <= main on the same combination line
def feasible(q,s,Delta,dp,n_eta=401,n_eps=41,**kw):
    H=2*q-1
    def gen(g): return 2-q+(H-g)/4+lXi(g,q)
    def tail_ok(g,e): return 2-q+(H-g)/4-e+lXi(g,q-e)<=s+1e-12
    etas=np.linspace(0,H,n_eta)
    # generic cost is convex in eta (max of affine): sup over [k0,H] = max(gen(k0),gen(H))
    ok_mds=[]
    for g in etas:
        emax=(H-g)/4
        # e-tail: need tail_ok for all e in [e0,emax]; tail cost is piecewise affine in e -> check grid incl. endpoints
        es=np.linspace(0,emax,n_eps)
        e0=emax
        for t in es:
            if all(tail_ok(g,u) for u in np.linspace(t,emax,n_eps)): e0=t; break
        xe=max(X(q,g,e,Delta,dp,**kw) for e in np.linspace(0,e0,n_eps))
        ok_mds.append((1+(1-q)+g+(H-g)/4+2*xe)/2<=s+1e-12)
    ok_mds=np.array(ok_mds)
    # k0 continuous: MDS on [0,k0] (grid), generic on [k0,H]
    for i,k0 in enumerate(etas):
        if not ok_mds[:i+1].all(): break
        if max(gen(k0),gen(H))<=s+1e-12: return True
    return False
def s_min(q,Delta,dp,**kw):
    lo,hi=2.0,3.0
    for _ in range(40):
        m=(lo+hi)/2
        if feasible(q,m,Delta,dp,**kw): hi=m
        else: lo=m
    return hi
def step(Delta,dp,qs=np.linspace(2/3,5/6,21),**kw):
    ex=[s_min(q,Delta,dp,**kw)-(q+4/3) for q in qs]
    return max(ex),ex
if __name__=="__main__":
    dp=1e-4
    D=1/12
    for n in range(8):
        Dn,ex=step(D,dp)
        print("n=%d Delta=%.6f -> %.6f   excess at q=2/3,0.7,0.75,0.8,5/6: %s"%(n,D,Dn,np.round([ex[0],ex[4],ex[10],ex[16],ex[20]],6)),flush=True)
        D=max(Dn,0.0)
