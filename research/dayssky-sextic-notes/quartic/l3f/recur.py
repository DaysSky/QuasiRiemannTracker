# Recursive S52: table s(q) for S(Q,Y) <= Y^{s(q)}, Q = Y^q, T*-family (fix (a)).
# Block sigma uses min(S49, current table); hull = 12 chambers with block regions (RMS + pointwise twist exponents).
import numpy as np, json, sys
from scipy.optimize import linprog
QS=np.round(np.arange(0.40,1.61,0.05),3)
def s49(q):   # S49 Theorem-2 shape for Y <= Q^{5/3}, i.e. q >= 0.6:  (Y + Q Y^{1/3} + Q^{1/2} Y^{3/4}) Y ; else reflect
    return 1+max(1,q+1/3,q/2+0.75)
def s_reflect(q,tab):   # q < 2/3: Y > Q^{3/2}; S(Q,Y) ~ (Y/Y') S(Q,Y'), Y' = Q^3/Y -> q' = q/(3q-1)
    if q<=1/3+1e-9: return 1.0+q    # beyond conductor: only diagonal-size (length >> conductor)
    qp=q/(3*q-1); Yp=3*q-1           # log_Y Y'
    return 1-Yp+Yp*S(qp,tab)
def S(q,tab):
    if q<2/3-1e-9: return s_reflect(q,tab)
    v=s49(q)
    if tab is not None: v=min(v,np.interp(q,QS,tab))
    return v
def sigma(tC,tD,tab):
    if tD<=0: return tC if tC<=0 else -1e9
    if tC>3*tD+1e-12: return -1e9
    cands=[tC+2.5*tD]
    q=tD/tC if tC>0 else 9
    if tC>0: cands.append(tC*S(q,tab))
    return min(cands)
def blk(x,y,tC,tD,tab,tK=None,rms=False):
    if tK is None:
        B1=-x*tC-y*tD+0.5*max(tC,tD,2/3*(tC+tD))
        B2=-y*tD-(0.5+x)*tC+0.5*sigma(tC,tD,tab)
    elif rms:
        B1=-x*tC-y*tD+0.5*(max(tK+tC,tD,2/3*(tK+tC+tD))-tK)
        B2=-y*tD-(0.5+x)*tC+0.5*sigma(tC,tD,tab)-0.5*min(tD,tK)
    else:
        B1=-x*tC-y*tD+0.5*max(tC,tD,2/3*(tC+tD))
        B2=-y*tD-(0.5+x)*tC+0.5*sigma(tC,tD,tab)
    return min(B1,B2)
def okd(x,y,tab): return all(blk(x,y,1.0,dl,tab)<-1e-9 for dl in np.linspace(0,0.5,61))
def okc(x,y,tab): return all(blk(x,y,r,1.0,tab)<-1e-9 for r in np.linspace(0,1.5,61))
def bnd(f,fixed,tab,lo,hi,which):
    if which=='d':
        if not okd(hi,fixed,tab): return None
        for _ in range(30):
            m=(lo+hi)/2
            if okd(m,fixed,tab): hi=m
            else: lo=m
    else:
        if not okc(fixed,hi,tab): return None
        for _ in range(30):
            m=(lo+hi)/2
            if okc(fixed,m,tab): hi=m
            else: lo=m
    return hi
def bexp(x,y,tab,rms):
    best=-1e9
    for tC in np.linspace(0,20,25):
        for tD in np.linspace(0,(1+tC)/2,25):
            best=max(best,blk(x,y,tC,tD,tab,tK=1.0,rms=rms))
    return max(best,0.0)
sc=np.array([[-1,0],[3,1]]); sd=np.array([[1,1],[0,-1]])
ch={(1,0,0,1):(np.eye(2,dtype=int),np.zeros(2))}; fr=[(1,0,0,1)]
while fr:
    nf=[]
    for k in fr:
        M,c=ch[k]
        for nm,g in (("c",sc),("d",sd)):
            N=g@M; kk=tuple(N.flatten()); cn=c-M[1] if nm=="d" else c.copy()
            if kk not in ch: ch[kk]=(N,cn); nf.append(kk)
    fr=nf
def local_points(lam,tab):
    L0=13/8
    def Lam(x):
        if x>0: return max(1.0,L0-3*x)
        if x==0: return lam
        return max(1-6*x,L0-3*x,lam-6*x)
    g=lambda x: Lam(x)/2 if x<=0.5 else 0.5
    e=1e-6; loc=[]
    for yp in np.linspace(-3,0,61): loc.append(((0.5-yp+e,yp),-yp,1/6-yp))
    for xp in np.linspace(-1.5,1.5,61): loc.append(((xp,g(xp)+e),0.0,0.0))
    for yp in np.linspace(0,1.6,17):
        xb=bnd(None,yp,tab,-2,2,'d')
        if xb is None: continue
        for dx in (1e-3,0.05):
            loc.append(((xb+dx,yp),bexp(xb+dx,yp,tab,True),bexp(xb+dx,yp,tab,False)))
    for xp in np.linspace(0,0.5,11):
        yb=bnd(None,xp,tab,-2,4,'c')
        if yb is not None: loc.append(((xp,yb+1e-3),0.0,0.0))
    pts=[]
    for k,(M,c) in ch.items():
        Mi=np.linalg.inv(M)
        for (vp,br,bp) in loc:
            v=Mi@np.array(vp); f=c@v; pts.append((v[0],v[1],br+f,bp+f))
    return np.array(pts)
def lpmin(P,obj,ylo,yhi):
    n=len(P); A_ub=np.vstack([P[:,1],-P[:,1],-(3*P[:,0]+P[:,1]),-(2*P[:,0]+P[:,1]),-(P[:,0]+P[:,1])])
    r=linprog(obj,A_ub=A_ub,b_ub=[yhi,-ylo,0,0,0],A_eq=np.ones((1,n)),b_eq=[1],bounds=(0,None),method="highs")
    return r.fun if r.status==0 else 9e9
def s52(P,q,n=13):
    # log_Y units, Q = Y^q, H = 2q-1 twists; X main = (3/4)q + 1/2 + x ; dual = 1/2 + q(1/2+y) + x
    H=2*q-1
    if H<=0: return 2.0
    def X(eta,eps):
        o=P[:,0]+eta*P[:,2]+2*eps*P[:,3]
        main=0.75*q+0.5+lpmin(P,o,0.25,0.25)
        dual=0.5+0.5*q+lpmin(P,o+q*P[:,1],0.0,0.2499)
        nos=0.5+0.5*q+lpmin(P,o+q*P[:,1],0.2501,5.0)
        return min(nos,max(main,dual))
    etas=np.linspace(0,H,n)
    eff=[max(X(g,t) for t in np.linspace(0,(H-g)/4,4)) for g in etas]
    best=9
    for i,k0 in enumerate(etas):
        offT=max((1-q)+etas[j]+(H-etas[j])/4+2*eff[j] for j in range(i+1))
        best=min(best,max(2+(H-k0)/4,(1+offT)/2,2))
    return best
tab=None; lam=1.5234
for it in range(5):
    P=local_points(lam,tab)
    new=[]
    for q in QS:
        v=S(q,tab)
        if 0.66<=q<=1.2: v=min(v,s52(P,q))
        new.append(v)
    tab=np.array(new)
    lam_new=1.5*(np.interp(2/3,QS,tab)-1)
    print("iter %d: lam=%.4f -> (L3)<<Q^%.4f ; s(q) at q=2/3,0.8,1,1.2: %s"%(it,lam,lam_new,np.round([np.interp(qq,QS,tab) for qq in (2/3,0.8,1.0,1.2)],4)),flush=True)
    json.dump({"QS":QS.tolist(),"tab":tab.tolist()},open('recur_tab.json','w'))
    if abs(lam_new-lam)<2e-4: break
    lam=min(lam,lam_new)
