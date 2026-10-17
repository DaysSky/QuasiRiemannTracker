# Recursive S52 (fix (a)), corrected: quartic large sieve Xi(A,B)=A+B+(AB)^{2/3} everywhere; joint twist exponent b(tK,tE)
# (tK = RMS twist h0 ~ Y^eta, tE = pointwise twist e^2 ~ Y^{2 eps}); S(q) table recursion.
import numpy as np, json, sys
from scipy.optimize import linprog
QS=np.array(sorted(set(np.round(np.arange(0.40,1.61,0.05),3).tolist()+[round(2/3,6),0.68,0.72])))
def s49(q): return 1+max(1,q+1/3,q/2+0.75)          # S(Q,Y) exponent (log_Y), Q=Y^q, valid q>=0.6
def S(q,tab):
    if q<2/3-1e-9:
        if q<=1/3+1e-9: return 1.0+q
        qp=q/(3*q-1); Yp=3*q-1
        return 1-Yp+Yp*S(qp,tab)
    v=s49(q)
    if tab is not None: v=min(v,float(np.interp(q,QS,tab)))
    return v
def sigma_vec(tC,tD,tab):
    out=np.full(tC.shape,-1e9)
    for idx in np.ndindex(tC.shape):
        c=tC[idx]; d=tD[idx]
        if d<=0: out[idx]= (0.0 if c<=0 else -1e9); continue
        if c>3*d+1e-12: continue
        v=c+2.5*d
        if c>0: v=min(v,c*S(d/c,tab))
        out[idx]=v
    return out
def blocks_val(x,y,tC,tD,tK,sig):
    B1=-x*tC-y*tD+0.5*(np.maximum.reduce([tK+tC,tD,2/3*(tK+tC+tD)])-tK)
    B2=-y*tD-(0.5+x)*tC+0.5*sig+0.5*(np.maximum.reduce([np.full(tD.shape,tK),tD,2/3*(tK+tD)])-tK-tD)
    return np.minimum(B1,B2)
# precomputed grids
GC=np.concatenate([np.linspace(0,4,81),np.linspace(4.1,40,60)])
def b_twist(x,y,tK,tE,tab,cache={}):
    key=(round(x,6),round(y,6),round(tK,6),round(tE,6),id(tab))
    if key in cache: return cache[key]
    best=0.0
    for c in GC:
        tD=np.linspace(0,(tK+tE+c)/2,81); tC=np.full(tD.shape,c)
        sig=sigma_vec(tC,tD,tab)
        best=max(best,float(blocks_val(x,y,tC,tD,tK,sig).max()))
    cache[key]=best; return best
def ok_dir(x,y,tC,tD,tab):
    sig=sigma_vec(tC,tD,tab); return bool((blocks_val(x,y,tC,tD,0.0,sig)<-1e-9).all())
DL=np.linspace(0,0.5,61); RH=np.linspace(0,1.5,61)
def okd(x,y,tab): return ok_dir(x,y,np.ones_like(DL),DL,tab)
def okc(x,y,tab): return ok_dir(x,y,RH,np.ones_like(RH),tab)
def bis(f,lo,hi):
    if not f(hi): return None
    for _ in range(28):
        m=(lo+hi)/2
        if f(m): hi=m
        else: lo=m
    return hi
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
    e=1e-6; loc=[]   # entries: ((x',y'), kind, extra)
    for yp in np.linspace(-3,0,61): loc.append(((0.5-yp+e,yp),'ldfe',-yp))
    for yp in np.linspace(0,3,31): loc.append(((0.5+e,yp),'ld',0.0))
    for xp in np.linspace(-1.5,1.5,61): loc.append(((xp,g(xp)+e),'zero',0.0))
    for yp in np.linspace(0,1.6,17):
        xb=bis(lambda xx: okd(xx,yp,tab),-2,2)
        if xb is None: continue
        for dx in (1e-3,0.05):
            loc.append(((xb+dx,yp),'blk',0.0))
    for xp in np.linspace(0,0.5,11):
        yb=bis(lambda yy: okc(xp,yy,tab),-2,4)
        if yb is not None: loc.append(((xp,yb+1e-3),'zero',0.0))
    return loc
def local_b(item,tK,tE,tab):
    (vp,kind,ex)=item
    if kind=='zero': return 0.0
    if kind=='ld': return tE/6
    if kind=='ldfe': return tE/6+ex*(tK+tE)
    return b_twist(vp[0],vp[1],tK,tE,tab)
def points(loc,tK,tE,tab):
    lb=[local_b(it,tK,tE,tab) for it in loc]; pts=[]
    for k,(M,c) in ch.items():
        Mi=np.linalg.inv(M)
        for it,b in zip(loc,lb):
            v=Mi@np.array(it[0]); pts.append((v[0],v[1],b+(c@v)*(tK+tE)))
    return np.array(pts)
def lpmin(P,obj,ylo,yhi):
    n=len(P); A_ub=np.vstack([P[:,1],-P[:,1],-(3*P[:,0]+P[:,1]),-(2*P[:,0]+P[:,1]),-(P[:,0]+P[:,1])])
    r=linprog(obj,A_ub=A_ub,b_ub=[yhi,-ylo,0,0,0],A_eq=np.ones((1,n)),b_eq=[1],bounds=(0,None),method="highs")
    return r.fun if r.status==0 else 9e9
def Xq(loc,q,eta,eps,tab):
    P=points(loc,eta,2*eps,tab); o=P[:,0]+P[:,2]
    main=0.75*q+0.5+lpmin(P,o,0.25,0.25)
    dual=0.5+0.5*q+lpmin(P,o+q*P[:,1],0.0,0.2499)
    nos=0.5+0.5*q+lpmin(P,o+q*P[:,1],0.2501,5.0)
    return min(nos,max(main,dual))
def s52(loc,q,tab,n=17):
    H=2*q-1
    if H<=0: return 2.0
    etas=np.linspace(0,H,n)
    eff=[max(Xq(loc,q,g,t,tab) for t in np.linspace(0,(H-g)/4,5)) for g in etas]
    best=9
    for i,k0 in enumerate(etas):
        offT=max((1-q)+etas[j]+(H-etas[j])/4+2*eff[j] for j in range(i+1))
        best=min(best,max(2+(H-k0)/4,(1+offT)/2,2))
    return best,eff
if __name__=="__main__":
    d=json.load(open('recur2_tab.json')); tab=np.interp(QS,np.array(d['QS']),np.array(d['tab'])); lam=d['lam']
    for it in range(2):
        loc=local_points(lam,tab)
        new=[]; info=""
        for q in QS:
            v=S(q,tab)
            if q in (round(2/3,6),0.7,0.8,0.9,1.0):
                r,eff=s52(loc,q,tab); v=min(v,r)
                if abs(q-2/3)<1e-3: info=str(np.round(eff,4))
            new.append(v)
        tab=np.array(new)
        lam_new=1.5*(float(np.interp(2/3,QS,tab))-1)
        print("iter %d: lam %.4f -> (L3)<<Q^%.4f | s(2/3,0.8,1,1.2)=%s | X_eff(q=2/3) %s"%(it,lam,lam_new,np.round([np.interp(qq,QS,tab) for qq in (2/3,0.8,1.0,1.2)],4),info),flush=True)
        json.dump({"QS":QS.tolist(),"tab":tab.tolist(),"lam":lam_new},open('recur2_fine.json','w'))
        if abs(lam_new-lam)<2e-4: break
        lam=min(lam,lam_new)
