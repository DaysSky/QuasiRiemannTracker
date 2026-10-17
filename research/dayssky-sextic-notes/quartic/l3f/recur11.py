# recur9: recur7 logic, vectorised, with three fixes found in the audit:
#  (1) the recursion starts from the PROVED lambda = 13/8 (S49), not from an unproved value;
#  (2) S(q) uses the table only inside [QS[0], QS[-1]]; outside it uses s49 (np.interp clamps, which under-estimated S for q>1.6);
#  (3) FE-reflected R_d points ('ldfe') in polar chambers carry the polar exponent of their image point when the image has y'<1/4.
import numpy as np, json, sys
from scipy.optimize import linprog
QS=np.array([2/3]+sorted(set(np.round(np.arange(0.67,0.905,0.01),3).tolist()+np.round(np.arange(0.95,1.61,0.05),3).tolist())))
def s49v(q): return 1+np.maximum.reduce([np.ones_like(q),q+1/3,q/2+0.75])
def S_high(q,tab):
    v=s49v(q)
    if tab is not None:
        inr=(q>=QS[0]-1e-9)&(q<=QS[-1])
        j=np.clip(np.searchsorted(QS,q-1e-9,side='left'),0,len(QS)-1)   # smallest grid point >= q (step-up: S is monotone in Q)
        v=np.where(inr,np.minimum(v,tab[j]),v)
    return v
def Sv(q,tab):
    q=np.asarray(q,float); out=np.empty_like(q)
    low=q<=1/3+1e-9; mid=(~low)&(q<2/3-1e-9); high=q>=2/3-1e-9
    out[low]=1.0+q[low]
    if mid.any():
        Yp=3*q[mid]-1; qp=q[mid]/Yp
        out[mid]=1-Yp+Yp*S_high(qp,tab)
    out[high]=S_high(q[high],tab)
    return out
def S(q,tab): return float(Sv(np.array([q]),tab)[0])
def sigma_vec(tC,tD,tab):
    tC=np.asarray(tC,float); tD=np.asarray(tD,float)
    out=np.full(tC.shape,-1e9)
    zero=tD<=0
    out[zero&(tC<=0)]=0.0
    act=(~zero)&(tC<=3*tD+1e-12)
    v=tC+2.5*tD
    pos=act&(tC>0)
    if pos.any():
        r=np.ones_like(tC); r[pos]=tD[pos]/tC[pos]
        sv=np.full(tC.shape,np.inf); sv[pos]=tC[pos]*Sv(r[pos],tab)
        v=np.minimum(v,sv)
    out[act]=v[act]
    return out
def blocks_val(x,y,tC,tD,tK,sig):
    B1=-x*tC-y*tD+0.5*(np.maximum.reduce([tK+tC,tD,2/3*(tK+tC+tD)])-tK)
    B2=-y*tD-(0.5+x)*tC+0.5*sig+0.5*(np.maximum.reduce([np.full(tD.shape,float(tK)),tD,2/3*(tK+tD)])-tK-tD)
    return np.minimum(B1,B2)
GC=np.linspace(0,20,41)
U=np.linspace(0,1,31)
def b_twist(x,y,tK,tE,tab,cache,th=0.5,polar=False):
    key=(round(x,6),round(y,6),round(tK,6),round(tE,6),round(th,4),polar)
    if key in cache: return cache[key]
    C=GC[:,None]*np.ones_like(U)[None,:]
    best=0.0
    tD=th*(tK+tE+C)*U[None,:]
    best=max(best,float(blocks_val(x,y,C,tD,tK,sigma_vec(C,tD,tab)).max()))
    if th<1:
        tD2=(1-th)*(tK+tE+C)*U[None,:]
        best=max(best,float(blocks_val(x+y,-y,C,tD2,tK,sigma_vec(C,tD2,tab)).max())-y*(tK+tE))
    if polar:
        best=max(best,float((GC*(3/8-x+th*(0.25-y))).max())+tK*(th*(0.25-y)-1/8)+tE*(th*(0.25-y)+3/8))
    cache[key]=best; return best
def ok_dir(x,y,tC,tD,tab):
    return bool((blocks_val(x,y,tC,tD,0.0,sigma_vec(tC,tD,tab))<-1e-9).all())
RH=np.linspace(0,1.5,61)
THS=np.linspace(0.1,1.0,10)
def okd_theta(x,y,th,tab,polar):
    if polar and x<=3/8+th*(0.25-y)+1e-9: return False
    dl=np.linspace(0,th,41)
    if not ok_dir(x,y,np.ones_like(dl),dl,tab): return False
    dl2=np.linspace(0,1-th,41)
    if 1-th>1e-9 and not ok_dir(x+y,-y,np.ones_like(dl2),dl2,tab): return False
    return True
def okd(x,y,tab,polar=False): return any(okd_theta(x,y,th,tab,polar) for th in THS)
def best_theta(x,y,tab,polar): return [th for th in THS if okd_theta(x,y,th,tab,polar)]
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
POLAR={k:(int((M@np.array([1,0]))[1])==0) for k,(M,c) in ch.items()}
def ld_polar(y,tK,tE):   # polar exponent of an R_d point (x=1/2, y), theta=1/2
    return 0.5*(0.25-y)*(tK+tE)-tK/8+3/8*tE
def local_points(lam,tab):
    L0=13/8
    def Lam(x):
        if x>0: return max(1.0,L0-3*x)
        if x==0: return lam
        return max(1-6*x,L0-3*x,lam-6*x)
    g=lambda x: Lam(x)/2 if x<=0.5 else 0.5
    e=1e-6; loc={False:[],True:[]}
    for pol in (False,True):
        L=loc[pol]
        for yp in np.linspace(-3,0,61): L.append(((0.5-yp+e,yp),'ldfe',-yp,None))
        for yp in np.linspace(0,3,31):
            L.append(((0.5+e if not pol else max(0.5+e,3/8+0.5*(0.25-yp)+e),yp),'ld',0.0,None))
        for xp in np.linspace(-1.5,1.5,61): L.append(((xp,g(xp)+e),'zero',0.0,None))
        for yp in np.linspace(0,1.6,17):
            xb=bis(lambda xx: okd(xx,yp,tab,pol),-2,2)
            if xb is None: continue
            for dx in (1e-4,1e-3,0.05):
                ths=best_theta(xb+dx,yp,tab,pol)
                if ths: L.append(((xb+dx,yp),'blk',0.0,ths))
        for xp in np.linspace(0,0.5,11):
            yb=bis(lambda yy: okc(xp,yy,tab),-2,4)
            if yb is not None: L.append(((xp,yb+1e-3),'zero',0.0,None))
    return loc
def local_b(item,tK,tE,tab,pol,cache):
    (vp,kind,ex,ths)=item
    if kind=='zero': return 0.0
    if kind=='ld':
        b=tE/6
        if pol and vp[1]<0.25: b=max(b,ld_polar(vp[1],tK,tE))
        return b
    if kind=='ldfe':
        b=tE/6
        if pol and ex<0.25: b=max(b,ld_polar(ex,tK,tE))      # image point (1/2, -y) needs the polar term
        return b+ex*(tK+tE)
    return min(b_twist(vp[0],vp[1],tK,tE,tab,cache,th=th,polar=pol) for th in ths)
def points(loc,tK,tE,tab,cache):
    pts=[]
    lb={pol:[local_b(it,tK,tE,tab,pol,cache) for it in loc[pol]] for pol in (False,True)}
    for k,(M,c) in ch.items():
        pol=POLAR[k]; Mi=np.linalg.inv(M)
        for it,b in zip(loc[pol],lb[pol]):
            v=Mi@np.array(it[0]); pts.append((v[0],v[1],b+(c@v)*(tK+tE)))
    return np.array(pts)
def lpmin(P,obj,ylo,yhi):
    n=len(P); A_ub=np.vstack([P[:,1],-P[:,1],-(3*P[:,0]+P[:,1]),-(2*P[:,0]+P[:,1]),-(P[:,0]+P[:,1])])
    r=linprog(obj,A_ub=A_ub,b_ub=[yhi,-ylo,0,0,0],A_eq=np.ones((1,n)),b_eq=[1],bounds=(0,None),method="highs")
    return r.fun if r.status==0 else 9e9
def Xq(loc,q,eta,eps,tab,cache):
    P=points(loc,eta,2*eps,tab,cache); o=P[:,0]+P[:,2]
    main=0.75*q+0.5+lpmin(P,o,0.25,0.25)
    dual=0.5+0.5*q+lpmin(P,o+q*P[:,1],0.0,0.2499)
    nos=0.5+0.5*q+lpmin(P,o+q*P[:,1],0.2501,5.0)
    return min(nos,max(main,dual))
def lXi(a,b): return max(a,b,2/3*(a+b))
def s52(loc,q,tab,bcache,n=13):
    H=2*q-1
    if H<=0: return 2.0
    etas=np.linspace(0,H,n)
    cache={}
    def X(g,t):
        k=(round(g,6),round(t,6))
        if k not in cache: cache[k]=Xq(loc,q,g,t,tab,bcache)
        return cache[k]
    def gen(g): return 2-q+(H-g)/4+lXi(g,q)
    def e0(g,s):
        emax=(H-g)/4
        for t in np.linspace(0,emax,41):
            if all(2-q+(H-g)/4-u+lXi(g,q-u)<=s+1e-9 for u in np.linspace(t,emax,9) if u>0): return t
        return emax
    def feasible(s):
        best=9
        for i,k0 in enumerate(etas):
            m=2.0
            for j in range(i+1):
                g=etas[j]; ee=e0(g,s)
                xe=max(X(g,t) for t in np.linspace(0,ee,3))
                m=max(m,(1+(1-q)+g+(H-g)/4+2*xe)/2)
            for g in etas[i:]: m=max(m,gen(g))
            best=min(best,m)
        return best<=s+1e-9
    lo,hi=2.0,3.0
    if not feasible(hi): return 3.0
    for _ in range(20):
        mm=(lo+hi)/2
        if feasible(mm): hi=mm
        else: lo=mm
    return hi
if __name__=="__main__":
    import time
    tab=None; lam=13/8; t0=time.time()
    for it in range(60):
        loc=local_points(lam,tab); bcache={}
        new=[]
        for q in QS:
            v=S(q,tab)
            if 0.666<=q<=1.2: v=min(v,s52(loc,q,tab,bcache))
            new.append(v)
        tab=np.array(new) if tab is None else np.minimum(tab,np.array(new))
        lam_new=1.5*(float(tab[int(np.argmin(abs(QS-2/3)))])-1)
        print("iter %d (%.0fs): lam %.4f -> (L3)<<Q^%.4f | s(2/3,0.8,1,1.2)=%s"%(it,time.time()-t0,lam,lam_new,np.round([np.interp(qq,QS,tab) for qq in (2/3,0.8,1.0,1.2)],4)),flush=True)
        json.dump({"QS":QS.tolist(),"tab":tab.tolist(),"lam":min(lam,lam_new)},open('recur11_tab.json','w'))
        if it>=1 and abs(lam_new-lam)<2e-5: break
        lam=min(lam,lam_new)
