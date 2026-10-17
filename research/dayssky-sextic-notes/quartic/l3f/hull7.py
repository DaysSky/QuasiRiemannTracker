# Full W(G2) hull with block-optimal regions (fix (a)), j-corrected S50 bookkeeping, bootstrap.
import numpy as np, sys
from scipy.optimize import linprog
from blocks import boundary_d, boundary_c, b_d
sc=np.array([[-1,0],[3,1]]); sd=np.array([[1,1],[0,-1]])
ch={(1,0,0,1):(np.eye(2,dtype=int),np.zeros(2),"")}; fr=[(1,0,0,1)]
while fr:
    nf=[]
    for k in fr:
        M,c,word=ch[k]
        for nm,g in (("c",sc),("d",sd)):
            N=g@M; kk=tuple(N.flatten()); cn=c.copy()
            if nm=="d": cn=c-M[1]
            if kk not in ch: ch[kk]=(N,cn,nm+word); nf.append(kk)
    fr=nf
def local_points(lam):
    L0=13/8
    def Lam(x):
        if x>0: return max(1.0,L0-3*x)
        if x==0: return lam
        return max(1-6*x,L0-3*x,lam-6*x)
    g=lambda x: Lam(x)/2 if x<=0.5 else 0.5
    e=1e-6; loc=[]
    for yp in np.linspace(-3,0,121): loc.append(((0.5-yp+e,yp),1/6-yp))           # old R_d, y'<0 (FE)
    for yp in np.linspace(0,3,121): loc.append(((max(0.0,0.5-0.4*yp)+e,yp),1/6))    # old R_d'
    for xp in np.linspace(-1.5,1.5,121): loc.append(((xp,g(xp)+e),0.0))            # old R_c
    for yp in np.linspace(0,1.6,33):                                                  # block Region_d
        xb=boundary_d(yp,lam)
        if xb is not None: loc.append(((xb+1e-3,yp),b_d(xb+1e-3,yp,lam,n=30)))
    for xp in np.linspace(0,0.5,21):                                                  # block Region_c
        yb=boundary_c(xp,lam)
        if yb is not None: loc.append(((xp,yb+1e-3),0.0))
    return loc
def sample(lam):
    loc=local_points(lam); pts=[]
    for k,(M,c,word) in ch.items():
        Mi=np.linalg.inv(M)
        for (vp,b) in loc:
            v=Mi@np.array(vp); pts.append((v[0],v[1],b+c@v))
    return np.array(pts)
def lpmin2(P,obj,ylo,yhi):
    n=len(P); A_ub=np.vstack([P[:,1],-P[:,1],-(3*P[:,0]+P[:,1]),-(2*P[:,0]+P[:,1]),-(P[:,0]+P[:,1])])
    r=linprog(obj,A_ub=A_ub,b_ub=[yhi,-ylo,0,0,0],A_eq=np.ones((1,n)),b_eq=[1],bounds=(0,None),method="highs")
    return r.fun if r.status==0 else 9e9
def Xexp(P,eta):
    main=1.0+lpmin2(P,P[:,0]+eta*P[:,2],0.25,0.25)
    dual=5/6+lpmin2(P,P[:,0]+(2/3)*P[:,1]+eta*P[:,2],0.0,0.2499)
    nos=5/6+lpmin2(P,P[:,0]+(2/3)*P[:,1]+eta*P[:,2],0.2501,5.0)
    return min(nos,max(main,dual))
H=1/3
def book(P,n=25,jcorr=True):
    etas=np.linspace(0,H,n); ex=[Xexp(P,e) for e in etas]
    eX=lambda t: np.interp(t,etas,ex)
    eff=(lambda eta: max(eX(eta+2*t) for t in np.linspace(0,(H-eta)/4,9))) if jcorr else eX
    best=(9,None)
    for i,k0 in enumerate(etas):
        offT=max(1/3+g+(H-g)/4+2*eff(g) for g in etas[:i+1])
        s=max(2+(H-k0)/4,(1+offT)/2,2)
        if s<best[0]: best=(s,k0)
    return best,ex
if __name__=="__main__":
    lam=1.5573
    for it in range(6):
        P=sample(lam); (s,k0),ex=book(P)
        print("lam=%.4f: X(eta=0)=%.4f X(1/6)=%.4f X(1/3)=%.4f | S<<Z^%.4f (L3)<<Q^%.4f k0=%.3f"%(lam,ex[0],ex[12],ex[-1],s,1.5*(s-1),k0),flush=True)
        if abs(lam-1.5*(s-1))<1e-4: break
        lam=min(lam,1.5*(s-1))
