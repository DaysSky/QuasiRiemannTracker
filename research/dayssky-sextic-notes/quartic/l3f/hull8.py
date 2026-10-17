# Two twist exponents: b_rms (RMS over h0 ~ K, vector-valued Bochner) and b_pt (pointwise, for the e^2 part of the j-pieces).
import numpy as np
from scipy.optimize import linprog
from blocks import boundary_d, boundary_c, b_d
from blocks_rms import b_rms as brms
exec(open('hull7.py').read().split("def local_points")[0])
def local_points(lam):
    L0=13/8
    def Lam(x):
        if x>0: return max(1.0,L0-3*x)
        if x==0: return lam
        return max(1-6*x,L0-3*x,lam-6*x)
    g=lambda x: Lam(x)/2 if x<=0.5 else 0.5
    e=1e-6; loc=[]
    for yp in np.linspace(-3,0,121): loc.append(((0.5-yp+e,yp),-yp,1/6-yp))          # old R_d (FE), (b_rms, b_pt)
    for yp in np.linspace(0,3,121): loc.append(((max(0.0,0.5-0.4*yp)+e,yp),0.0,1/6))   # old R_d'
    for xp in np.linspace(-1.5,1.5,121): loc.append(((xp,g(xp)+e),0.0,0.0))           # old R_c
    for yp in np.linspace(0,1.6,33):
        xb=boundary_d(yp,lam)
        if xb is None: continue
        for dx in (1e-3,0.05,0.1):
            loc.append(((xb+dx,yp),brms(xb+dx,yp,lam,n=40),b_d(xb+dx,yp,lam,n=30)))
    for xp in np.linspace(0,0.5,21):
        yb=boundary_c(xp,lam)
        if yb is not None: loc.append(((xp,yb+1e-3),0.0,0.0))
    return loc
def sample(lam):
    loc=local_points(lam); pts=[]
    for k,(M,c,word) in ch.items():
        Mi=np.linalg.inv(M)
        for (vp,br,bp) in loc:
            v=Mi@np.array(vp); f=c@v; pts.append((v[0],v[1],br+f,bp+f))
    return np.array(pts)
def lpmin(P,obj,ylo,yhi):
    n=len(P); A_ub=np.vstack([P[:,1],-P[:,1],-(3*P[:,0]+P[:,1]),-(2*P[:,0]+P[:,1]),-(P[:,0]+P[:,1])])
    r=linprog(obj,A_ub=A_ub,b_ub=[yhi,-ylo,0,0,0],A_eq=np.ones((1,n)),b_eq=[1],bounds=(0,None),method="highs")
    return r.fun if r.status==0 else 9e9
def Xexp(P,eta,eps):   # RMS twist Z^eta, pointwise twist Z^{2 eps}
    o=P[:,0]+eta*P[:,2]+2*eps*P[:,3]
    main=1.0+lpmin(P,o,0.25,0.25)
    dual=5/6+lpmin(P,o+(2/3)*P[:,1],0.0,0.2499)
    nos=5/6+lpmin(P,o+(2/3)*P[:,1],0.2501,5.0)
    return min(nos,max(main,dual))
H=1/3
def book(P,n=17):
    etas=np.linspace(0,H,n)
    eff=[max(Xexp(P,g,t) for t in np.linspace(0,(H-g)/4,5)) for g in etas]
    best=(9,None)
    for i,k0 in enumerate(etas):
        offT=max(1/3+etas[j]+(H-etas[j])/4+2*eff[j] for j in range(i+1))
        s=max(2+(H-k0)/4,(1+offT)/2,2)
        if s<best[0]: best=(s,k0)
    return best,eff
lam=1.5573
for it in range(6):
    P=sample(lam); (s,k0),eff=book(P)
    print("lam=%.4f eff X: %s | S<<Z^%.4f (L3)<<Q^%.4f k0=%.3f"%(lam,np.round(eff[::4],4),s,1.5*(s-1),k0),flush=True)
    if abs(lam-1.5*(s-1))<1e-4: break
    lam=min(lam,1.5*(s-1))
