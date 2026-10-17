# Full W(G2) mean-value Bochner hull for the first moment, then S50 bookkeeping with bootstrapping.
# Chamber w: v' = M_w v (v=(x,y)); twist t_w = M_w (1,0); h0-factor exponent phi_w(v) = c_w . v (affine part from d-FEs).
# Regions in chamber coordinates (x',y'):
#   Rd: x' > 1/2 + max(0,-y')    b = 1/6 + max(0,-y')     (LoA for the d-node Kubota family, index h0*c)
#   Rc: y' > g(x')               b = 0                     ((L3)-type mean square, exponent lam)
import numpy as np, sys
BD=float(sys.argv[2]) if len(sys.argv)>2 else 0.0
from scipy.optimize import linprog
sc=np.array([[-1,0],[3,1]]); sd=np.array([[1,1],[0,-1]])
# enumerate chambers with phi: Z(v) = N^{phi_w(v)} Z_w(M_w v); applying s_d at chamber w: Z_w(v') = N^{-y'} Z_{dw}(sd v')
ch={(1,0,0,1):(np.eye(2,dtype=int),np.zeros(2),"")}; fr=[(1,0,0,1)]
while fr:
    nf=[]
    for k in fr:
        M,c,word=ch[k]
        for nm,g in (("c",sc),("d",sd)):
            N=g@M; kk=tuple(N.flatten())
            cn=c.copy()
            if nm=="d": cn=c-M[1]          # phi gains -(M v)_y
            if kk not in ch: ch[kk]=(N,cn,nm+word); nf.append(kk)
    fr=nf
assert len(ch)==12
def sample(lam,which):
    g=lambda x: lam/2-3*x if x<=0 else (lam/2 if x<=0.5 else 0.5)
    h=lambda y: 0.5+max(0,-y)
    e=1e-7; pts=[]
    for k,(M,c,word) in ch.items():
        if word not in which and which!="all": continue
        Mi=np.linalg.inv(M)
        loc=[]
        for yp in np.linspace(-3,3,301): loc.append(((h(yp)+e,yp),BD+max(0,-yp)))
        for xp in np.linspace(-1.5,1.5,301): loc.append(((xp,g(xp)+e),0.0))
        for (vp,b) in loc:
            v=Mi@np.array(vp); pts.append((v[0],v[1],b+c@v))
    return np.array(pts)
def lpmin(P,obj,ylo,yhi,xlo=0.0):
    # min sum t_i obj_i s.t. sum t=1, ylo<=sum t y<=yhi, sum t x>=xlo, t>=0
    n=len(P)
    A_eq=np.ones((1,n)); b_eq=[1]
    A_ub=np.vstack([P[:,1],-P[:,1],-P[:,0]]); b_ub=[yhi,-ylo,-xlo]
    r=linprog(obj,A_ub=A_ub,b_ub=b_ub,A_eq=A_eq,b_eq=b_eq,bounds=(0,None),method="highs")
    return r.fun if r.status==0 else 9e9
def argmin_combo(P,obj,ylo,yhi,xlo=0.0):
    n=len(P); A_eq=np.ones((1,n)); A_ub=np.vstack([P[:,1],-P[:,1],-P[:,0]])
    r=linprog(obj,A_ub=A_ub,b_ub=[yhi,-ylo,-xlo],A_eq=A_eq,b_eq=[1],bounds=(0,None),method="highs")
    idx=np.where(r.x>1e-9)[0]; return [(round(r.x[i],4),np.round(P[i],4)) for i in idx]
def Xexp(P,eta):
    main=1.0+lpmin(P,P[:,0]+eta*P[:,2],0.25,0.25)
    dual=5/6+lpmin(P,P[:,0]+(2/3)*P[:,1]+eta*P[:,2],0.0,0.2499)
    nos=5/6+lpmin(P,P[:,0]+(2/3)*P[:,1]+eta*P[:,2],0.2501,5.0)
    return min(nos,max(main,dual)),main,dual,nos
def book(P,H=1/3,n=25):
    etas=np.linspace(0,H,n); ex=[Xexp(P,e)[0] for e in etas]
    best=(9,None)
    for i,k0 in enumerate(etas):
        offT=max(1/3+etas[j]+(H-etas[j])/4+2*ex[j] for j in range(i+1))
        s=max(2+(H-k0)/4,(1+offT)/2,2)
        if s<best[0]: best=(s,k0)
    return best
which=sys.argv[1] if len(sys.argv)>1 else "all"
words=[w for (_,_,w) in ch.values()]
print("chambers:",sorted(words,key=len))
lam=13/8
for it in range(8):
    P=sample(lam,which if which!="all" else "all")
    s,k0=book(P)
    x0=Xexp(P,0.0)
    print("iter %d lam=%.4f  X(1): main %.4f dual %.4f noshift %.4f | S<<Z^%.4f  (L3)<<Q^%.4f  k0=%.3f"%(it,lam,x0[1],x0[2],x0[3],s,1.5*(s-1),k0),flush=True)
    new=1.5*(s-1)
    if lam-new<1e-4: break
    lam=new
