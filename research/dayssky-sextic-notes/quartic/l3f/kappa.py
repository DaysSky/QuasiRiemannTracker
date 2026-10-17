# kappa-exponent of the X(h0) bound for the kappa-twisted c-node family (S52 §7).
import numpy as np, sys
from scipy.optimize import linprog
sc=np.array([[-1,0],[3,1]]); sd=np.array([[1,1],[0,-1]])
ch={(1,0,0,1):(np.eye(2,dtype=int),np.zeros(2),np.zeros(2),"")}; fr=[(1,0,0,1)]
while fr:
    nf=[]
    for k in fr:
        M,c,ck,word=ch[k]
        for nm,g in (("c",sc),("d",sd)):
            N=g@M; kk=tuple(N.flatten()); cn=c.copy(); ckn=ck.copy()
            if nm=="d": cn=c-M[1]
            else: ckn=ck-M[0]          # c-FE factor N(kappa)^{-(M v)_x}
            if kk not in ch: ch[kk]=(N,cn,ckn,nm+word); nf.append(kk)
    fr=nf
lam=float(sys.argv[1]) if len(sys.argv)>1 else 1.5417
th=float(sys.argv[2]) if len(sys.argv)>2 else 0.5     # mean-square kappa exponent of the c-family at the centre
L0=13/8
def Lam(x):
    if x>0: return max(1.0,L0-3*x)
    if x==0: return lam
    return max(1-6*x,L0-3*x,lam-6*x)
def thm(x): return th if x>=0 else th-2*x            # mean-square kappa-exponent (FE adds N(kappa)^{-2x})
g=lambda x: Lam(x)/2 if x<=0.5 else 0.5
h=lambda y: (0.5-y) if y<0 else max(0.0,0.5-0.4*y)
e=1e-7; pts=[]
for k,(M,c,ck,word) in ch.items():
    Mi=np.linalg.inv(M)
    loc=[]
    for yp in np.linspace(-3,3,601):
        # R_d / R_d': kappa enters only through the convexity block (3/4 per unit of the c-series index); take 3/4 on the R_d' branch with y'>0
        loc.append(((h(yp)+e,yp),max(0,-yp),0.75 if yp>0 else 0.0))
    for xp in np.linspace(-1.5,1.5,301):
        loc.append(((xp,g(xp)+e),0.0,thm(xp)/2 if xp<=0.5 else 0.0))
    for (vp,b,bk) in loc:
        v=Mi@np.array(vp); pts.append((v[0],v[1],b+c@v,bk+ck@v))
P=np.array(pts)
def opt(obj,ylo,yhi):
    n=len(P); A_ub=np.vstack([P[:,1],-P[:,1],-P[:,0]])
    r=linprog(obj,A_ub=A_ub,b_ub=[yhi,-ylo,0],A_eq=np.ones((1,n)),b_eq=[1],bounds=(0,None),method="highs")
    return r
for eta in (0.0,0.1,0.2):
    for wk in (0.0,0.05,0.2):     # penalise kappa slightly to find the kappa-cheapest optimum
        r=opt(P[:,0]+eta*P[:,2]+wk*P[:,3],0.25,0.25)
        x=r.x; xv=x@P[:,0]; bv=x@P[:,2]; kv=x@P[:,3]
        print("eta=%.2f wk=%.2f: main x=%.4f b=%.4f kappa-exp=%.4f  (X exponent %.4f)"%(eta,wk,xv,bv,kv,1+xv+eta*bv))
