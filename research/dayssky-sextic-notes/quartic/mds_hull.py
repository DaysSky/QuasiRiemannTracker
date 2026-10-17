# Second-moment multiple Dirichlet series  Z(s1,s3,w) = sum_m N(m)^{-w} D_m(s1) conj D_m(s3)
# centered coords x1=s1-1/2, x3=s3-1/2, y=w-1/2.  Reflections:
#   tau1: x1->-x1, y->y+k*x1   (row FE, conductor N(m)^k ; k=3 quartic, k=2 cubic, k=1 quadratic analogue)
#   tau3: same in x3;   tauw: y->-y, x1->x1+y, x3->x3+y  (Hecke L(w, conj chi_c1 chi_c3), conductor N(c1 c3))
# Region of absolute convergence R = {x1>1/2, x3>1/2, y>1/2}.
# Compute min y such that (0,0,y) lies in conv hull of union of g(R), g in words of length <= L (Bochner).
import numpy as np, itertools, sys
from scipy.optimize import linprog
def gens(k):
    t1=np.array([[-1,0,0],[0,1,0],[k,0,1]],float)
    t3=np.array([[1,0,0],[0,-1,0],[0,k,1]],float)
    tw=np.array([[1,0,1],[0,1,1],[0,0,-1]],float)
    return [t1,t3,tw]
def orbit(k,L):
    G=gens(k); els={tuple(np.eye(3).flatten())}; frontier=[np.eye(3)]
    for _ in range(L):
        new=[]
        for g in frontier:
            for s in G:
                h=s@g; key=tuple(np.round(h,9).flatten())
                if key not in els: els.add(key); new.append(h)
        frontier=new
    return [np.array(e).reshape(3,3) for e in els]
def miny(k,L,eps=0.5):
    gs=orbit(k,L); v0=np.array([eps,eps,eps]); E=np.eye(3)
    pts=[g@v0 for g in gs]; dirs=[g@E[:,i] for g in gs for i in range(3)]
    n,m=len(pts),len(dirs)
    # variables: lambda (n, sum=1, >=0), mu (m, >=0); constraint: sum lam p + sum mu d = (0,0,y); minimize y
    c=np.concatenate([np.array([p[2] for p in pts]),np.array([d[2] for d in dirs])])
    Aeq=np.zeros((3,n+m)); beq=np.zeros(3)
    for j,p in enumerate(pts): Aeq[0,j]=p[0]; Aeq[1,j]=p[1]; Aeq[2,j]=1
    for j,d in enumerate(dirs): Aeq[0,n+j]=d[0]; Aeq[1,n+j]=d[1]
    beq[2]=1
    r=linprog(c,A_eq=Aeq,b_eq=beq,bounds=[(0,None)]*(n+m),method='highs')
    return len(gs), (r.fun if r.status==0 else r.status)
for k in (1,2,3):
    for L in (2,4,6,8,10,12):
        ng,y=miny(k,L)
        print(f"k={k} L={L:2d} |G_L|={ng:5d}  min y at center = {y}  (w = y+1/2)")
