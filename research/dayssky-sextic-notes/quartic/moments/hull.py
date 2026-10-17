# Bochner-hull bookkeeping for the second-moment triple series
#   Z(s1,s2,w) = sum_m psi(m^3,s1) conj psi(m^3,conj s2) N(m)^{-w}
# coordinates v=(x1,x2,y): x_i = Re s_i - 1/2, y = Re w - 1/2.  Target: (L^3) <=> holomorphy at (0,0,1/2+).
# tau_i : Kubota FE of psi(m^3,s_i) (conductor N(m)^3);  tau_w : FE of L(w, conj chi_c1 chi_c2) (conductor N(c1c2)).
# tau_w maps Z to Z' = sum_m N(m)^{-w'} M_{u1}(m) conj M_{u2}(m), M_u(m) = sum_c mu(c) chi_c(m) N(c)^{-u}  (~1/L), because
# a_c gt(conj chi_c) = mu(c)  (verified in ident.py).  Z' has NO usable FE in u_i (1/L), so only <tau1,tau2> acts.
import numpy as np, itertools
from scipy.optimize import linprog
T1=np.array([[-1,0,0],[0,1,0],[3,0,1]],float)
T2=np.array([[1,0,0],[0,-1,0],[0,3,1]],float)
TW=np.array([[1,0,1],[0,1,1],[0,0,-1]],float)
E=np.eye(3)
def poly(verts,rays): return ([np.array(v,float) for v in verts],[np.array(r,float) for r in rays])
def image(M,P): return ([M@v for v in P[0]],[M@r for r in P[1]])
# --- base regions for Z (unconditional) ---
RA=poly([(.5,.5,.5)],[E[0],E[1],E[2]])                                   # absolute convergence
Rm=poly([(.5,.5,.5),(-.5,.5,2),(.5,-.5,2),(-.5,-.5,3.5)],[E[0],E[1],E[2],(-1,0,3),(0,-1,3)])  # m-sum, psi convexity
Rc=poly([(.5,.5,.5),(1,1,-.5)],[E[0],E[1],E[2],(1,1,-1)])               # (c1,c2)-sum, L convexity
# --- base regions for Z' (primed coords), transported by tau_w ---
def Rp_m(sig0): return poly([(sig0-.5,sig0-.5,.5)],[E[0],E[1],E[2]])     # m-sum of 1/L * 1/L, zero-free for Re u > sig0
Rp_c=Rc
def Rp_hyp(d0,th):  # hypothetical: sum_{N k<=H} M_{u1} conj M_{u2} = main + O(H^th) for Re u_i >= 1/2+d0
    return poly([(d0,d0,th-.5)],[E[0],E[1],E[2]])
def orbit_group(gens,L):
    els={tuple(np.eye(3).ravel())}; fr=[np.eye(3)]
    for _ in range(L):
        new=[]
        for g in fr:
            for s in gens:
                h=s@g; k=tuple(np.round(h,9).ravel())
                if k not in els: els.add(k); new.append(h)
        fr=new
    return [np.array(k).reshape(3,3) for k in els]
def center_min(regionsZ,regionsZp,group):
    V=[];R=[]
    base=list(regionsZ)+[image(TW,P) for P in regionsZp]
    for g in group:
        for P in base:
            Q=image(g,P); V+=Q[0]; R+=Q[1]
    n,m=len(V),len(R)
    c=np.array([v[2] for v in V]+[r[2] for r in R])
    A=np.zeros((3,n+m)); b=np.array([0,0,1.])
    for j,v in enumerate(V): A[0,j]=v[0];A[1,j]=v[1];A[2,j]=1
    for j,r in enumerate(R): A[0,n+j]=r[0];A[1,n+j]=r[1]
    res=linprog(c,A_eq=A,b_eq=b,bounds=[(0,None)]*(n+m),method='highs')
    return res.fun if res.status==0 else ('unbounded' if res.status==3 else res.status)
lam=lambda v: v[2]+1.5*(v[0]+v[1])
G4=orbit_group([T1,T2],4)
print("|<tau1,tau2>| =",len(G4))
y=center_min([RA,Rm,Rc],[Rp_m(1.0),Rp_c],G4)
print(f"unconditional, valid FEs only:   min y at center = {y:.4f}  -> theta = y+1/2 = {y+.5:.4f}  (convexity = 5/2)")
for s0 in (0.5,0.6,2/3,0.75,5/6,0.9,1.0):
    y=center_min([RA,Rm,Rc],[Rp_m(s0),Rp_c],G4)
    print(f"zero-free Re u>{s0:.4f}:  min y = {y:.4f}  theta = {y+.5:.4f}   predicted min(5/2,3*sig0) = {min(2.5,3*s0):.4f}")
for d0,th in ((0.0,0.75),(0.0,0.5),(0.05,0.70),(0.1,0.6),(0.0,1.0),(0.0,0.875)):
    y=center_min([RA,Rm,Rc],[Rp_m(0.5),Rp_c,Rp_hyp(d0,th)],G4)
    print(f"GRH + hyp. ratio asymptotic error H^{th:.3f} at Re u=1/2+{d0:.2f}: min y={y:.4f} theta={y+.5:.4f}  pred 3d0+2th-1/2={min(1.5,3*d0+2*th-.5):.4f}")
# hypothetical full Coxeter group (if the 1/L-series had its FEs): still convexity, by the imaginary root lambda
for L in (4,8,12):
    G=orbit_group([T1,T2,TW],L)
    y=center_min([RA,Rm,Rc],[],G)
    print(f"full Coxeter group (words<= {L}, |G|={len(G)}): min y = {y:.4f}")
# Kac-lemma check: lambda o g - lambda has nonnegative coefficients in the basis (x1,x2,y) for all g in the Coxeter group
G=orbit_group([T1,T2,TW],14); lv=np.array([1.5,1.5,1.0]); worst=1e9
for g in G:
    d=lv@g-lv; worst=min(worst,d.min())
print("words<=14:",len(G),"elements; min coefficient of (lambda o g - lambda) =",worst)
