# Mean-value Bochner hull + S50 bookkeeping, with bootstrapping in lam.
# F = L(4s-1,a^4) zeta(4w-1) (w-1/2)(w-3/4) Z_full is holomorphic on each region below (normalisers make the FEs clean).
# Evaluation of Z_full = F / normalisers needs x >= 0 (1/L(1+4x,a^4)) and y >= 0 (1/zeta(1+4y)).
import numpy as np, sys
def regions(lam, use):
    g=lambda x: lam/2-3*x if x<=0 else (lam/2 if x<=0.5 else 0.5)
    h=lambda y: 0.5+max(0,-y)
    e=1e-7; P=[]
    for y in np.linspace(-4,4,1601): P.append((h(y)+e,y,1/6+max(0,-y)))           # R_d
    for x in np.linspace(-2,2,1601): P.append((x,g(x)+e,0.0))                       # R_c
    if 'D' in use:
        for u in np.linspace(-2,2,1601):
            y=-(g(u)+e); P.append((u-y,y,-y))                                         # Sd^{-1} R_c, b = -y
    if 'C' in use:
        for v in np.linspace(-4,4,1601):
            x=-(h(v)+e); P.append((x,v-3*x,1/6+max(0,-v)))                            # Sc^{-1} R_d
    return np.array(P)
def best_on_line(P,y0,eta,xmin=0.0):
    A=P[P[:,1]<=y0]; B=P[P[:,1]>=y0]; best=9e9
    for a in A[::2]:
        t=(B[:,1]-y0)/np.maximum(B[:,1]-a[1],1e-300)
        x=t*a[0]+(1-t)*B[:,0]; b=t*a[2]+(1-t)*B[:,2]
        ok=x>=xmin
        if ok.any():
            v=(x+eta*b)[ok].min()
            best=min(best,v)
    return best
def X_exponent(P,eta):
    main=1.0+best_on_line(P,0.25,eta)                       # Q^{3/4} Z^{1/2+x} h0^b with Q=Z^{2/3}
    dual=min(0.5+(2/3)*(0.5+y0)+best_on_line(P,y0,eta) for y0 in np.linspace(0,0.24,13))
    noshift=min(0.5+(2/3)*(0.5+y0)+best_on_line(P,y0,eta) for y0 in np.linspace(0.26,1.0,12))
    return min(noshift,max(main,dual)),main,dual
def bookkeeping(P,H=1/3):
    etas=np.linspace(0,H,25); ex=[X_exponent(P,e)[0] for e in etas]
    res=(9,None)
    for i,k0 in enumerate(etas):
        offT=max(1/3+etas[j]+(H-etas[j])/4+2*ex[j] for j in range(i+1))
        s=max(2+(H-k0)/4,(1+offT)/2,2)
        if s<res[0]: res=(s,k0)
    return res,ex
use=sys.argv[1] if len(sys.argv)>1 else "dcDC"
lam=13/8
for it in range(6):
    P=regions(lam,use)
    (s,k0),ex=bookkeeping(P)
    m0=X_exponent(P,0.0)
    print("iter %d lam=%.4f: X(h0=1) exp=%.4f (main %.4f, dual %.4f); S<<Z^%.4f -> (L3)<<Q^%.4f (k0=%.3f)"%(it,lam,m0[0],m0[1],m0[2],s,1.5*(s-1),k0))
    new=1.5*(s-1)
    if abs(new-lam)<1e-4: break
    lam=min(lam,new)
