# Quartic large sieve recursion with square-class bookkeeping (exponents, log_N units).
# g4(r): squarefree rows K=N^r, squarefree columns N (quartic symbol). Lower bounds: max(r,1), theta (r+1)/2.
# gall(eta): ALL rows h <= N^eta (h = h1 h2^2 h3^3 j^4), squarefree columns.
import numpy as np, itertools
R=np.linspace(0.0,4.0,401)
def ev(g,x):
    x=float(x)
    if x<=R[-1]: return float(np.interp(x,R,g))
    return x*float(np.interp(1/x,R,g)) if 1/x>=R[0] else max(x,1)
def gall(g,eta,step=0.05):
    worst=-1e9
    for x2 in np.arange(0,eta/2+1e-9,step):
        for x3 in np.arange(0,(eta-2*x2)/3+1e-9,step):
            for y in np.arange(0,(eta-2*x2-3*x3)/4+1e-9,step):
                    x1=max(0.0,eta-2*x2-3*x3-4*y)   # costs increase in x1: worst box saturates the norm constraint
                    c1=x2+x3+y+ev(g,x1)                 # sweep h1, quartic LS
                    c2=x1+x3+y+max(x2,1.0)              # sweep h2, quadratic LS (optimal)
                    c3=x1+x2+y+ev(g,x3)                 # sweep h3, conjugate quartic LS
                    c4=y+max(2*(x1+x2+x3),1.0)          # classical LS over all chars of conductor rad
                    worst=max(worst,min(c1,c2,c3,c4))
    return worst
g=np.array([min(max(r,1,0.75*(r+1)),max(2*r,1)) for r in R])  # GZ + classical
lower=np.maximum(np.maximum(R,1.0),(R+1)/2)
for it in range(12):
    old=g.copy()
    new=g.copy()
    for i,r in enumerate(R):
        cands=[g[i]]
        if r>0: cands.append(r*ev(g,1/r))               # duality (both sqfree)
        if r<2: cands.append(max(r, r-1+gall(g,2-r,step=0.025)))  # Poisson in rows -> all dual rows
        new[i]=min(cands)
    g=np.maximum(new,lower)
    viol=np.max(lower-new)
    print("iter",it,"max change",np.max(np.abs(g-old)),"lower-bound violation before clamp %.3f"%viol, " g(2/3)=%.4f g(1)=%.4f g(4/3)=%.4f"%(ev(g,2/3),ev(g,1),ev(g,4/3)))
    if np.max(np.abs(g-old))<1e-6: break
np.save('g4_fixed.npy',g)
for r in (0.5,2/3,0.8,1.0,1.25,4/3,1.5,2.0):
    print("r=%.3f  g4=%.4f  GZ=%.4f  lower=%.4f"%(r,ev(g,r),max(r,1,0.75*(r+1)),max(r,1,(r+1)/2)))
