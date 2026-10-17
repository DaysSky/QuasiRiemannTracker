# Fixed point of large-sieve exponent rules for the quartic family (general |a_c|<=1 coefficients).
# g(r) = log_D Delta(D^r, D): rows H = D^r, columns D. Rules (all in log units, homogeneous):
#  GZ:        g(r) <= max(r, 1, 3(r+1)/4)                  [Gao-Zhao, squarefree]
#  classical: g(r) <= max(2r, 1)
#  symmetry:  g(r) = r * g(1/r)                            [rows <-> columns via quartic reciprocity]
#  Poisson:   g(r) <= max(r, r - 1 + g(2 - r))   (r<=2)    [Poisson in rows, Gauss sums of modulus c c']
#  Hoelder:   g(r) <= r/2 + g(r/2)                         [square the column sum: length D^2]
import numpy as np
R=np.linspace(0.02,4.0,1991)
def interp(g,x):
    x=np.asarray(x,float)
    # use homogeneity outside range: g(r)=r*g(1/r)
    out=np.empty_like(x)
    for i,v in np.ndenumerate(x):
        if v< R[0]: out[i]=max(1.0,v)  # tiny row count: Delta = D
        elif v>R[-1]: out[i]=v*np.interp(1/v,R,g)
        else: out[i]=np.interp(v,R,g)
    return out
g=np.array([min(max(r,1,0.75*(r+1)),max(2*r,1)) for r in R])
for it in range(200):
    old=g.copy()
    sym=R*interp(g,1/R)
    poi=np.where(R<=2, np.maximum(R, R-1+interp(g,2-R)), np.inf)
    hol=R/2+interp(g,R/2)
    g=np.minimum.reduce([g,sym,poi,hol])
    g=np.maximum(g,np.maximum(R,1.0))
    if np.max(np.abs(g-old))<1e-12: break
print("iterations",it)
for r in (0.5,2/3,0.75,1.0,4/3,1.5,2.0):
    v=np.interp(r,R,g); print("r=%.4f  g=%.4f  excess over max(r,1) = %.4f"%(r,v,v-max(r,1)))
# (L3) exponent: S/Z with Mobius side at H=Z^{4/3}: lambda = g(4/3)-4/3 in units of D=Z; (L3) sum ~ Q^{3/2 + 2*lambda*... }
lam=np.interp(4/3,R,g)-4/3
print("Mobius-side excess at H=D^{4/3}: lambda=%.4f -> (L3) sum << Q^{%.4f}"%(lam, 1.5+1.5*(4/3)*lam*3/4*2/1.5 if False else 1.5+2*lam*0.75))
