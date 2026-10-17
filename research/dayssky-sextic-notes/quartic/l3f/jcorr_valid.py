# S52 bookkeeping with the j-coprimality correction: X(h0,j) pieces = first moments at twist h0*e^2 (e | j, N(e) <= J), c-twist kappa=e
# (kappa-cost ~0.75*eps cancels the Q/N(e) shortening ~0.75*eps on the main term), so e_X^eff(eta) = max_{eps<=(H-eta)/4} e_X(eta+2eps).
import sys
fix=sys.argv[1] if len(sys.argv)>1 else "a"
if fix=="a":
    exec(open('fixa6.py').read().split("lam=13/8")[0]); X=lambda P,e: Xexp2(P,e)[0]
else:
    import sys as _s; _s.argv=['x','all','0.1667']; exec(open('hull6.py').read().split("which=sys.argv")[0]); X=lambda P,e: Xexp(P,e)[0]
import numpy as np
H=1/3
def book(P,n=25):
    etas=np.linspace(0,2*H,2*n)          # need e_X up to eta+2eps <= (eta+H)/2+... <= H; sample to 2H for safety
    ex={round(e,6):X(P,e) for e in etas}
    def eX(eta):  # interpolate
        ks=sorted(ex); return np.interp(eta,ks,[ex[k] for k in ks])
    eff=lambda eta: max(eX(eta+2*t) for t in np.linspace(0,(H-eta)/4,9))
    grid=np.linspace(0,H,n); best=(9,None)
    for i,k0 in enumerate(grid):
        offT=max(1/3+g+(H-g)/4+2*eff(g) for g in grid[:i+1])
        s=max(2+(H-k0)/4,(1+offT)/2,2)
        if s<best[0]: best=(s,k0)
    return best
lam=13/8
for it in range(5):
    P=sample(lam,"all"); s,k0=book(P)
    print("fix=%s lam=%.4f: S<<Z^%.4f  (L3)<<Q^%.4f  k0=%.3f"%(fix,lam,s,1.5*(s-1),k0),flush=True)
    if abs(lam-1.5*(s-1))<1e-4: break
    lam=1.5*(s-1)
