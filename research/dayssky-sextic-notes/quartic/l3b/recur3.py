# Exact-ish check: apply the recursion rules to candidate functions and see what they prove.
import numpy as np
def L(r):  # dFDH lower bound (log_B), A=B^r rows
    return max(r,1.0,0.75*r+0.5,0.5*r+0.75)
def BGL(r): return max(r,1.0,(2/3)*(r+1))
def gall(g,eta,n=60):
    worst=-1e9; best_box=None
    for x2 in np.linspace(0,eta/2,n):
        for x3 in np.linspace(0,(eta-2*x2)/3,n//2):
            for y in np.linspace(0,(eta-2*x2-3*x3)/4,n//2):
                x1=max(0.0,eta-2*x2-3*x3-4*y)
                c1=x2+x3+y+g(x1); c2=x1+x3+y+max(x2,1.0); c3=x1+x2+y+g(x3); c4=y+max(2*(x1+x2+x3),1.0)
                v=min(c1,c2,c3,c4)
                if v>worst: worst=v; best_box=(x1,x2,x3,y)
    return worst,best_box
def apply(g,r):
    out={'self':g(r)}
    if r>0: out['dual']=r*g(1/r)
    if r<2:
        ga,box=gall(g,2-r); out['poisson']=max(r,r-1+ga); out['box']=box
    return out
for name,g in (("dFDH lower bound L",L),("BGL upper",BGL)):
    print("== rules applied to",name)
    for r in (0.5,2/3,0.8,1.0,1.25,4/3,1.5,1.8):
        o=apply(g,r)
        best=min(v for k,v in o.items() if k in('self','dual','poisson'))
        print("  r=%.3f  self=%.4f dual=%.4f poisson=%.4f  -> best %.4f  (L=%.4f)%s  worst dual box %s"%(r,o['self'],o.get('dual',9),o.get('poisson',9),best,L(r)," <-- BELOW proven lower bound" if best<L(r)-1e-6 else "", tuple(round(x,3) for x in o.get('box',())) ))
