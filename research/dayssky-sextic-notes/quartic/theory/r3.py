# Exponent bookkeeping for the Poisson/duality recursion (all exponents = log of size in a common base).
# GZ(x,y): Gao-Zhao squarefree x squarefree; Bsf: best squarefree bound; Ball(r,c): all rows<=L^r, sqfree cols<=L^c;
# T(q,c): all rows<=L^q (dual rows), sqfree cols c, decomposed h=h1 h2^2 h3^3 h4^4.
import sys, functools, itertools
import numpy as np
TH=0.1;KA=0.25
STEP=float(sys.argv[1]) if len(sys.argv)>1 else 1/48
DEPTH=int(sys.argv[2]) if len(sys.argv)>2 else 3
def rd(x): return round(x/STEP)*STEP
def GZ(x,y): return max(x,y,0.75*(x+y))
@functools.lru_cache(None)
def Bsf(x,y,d):
    x,y=rd(x),rd(y)
    b=GZ(x,y)
    if d<=0 or min(x,y)<=0: return b if min(x,y)>0 else max(x,y)
    r,c=max(x,y),min(x,y)
    return min(b,Ball(r,c,d-1))
@functools.lru_cache(None)
def Ball(r,c,d,mu=False):
    q=2*c-r
    if q<0: return max(r,c)
    return min(max(r,c, r-c+T(q,c,d,mu)), max(r,c,0.75*(r+c), r/4+c) )
@functools.lru_cache(None)
def T(q,c,d,mu=False):
    q=rd(q)
    best=-1e9
    n=int(round(q/STEP))
    for i4 in range(0,n//4+1):
        for i3 in range(0,(n-4*i4)//3+1):
            for i2 in range(0,(n-4*i4-3*i3)//2+1):
                i1=n-4*i4-3*i3-2*i2   # take x1 maximal (bounds monotone in x1)
                x1,x2,x3,x4=i1*STEP,i2*STEP,i3*STEP,i4*STEP
                opts=[x2+x3+Bsf(x1,c,d), x1+x3+max(x2,c), x1+x2+Bsf(x3,c,d), max(2*(x1+x2+x3),c)]
                if mu: opts.append((x1+x2+x3)+(1-2*TH)*c+2*KA*(x1+2*x2+3*x3))   # conditional (K): |Y(h0)|<<D^{1/2}N(h0)^{3/4}
                val=x4+min(opts)
                if val>best: best=val; arg=(x1,x2,x3,x4)
    return best
for r in (4/3,):
    out=[]
    for d in range(0,DEPTH+1):
        b=Ball(r,1.0,d); bm=Ball(r,1.0,d,True)
        out.append(f"d={d}: lam_gen={b-r:.4f} lam_mu(K)={bm-r:.4f}")
    print(f"r={r:.4f}:", " | ".join(out))
