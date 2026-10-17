import sys; sys.path.insert(0,'../c61')
import numpy as np
from fastR import R, short
th=0.8748389; kap=2*th-1
for dl in (0.44,0.455,0.466,0.475,0.49):
    p=0.5; q=dl*p
    lng=lambda t: 1-dl+(5/6-dl)*(t-1)
    ts=np.linspace(1,1.5,5001)
    vals=[max(short(t,dl,p,kap,'new'),lng(t)) for t in ts]
    i=int(np.argmin(vals)); t=ts[i]
    # find r achieving short
    lo,hi=max(t-.5,0),min(t,1.0); rs=np.linspace(lo,hi,20001)
    AI=1-q-(dl-q)*rs; m=t-rs
    B=14200/14361-1e-5; Lam=19312/4787
    zold=np.maximum(1-2*m,0)/(6*kap); znew=(B-2*m)/Lam
    inW=(m>=(B-.04*Lam)/2)&(m<=min(.44,(B-.034*Lam)/2,(1-.034*6*kap)/2))
    z=np.where(inW,np.maximum(zold,znew),zold)
    AP=1-2*dl*m-2*q*z
    j=int(np.argmax(np.minimum(AI,AP)))
    print(f"dl={dl} R={R(dl,p,kap,'new'):.6f} t*={t:.4f} short={short(t,dl,p,kap,'new'):.6f} long={lng(t):.6f} r*={rs[j]:.4f} m*={m[j]:.4f} AI={AI[j]:.5f} AP={AP[j]:.5f} z={z[j]:.4f}")
