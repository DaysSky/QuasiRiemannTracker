# RMS-over-h0 (vector-valued Bochner) h0-exponent for the block Region_d.
import numpy as np
from blocks import sigma, boundary_d
def blk_rms(x,y,tC,tD,lam):   # tK = 1
    B1=-x*tC-y*tD+0.5*(max(1+tC,tD,2/3*(1+tC+tD))-1)
    B2=-y*tD-(0.5+x)*tC+0.5*sigma(tC,tD,lam)-0.5*min(tD,1.0)
    return min(B1,B2)
def b_rms(x,y,lam,n=80):
    best=-1e9
    for tC in np.linspace(0,30,n):
        for tD in np.linspace(0,(1+tC)/2,n):
            best=max(best,blk_rms(x,y,tC,tD,lam))
    return max(best,0.0)
lam=1.5573
for y in (0,0.1,0.25,0.4,0.5,0.75,1.0):
    xb=boundary_d(y,lam)
    print("y=%.2f x>%.4f  b_rms=%.3f   (b at x+0.05: %.3f, x+0.1: %.3f)"%(y,xb,b_rms(xb+1e-3,y,lam),b_rms(xb+0.05,y,lam),b_rms(xb+0.1,y,lam)))
