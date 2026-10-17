import sys, importlib.util, numpy as np
sys.path.insert(0,'../combo/gen')
import fmodel as F
from scipy.optimize import minimize_scalar
def theta_for(B, reg=F.REG):
    F.Bc=B
    G0=F.build(reg)
    def best(th):
        r=minimize_scalar(lambda h:-F.worst(th,h,G0)[0],bounds=(.785,.80),method='bounded',options={'xatol':1e-10})
        return -r.fun,r.x
    lo,hi=.8745,.87486
    for _ in range(36):
        mid=(lo+hi)/2; mm,h=best(mid)
        if mm>=0: hi=mid
        else: lo=mid
    return hi
B0=F.Bc
for dB in (0,1e-4,5e-4,1e-3,2e-3,5e-3):
    print("dB=%.4f  theta=%.9f"%(dB,theta_for(B0+dB)))
