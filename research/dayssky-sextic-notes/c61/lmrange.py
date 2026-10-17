# theta when LM (slope a1) holds on s in [1, s0] and copies slope 5/6 beyond (capc4 grid model, current slot costs)
import numpy as np, capc4
ND, NP, NT, NR = capc4.ND, capc4.NP, capc4.NT, capc4.NR
def make_Rgrid(s0, a1):
    def Rgrid(kap, cP, cI, sup, al):
        dl = np.linspace(0, kap, ND)[:,None,None,None]; p = np.linspace(0, .5, NP)[None,:,None,None]
        t = np.linspace(1, 1.5, NT)[None,None,:,None]; a = np.linspace(0, 1, NR)[None,None,None,:]
        lo = np.maximum(t - .5, 0); hi = np.minimum(t, 1.0); r = lo + (hi - lo)*a
        q = dl*p; m = t - r
        AI = 1 - dl*r - 2*q*np.minimum((1 - r)/cI, sup)
        AP = 1 - 2*dl*m - 2*q*np.minimum(np.maximum(1 - 2*m, 0)/(cP*kap), sup)
        short = np.minimum(AI, AP).max(axis=3)
        tt = t[...,0]
        slope_pt = np.where(tt <= s0, a1*(tt-1), a1*(s0-1) + (5/6)*(tt-s0))   # integrated LM loss
        lng = 1 + slope_pt - dl[...,0]*tt
        R = np.maximum(short, lng).min(axis=2)
        return dl[:,:,0,0], p[:,:,0,0], R
    return Rgrid
if __name__ == "__main__":
    for s0 in (1.0, 1.05, 1.1, 1.16, 1.2, 1.3, 1.5):
        capc4.Rgrid = make_Rgrid(s0, 0.0)
        print(f"LM (slope 0) on [1, {s0}], 5/6 beyond: theta = {capc4.theta(4.498/0.7497, 2, 5/6):.5f}", flush=True)
