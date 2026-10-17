# saving eta over copies only for t in [t1, t2]
import numpy as np, capc4, sys
ND, NP, NT, NR = capc4.ND, capc4.NP, capc4.NT, capc4.NR
def make(eta, t1, t2):
    def Rgrid(kap, cP, cI, sup, al):
        dl = np.linspace(0, kap, ND)[:,None,None,None]; p = np.linspace(0, .5, NP)[None,:,None,None]
        t = np.linspace(1, 1.5, NT)[None,None,:,None]; a = np.linspace(0, 1, NR)[None,None,None,:]
        lo = np.maximum(t - .5, 0); hi = np.minimum(t, 1.0); r = lo + (hi - lo)*a
        q = dl*p; m = t - r
        AI = 1 - dl*r - 2*q*np.minimum((1 - r)/cI, sup)
        AP = 1 - 2*dl*m - 2*q*np.minimum(np.maximum(1 - 2*m, 0)/(cP*kap), sup)
        short = np.minimum(AI, AP).max(axis=3)
        tt = t[...,0]
        loss = (5/6)*(tt - 1) - np.where((tt >= t1) & (tt <= t2), eta, 0.0)
        lng = 1 + np.maximum(loss, 0) - dl[...,0]*tt
        R = np.maximum(short, lng).min(axis=2)
        return dl[:,:,0,0], p[:,:,0,0], R
    return Rgrid
print("NT grid points in t:", capc4.NT, flush=True)
for eta, t1, t2 in [(0.003, 1.14, 1.18), (0.003, 1.10, 1.25), (0.003, 1.05, 1.5)]:
    capc4.Rgrid = make(eta, t1, t2)
    print(f"saving {eta} only on t in [{t1},{t2}]: theta = {capc4.theta(4.498/0.7497, 2, 5/6):.5f}", flush=True)
