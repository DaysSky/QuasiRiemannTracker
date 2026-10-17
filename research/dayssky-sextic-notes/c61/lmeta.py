# theta when the copies bound is beaten by a uniform factor U^{-eta} beyond the conductor:
# e(t) = max(1, 1 + (5/6)(t-1) - eta)   (capc4 grid model, current slot costs; baseline 0.87492)
import numpy as np, capc4, sys
ND, NP, NT, NR = capc4.ND, capc4.NP, capc4.NT, capc4.NR
def make(eta):
    def Rgrid(kap, cP, cI, sup, al):
        dl = np.linspace(0, kap, ND)[:,None,None,None]; p = np.linspace(0, .5, NP)[None,:,None,None]
        t = np.linspace(1, 1.5, NT)[None,None,:,None]; a = np.linspace(0, 1, NR)[None,None,None,:]
        lo = np.maximum(t - .5, 0); hi = np.minimum(t, 1.0); r = lo + (hi - lo)*a
        q = dl*p; m = t - r
        AI = 1 - dl*r - 2*q*np.minimum((1 - r)/cI, sup)
        AP = 1 - 2*dl*m - 2*q*np.minimum(np.maximum(1 - 2*m, 0)/(cP*kap), sup)
        short = np.minimum(AI, AP).max(axis=3)
        tt = t[...,0]
        loss = np.maximum(0, (5/6)*(tt - 1) - eta)
        lng = 1 + loss - dl[...,0]*tt
        R = np.maximum(short, lng).min(axis=2)
        return dl[:,:,0,0], p[:,:,0,0], R
    return Rgrid
for eta in [float(e) for e in sys.argv[1:]]:
    capc4.Rgrid = make(eta)
    print(f"uniform saving U^-{eta} over copies: theta = {capc4.theta(4.498/0.7497, 2, 5/6):.5f}", flush=True)
