# Zero-density for the quartic family {chi_k: N(k)<=H} from large-sieve inputs, then the (L3) exponent it implies.
# Everything in log_H units. Delta(n) = log_H of the large-sieve constant at length H^n.
import numpy as np
def Delta_GZ(n):  # min of classical multiplicative LS (H^2+N) and Gao-Zhao (H+N+(HN)^{3/4})
    return min(max(2.0, n), max(1.0, n, 0.75*(1+n)))
def Delta_opt(n): return max(1.0, n)
def typeI(sig, nlo, nhi, D, jmax=12):
    # max over piece lengths n in [nlo,nhi] of min over powers j of Delta(j n) + j n (1-2 sig)
    worst=-1e9
    for n in np.linspace(nlo, nhi, 200):
        best=min(D(j*n)+j*n*(1-2*sig) for j in range(1,jmax+1))
        worst=max(worst,best)
    return worst
def N_sigma(sig, D):
    best=9e9
    for x in np.linspace(0.0, 1.0, 41):          # mollifier X = H^x
        for y in np.linspace(max(x,0.01), 6.0, 120):   # smoothing Y = H^y
            tI = typeI(sig, x, y, D) if y>x else -9e9
            DII = 0.75*1.5*0 # placeholder
            # type II: Y^{1-2sig} * (sum |L|^4)^{1/2} (sum |M_X|^4)^{1/2}
            L4 = D(1.0)            # fourth moment of L: length H
            M4 = D(2*x)
            tII = y*(1-2*sig) + 0.5*L4 + 0.5*M4
            val=max(tI,tII)
            best=min(best,val)
    return min(best,1.0)
if __name__=="__main__":
    for name,D in (("optimal LS",Delta_opt),("GZ+classical",Delta_GZ)):
        print(name)
        for sig in (0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9):
            n=N_sigma(sig,D); print("  sigma=%.2f  N<=H^%.4f  A=%.3f"%(sig,n,n/(1-sig)))
