from ser import *
s=Data(sys.argv[1])
Zs=np.unique(np.round(10**np.arange(2.0,math.log10(s.Zmax)+1e-9,0.0625)).astype(int))
sig=s.sigma(Zs)
xs=np.linspace(0.5,1,200001); Wx=W(xs)
def Wh(u): return np.trapezoid(Wx*xs**(u-1),xs)
for name,cf in (("noalpha",s.Gt),("alpha",s.alpha*s.Gt),("conjalpha",np.conj(s.alpha)*s.Gt)):
    S=s.smooth([(s.Nf,cf)],Zs)
    for pole in (0.75,):
        X=Zs**pole*Wh(pole)
        sel=Zs>=1e4
        c=np.sum(S[sel]*X[sel])/np.sum(X[sel]**2)
        res=S-c*X
        print(name,"fitted c=",c, " |c|=",abs(c))
        for i in range(0,len(Zs),4):
            print(f"   Z={Zs[i]:8d}  |S|/sig={abs(S[i])/sig[i]:7.3f}  |S-cZ^.75|/sig={abs(res[i])/sig[i]:7.3f}")
