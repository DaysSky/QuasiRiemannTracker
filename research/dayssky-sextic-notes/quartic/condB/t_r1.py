from ser import *
s=Data(sys.argv[1])
Zs=np.unique(np.round(10**np.arange(1.0,math.log10(s.Zmax)+1e-9,0.125)).astype(int))
sig=s.sigma(Zs)
pc=s.series(0,0)
SB=s.smooth(pc,Zs)
out={'B':SB}
for tw in (4,-4,0):
    out[f'Z*tw{tw}']=zstar(s,pc,Zs,tw)
# also untwisted (no alpha): coefficients Gt only
pu=[(s.Nf,s.Gt)]
out['noalpha']=s.smooth(pu,Zs)
print("Z      sigma   "+"  ".join(f"{k:>9s}" for k in out))
for i,Z in enumerate(Zs):
    print(f"{Z:7d} {sig[i]:7.2f}  "+"  ".join(f"{abs(v[i])/sig[i]:9.2e}" for v in out.values()))
