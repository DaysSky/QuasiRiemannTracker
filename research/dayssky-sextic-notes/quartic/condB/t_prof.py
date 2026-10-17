from ser import *
s=Data(sys.argv[1])
Zs=np.unique(np.round(10**np.arange(2.0,math.log10(s.Zmax)+1e-9,0.03125)).astype(int))
sig=s.sigma(Zs)
edges=10**np.arange(2.0,6.31,0.5)
def prof(S):
    out=[]
    for a,b in zip(edges[:-1],edges[1:]):
        sel=(Zs>=a)&(Zs<b); out.append(np.mean((np.abs(S[sel])/sig[sel])**2) if sel.any() else np.nan)
    return np.array(out)
print("window lower edges:", " ".join(f"{e:8.0e}" for e in edges[:-1]))
seen=set()
for k in range(s.K):
    (a,b),Nm,_=s.tp[k]
    if Nm>50 or (a,-b) in seen: continue
    seen.add((a,b))
    for j in (1,2,3):
        for kind in ('B','D'):
            p=prof(s.smooth(s.series(k,j,kind),Zs))
            print(f"m={str((a,b)):9s} N={Nm:3d} j={j} {kind}: "+" ".join(f"{x:8.2f}" for x in p))
p=prof(s.smooth(s.series(0,0),Zs)); print("r=1               : "+" ".join(f"{x:8.2f}" for x in p))
