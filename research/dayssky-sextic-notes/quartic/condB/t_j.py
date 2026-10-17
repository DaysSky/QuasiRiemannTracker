from ser import *
s=Data(sys.argv[1]); ks=[int(x) for x in sys.argv[2].split(',')]
Zs=np.unique(np.round(10**np.arange(1.5,math.log10(s.Zmax)+1e-9,0.125)).astype(int))
sig=s.sigma(Zs)
cols={}
for k in ks:
    for j in (1,2,3):
        for kind in ('B','D'):
            pc=s.series(k,j,kind)
            cols[f"m{s.tp[k][1]}j{j}{kind}"]=s.smooth(pc,Zs)
        cols[f"m{s.tp[k][1]}j{j}Z*"]=zstar(s,s.series(k,j,'D'),Zs,4)
names=list(cols)
print("Z       "+" ".join(f"{n:>9s}" for n in names))
for i,Z in enumerate(Zs):
    print(f"{Z:7d} "+" ".join(f"{abs(cols[n][i])/sig[i]:9.2e}" for n in names))
