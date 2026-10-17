from ser import *
s=Data(sys.argv[1])
Zs=np.unique(np.round(10**np.arange(3.0,math.log10(s.Zmax)+1e-9,0.03125)).astype(int))
sig=s.sigma(Zs)
wins=[(1e3,1e4),(1e4,1e5),(1e5,5e5),(5e5,2e6)]
def rms(S,a,b):
    sel=(Zs>=a)&(Zs<b); return np.mean((np.abs(S[sel])/sig[sel])**2)
rng=np.random.default_rng(7)
R=[]
for t in range(40):
    ph=np.exp(2j*np.pi*rng.random(s.n)); S=s.smooth([(s.Nf,ph)],Zs); R.append([rms(S,a,b) for a,b in wins])
R=np.array(R); print("random (40 draws): mean of R^2 per window",np.round(R.mean(0),3)," sd of single-draw",np.round(R.std(0),3))
for j in (1,2,3):
    rows=[]
    for k in range(s.K):
        S=s.smooth(s.series(k,j,'B'),Zs); rows.append((s.tp[k][1],[rms(S,a,b) for a,b in wins]))
    small=np.array([v for N,v in rows if N<=50]); big=np.array([v for N,v in rows if N>50])
    print(f"j={j}: N(m)<=50 ({len(small)} m) mean R^2:",np.round(small.mean(0),3),f"| N(m)>50 ({len(big)} m):",np.round(big.mean(0),3))
