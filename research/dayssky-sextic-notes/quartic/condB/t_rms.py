from ser import *
s=Data(sys.argv[1])
Zs=np.unique(np.round(10**np.arange(2.0,math.log10(s.Zmax)+1e-9,0.0625)).astype(int))
sig=s.sigma(Zs)
wins=[(1e2,1e3),(1e3,1e4),(1e4,1e5),(1e5,2e6)]
def rms(S,a,b):
    sel=(Zs>=a)&(Zs<b); return math.sqrt(np.mean((np.abs(S[sel])/sig[sel])**2))
rng=np.random.default_rng(5)
print("series            "+"  ".join(f"[{a:.0e},{b:.0e})" for a,b in wins))
# random control
for t in range(3):
    ph=np.exp(2j*np.pi*rng.random(s.n))
    S=s.smooth([(s.Nf,ph)],Zs); print(f"random{t}           "+"  ".join(f"{rms(S,a,b):12.3f}" for a,b in wins))
S=s.smooth(s.series(0,0),Zs); print(f"r=1               "+"  ".join(f"{rms(S,a,b):12.3f}" for a,b in wins))
agg={j:[] for j in (1,2,3)}
for k in range(s.K):
    Nm=s.tp[k][1]
    for j in (1,2,3):
        S=s.smooth(s.series(k,j,'B'),Zs)
        v=[rms(S,a,b) for a,b in wins]; agg[j].append(v)
        if k in (0,3,5,9,20,44): print(f"N(m)={Nm:4d} j={j} B   "+"  ".join(f"{x:12.3f}" for x in v))
for j in (1,2,3):
    a=np.array(agg[j]); print(f"all 45 m, j={j} RMS "+"  ".join(f"{x:12.3f}" for x in np.sqrt(np.mean(a**2,axis=0))))
