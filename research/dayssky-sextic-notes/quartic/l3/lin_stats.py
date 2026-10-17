import sys, numpy as np
def load(fn):
    with open(fn,'rb') as f:
        NR,J=np.frombuffer(f.read(8),dtype=np.int32)
        rows=np.frombuffer(f.read(8*6*NR),dtype=np.float64).reshape(NR,6)
        win=np.frombuffer(f.read(),dtype=np.float64).reshape(J,NR,3)
    return NR,J,rows,win
rng=np.random.default_rng(1)
for fn in sys.argv[1:]:
    NR,J,rows,win=load(fn)
    N,sq,a,b,gr,gi=rows.T
    ok=(sq==1)&(N>1); n=ok.sum()
    g=(gr+1j*gi)/np.sqrt(N); al=(a+1j*b)/np.sqrt(N)
    print(fn,"squarefree rows",n,"sqrt(rows)=%.1f"%np.sqrt(n))
    for j in range(J):
        R=(win[j,:,0]+1j*win[j,:,1])[ok]; S=np.sum(np.abs(R)**2)
        out=[]
        for name,w in (("1",np.ones(n)),("g",g[ok]),("gbar",g[ok].conj()),("g*al",g[ok]*al[ok]),("g*albar",g[ok]*al[ok].conj()),("gbar*al",g[ok].conj()*al[ok]),("gbar*albar",(g[ok]*al[ok]).conj()),("al",al[ok]),("albar",al[ok].conj()),("al^3",al[ok]**3),("g*al^3",g[ok]*al[ok]**3),("gbar*al^3",g[ok].conj()*al[ok]**3)):
            out.append("%s:%.2f"%(name,abs(np.sum(w*R))/np.sqrt(S)))
        # null: random phases
        nulls=[abs(np.sum(np.exp(2j*np.pi*rng.random(n))*R))/np.sqrt(S) for _ in range(200)]
        print("  win",j," ".join(out),"  null mean %.2f 99%% %.2f"%(np.mean(nulls),np.quantile(nulls,.99)))
