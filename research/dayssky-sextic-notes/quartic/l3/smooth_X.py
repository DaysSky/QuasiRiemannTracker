import sys, numpy as np
sys.path.insert(0,'..'); sys.path.insert(0,'.')
def load(fn):
    with open(fn,'rb') as f:
        NR,J=np.frombuffer(f.read(8),dtype=np.int32)
        rows=np.frombuffer(f.read(8*6*NR),dtype=np.float64).reshape(NR,6)
        win=np.frombuffer(f.read(),dtype=np.float64).reshape(J,NR,3)
    return NR,J,rows,win
def bump(x,a,b):
    y=np.zeros_like(x); m=(x>a)&(x<b); t=x[m]
    y[m]=np.exp(-1.0/((t-a)*(b-t))+4.0/(b-a)**2*0+0); return y
rng=np.random.default_rng(5)
for fn in sys.argv[1:]:
    NR,J,rows,win=load(fn)
    N,sq,a,b,gr,gi=rows.T
    ok=(sq==1)&(N>1); Qm=N[ok].max()
    g=((gr+1j*gi)/np.sqrt(N))[ok]; Nn=N[ok]
    print(fn,"rows",ok.sum(),"Qmax",Qm)
    for j in range(J):
        R=(win[j,:,0]+1j*win[j,:,1])[ok]
        for (lo,hi) in ((0.5,1.0),(0.25,1.0),(0.3,0.9)):
            V=bump(Nn/Qm,lo,hi); sh=((Nn/Qm>lo)&(Nn/Qm<=hi)).astype(float)
            for name,w in (("gbar",g.conj()),("g",g)):
                Xs=abs(np.sum(V*w*R)); Xh=abs(np.sum(sh*w*R))
                ns=np.sqrt(np.sum(V**2*np.abs(R)**2)); nh=np.sqrt(np.sum(sh*np.abs(R)**2))
                nulls=[abs(np.sum(V*np.exp(2j*np.pi*rng.random(len(R)))*R))/ns for _ in range(100)]
                print("  win %d V=(%.2f,%.2f) %-4s smooth |X|/null=%.3f   sharp |X|/null=%.3f   (random-phase smooth: mean %.2f, 5%% %.2f)"%(j,lo,hi,name,Xs/ns,Xh/nh,np.mean(nulls),np.quantile(nulls,.05)))
