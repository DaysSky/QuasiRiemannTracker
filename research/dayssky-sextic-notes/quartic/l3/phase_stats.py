import sys, numpy as np
def load(fn):
    with open(fn,'rb') as f:
        NR,J=np.frombuffer(f.read(8),dtype=np.int32)
        rows=np.frombuffer(f.read(8*6*NR),dtype=np.float64).reshape(NR,6)
        win=np.frombuffer(f.read(),dtype=np.float64).reshape(J,NR,3)
    return NR,J,rows,win
for fn in sys.argv[1:]:
    NR,J,rows,win=load(fn)
    N,sq,a,b,gr,gi=rows.T
    ok=(sq==1)&(N>1)
    g=(gr+1j*gi)/np.sqrt(N); al=(a+1j*b)/np.sqrt(N); T=al*g
    print(fn,"squarefree rows",ok.sum())
    for j in range(J):
        R=win[j,:,0]+1j*win[j,:,1]; R=R[ok]; d=win[j,ok,2]
        S=np.sum(np.abs(R)**2)
        out=[]
        for name,w in (("1",np.ones(ok.sum())),("g",g[ok]),("gbar",g[ok].conj()),("T",T[ok]),("Tbar",T[ok].conj()),("alpha",al[ok]),("alphabar",al[ok].conj()),("g^2",g[ok]**2),("gbar^2",g[ok].conj()**2)):
            C=np.sum(w*R**2)/S
            out.append("%s:%.3f"%(name,abs(C)))
        print("  window",j," |sum w R^2|/sum|R|^2 :"," ".join(out),"  null~%.3f"%(np.sqrt(np.sum(np.abs(R)**4))/S))
