import sys, numpy as np
def load(fn):
    with open(fn,'rb') as f:
        NR,J=np.frombuffer(f.read(8),dtype=np.int32)
        rows=np.frombuffer(f.read(8*6*NR),dtype=np.float64).reshape(NR,6)
        win=np.frombuffer(f.read(),dtype=np.float64).reshape(J,NR,3)
    return NR,J,rows,win
def bump(x,a=0.5,b=1.0):
    y=np.zeros_like(x); m=(x>a)&(x<b); t=x[m]; y[m]=np.exp(-1.0/((t-a)*(b-t))); return y
xs=np.linspace(0.5,1,20001); Vhat=np.trapezoid(bump(xs)*xs**(-0.25),xs)   # V^(3/4) = int V(x) x^{3/4-1} dx
for fn in sys.argv[1:]:
    NR,J,rows,win=load(fn)
    N,sq,a,b,gr,gi=rows.T
    ok=(sq==1)&(N>1); g=((gr+1j*gi)/np.sqrt(N))[ok]; Nn=N[ok]; Qm=Nn.max()
    print(fn,"Qmax",Qm)
    for j in range(J):
        R=(win[j,:,0]+1j*win[j,:,1])[ok]
        for Q0 in (Qm/8,Qm/4,Qm/2,Qm):
            V=bump(Nn/Q0); X=np.sum(V*g.conj()*R); noise=np.sqrt(np.sum(V**2*np.abs(R)**2))
            c=X/(Q0**0.75*Vhat)
            print("  win %d Q0=%7.0f  rows=%4d  X/(Q0^.75 Vhat) = %8.2f %+8.2fi  |.|=%7.2f   noise/(Q0^.75 Vhat)=%7.2f"%(j,Q0,(V>0).sum(),c.real,c.imag,abs(c),noise/(Q0**0.75*Vhat)))
