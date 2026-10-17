import sys, numpy as np
def load(fn):
    with open(fn,'rb') as f:
        NR,J=np.frombuffer(f.read(8),dtype=np.int32)
        rows=np.frombuffer(f.read(8*6*NR),dtype=np.float64).reshape(NR,6)
        win=np.frombuffer(f.read(),dtype=np.float64).reshape(J,NR,3)
    return NR,J,rows,win
NR,J,rows,w2=load("theta_L2e5.bin"); _,_,rows4,w4=load("theta_L4e5.bin")
N,sq,a,b,gr,gi=rows.T
Y2=w2[0,:,0]+1j*w2[0,:,1]; Y4=w4[0,:,0]+1j*w4[0,:,1]
t2=Y2/2e5**0.75; t4=Y4/4e5**0.75
ok=(sq==1)&(N>1)
big=ok&(np.abs(t4)>np.quantile(np.abs(t4[ok]),0.5))
print("theta extraction consistency (rows with |tau| above median): median |t4/t2| = %.3f, median |arg(t4/t2)| = %.3f rad"%(np.median(np.abs(t4[big]/t2[big])),np.median(np.abs(np.angle(t4[big]/t2[big])))))
tau=t4
g=(gr+1j*gi)/np.sqrt(N); al=(a+1j*b)/np.sqrt(N); T=al*g
nz=ok&(np.abs(tau)>1e-3*np.abs(tau[ok]).max())
print("squarefree rows:",ok.sum()," with non-negligible tau:",nz.sum(),"  (Suzuki support predicts density ~2^-omega)")
def bump(x,lo=0.5,hi=1.0):
    y=np.zeros_like(x); m=(x>lo)&(x<hi); t=x[m]; y[m]=np.exp(-1.0/((t-lo)*(hi-t))); return y
print(" Z0    #c   |C(conjT.tau)|  |C(T.tau)|  random-phase  Cauchy     |C|/rand  Cauchy/rand")
for Z0 in (125,250,500,1000,2000,4000):
    W=bump(N/Z0)*ok
    C1=np.sum(W*np.conj(T)*tau); C2=np.sum(W*T*tau)
    rand=np.sqrt(np.sum(W**2*np.abs(tau)**2)); cs=np.sqrt(np.sum(W**2))*rand
    print("%5d %5d   %10.4f   %10.4f   %10.4f   %9.4f   %7.3f  %7.3f"%(Z0,(W>0).sum(),abs(C1),abs(C2),rand,cs,abs(C1)/rand,cs/rand))
