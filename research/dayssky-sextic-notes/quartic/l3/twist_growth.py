import sys, numpy as np
sys.path.insert(0,'..'); sys.path.insert(0,'.')
import gauss4 as g4
def load(fn):
    with open(fn,'rb') as f:
        NR,J=np.frombuffer(f.read(8),dtype=np.int32)
        rows=np.frombuffer(f.read(8*6*NR),dtype=np.float64).reshape(NR,6)
        win=np.frombuffer(f.read(),dtype=np.float64).reshape(J,NR,3)
    return NR,J,rows,win
NR,J,rows,win=load(sys.argv[1])
N,sq,a,b,gr,gi=rows.T
ok=np.where((sq==1)&(N>1))[0]
g=(gr+1j*gi)/np.sqrt(N)
facs=[g4.gfactor((int(a[i]),int(b[i]))) for i in ok]
# twists h0: all nonzero Gaussian integers up to units with N(h0)<=Hmax, fourth-power-free (just take all, small)
Hmax=int(sys.argv[2]) if len(sys.argv)>2 else 100
hs=[]
for x in range(0,12):
    for y in range(0,12):
        n=x*x+y*y
        if 1<=n<=Hmax and (x>0 or y>0) and not (x==0 and y>0): hs.append((x,y))   # one per associate class roughly (first quadrant)
hs=sorted(set(hs),key=lambda h:h[0]**2+h[1]**2)
def chi_fac(h,fac):
    r=1
    for pi in fac:
        r*=g4.chi_prime(h,pi)
        if r==0: return 0
    return r
res=[]
for j in range(J):
    R=(win[j,:,0]+1j*win[j,:,1])[ok]; Sn=np.sqrt(np.sum(np.abs(R)**2))
    for h in hs:
        c=np.array([chi_fac(h,f) for f in facs])
        X=abs(np.sum(np.conj(g[ok])*np.conj(c)*R))/Sn
        res.append((j,h[0]**2+h[1]**2,X))
res=np.array(res)
for lo,hi in ((1,5),(5,20),(20,50),(50,101)):
    m=(res[:,1]>=lo)&(res[:,1]<hi)
    if m.sum(): print("N(h0) in [%3d,%3d): n=%3d  mean |X|/sqrt(S)=%.3f  rms=%.3f  max=%.3f"%(lo,hi,m.sum(),res[m,2].mean(),np.sqrt((res[m,2]**2).mean()),res[m,2].max()))
