import sys, time, json, math
import numpy as np
from pocore import *
D=int(float(sys.argv[1])); eta=float(sys.argv[2]); nsamp=int(float(sys.argv[3])) if len(sys.argv)>3 else 0
seed=int(sys.argv[4]) if len(sys.argv)>4 else 0
H=int(round(D**eta))
t0=time.time()
kinds,mods,rs,offs,norms,tab=build_primes(D)
ptr,idx,cn,om=build_columns(norms,D)
C=len(cn); W=Wfun(cn/D); w2=W*W
rng=np.random.default_rng(12345+seed)
mu=(-1.0)**om
# random completely multiplicative +-1 on primes
fpr=rng.choice([-1.0,1.0],size=len(kinds))
rmult=np.array([np.prod(fpr[idx[ptr[c]:ptr[c+1]]]) for c in range(C)])
from gsums import gens_for, prime_gauss, column_gauss
G=gens_for(kinds,mods,rs); gp=prime_gauss(kinds,mods,rs,offs,tab,G)
gt,_=column_gauss(ptr,idx,kinds,mods,rs,offs,tab,G,gp)
names=['mu','rand1','sqfree(=1)','conj(gt)','gt','mu*conj(gt)']
A=np.stack([mu+0j,rng.choice([-1.,1.],C)+0j,np.ones(C)+0j,np.conj(gt),gt,mu*np.conj(gt)],1)*W[:,None]
coef_re=np.ascontiguousarray(A.real); coef_im=np.ascontiguousarray(A.imag)
# rows
R=int(math.isqrt(H))+1
# special rows: k = u * lam^e * m^2, e in {0,1}, u unit  (quadratic-or-trivial character class); flag 4th-power class separately
spec={}
lam=(1,1)
def mul(a,b): return (a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0])
Rm=int(H**0.25)+2
for a in range(-Rm,Rm+1):
    for b in range(-Rm,Rm+1):
        if a==0 and b==0: continue
        m2=mul((a,b),(a,b)); nm2=m2[0]**2+m2[1]**2
        if nm2>H: continue
        m4=mul(m2,m2); nm4=nm2*nm2
        for u in ((1,0),(0,1),(-1,0),(0,-1)):
            for e,le in ((0,(1,0)),(1,lam)):
                k=mul(mul(u,le),m2)
                if k[0]**2+k[1]**2<=H: spec.setdefault(k,'sq')
        for u in ((1,0),(0,1),(-1,0),(0,-1)):
            for le in ((1,0),lam,(0,2),mul((0,2),lam)):
                k=mul(mul(u,le),m4)
                if k[0]**2+k[1]**2<=H: spec[k]='4th'
sx=np.array([k[0] for k in spec],np.int64); sy=np.array([k[1] for k in spec],np.int64)
stype=np.array([spec[k]=='4th' for k in spec])
# total number of lattice rows
ntot=sum(2*int(math.isqrt(H-x*x))+1 for x in range(-int(math.isqrt(H)),int(math.isqrt(H))+1))-1
ngen=ntot-len(spec)
specset=set(spec.keys())
if nsamp==0 or nsamp>=ngen:
    gx=[];gy=[]
    for x in range(-R,R+1):
        yy=int(math.isqrt(max(H-x*x,0))) if H-x*x>=0 else -1
        if yy<0: continue
        ys_=np.arange(-yy,yy+1)
        gx.append(np.full(len(ys_),x)); gy.append(ys_)
    gx=np.concatenate(gx); gy=np.concatenate(gy)
    keep=np.ones(len(gx),bool)
    keep&=~((gx==0)&(gy==0))
    # remove special
    sp=np.array(list(specset)) if specset else np.zeros((0,2),int)
    key=gx*(4*R+3)+gy; skey=sp[:,0]*(4*R+3)+sp[:,1]
    keep&=~np.isin(key,skey)
    gx=gx[keep];gy=gy[keep]; wgen=1.0
else:
    s=0; gx=[];gy=[]
    while s<nsamp:
        x=rng.integers(-R,R+1,size=2*nsamp); y=rng.integers(-R,R+1,size=2*nsamp)
        m=(x*x+y*y<=H)&((x!=0)|(y!=0))
        x=x[m];y=y[m]
        key=x*(4*R+3)+y
        sp=np.array(list(specset)); skey=sp[:,0]*(4*R+3)+sp[:,1]
        m=~np.isin(key,skey); x=x[m];y=y[m]
        gx.append(x);gy.append(y); s+=len(x)
    gx=np.concatenate(gx)[:nsamp].astype(np.int64); gy=np.concatenate(gy)[:nsamp].astype(np.int64)
    wgen=ngen/nsamp
t1=time.time()
Fs_re,Fs_im,ds=kernel(sx,sy,kinds,mods,rs,offs,tab,ptr,idx,coef_re,coef_im,w2)
Fg_re,Fg_im,dgn=kernel(gx,gy,kinds,mods,rs,offs,tab,ptr,idx,coef_re,coef_im,w2)
t2=time.time()
As=Fs_re**2+Fs_im**2; Ag=Fg_re**2+Fg_im**2
S_spec4=As[stype].sum(0); S_specq=As[~stype].sum(0); S_gen=wgen*Ag.sum(0)
Dg=wgen*dgn.sum()+ds.sum()
S=S_spec4+S_specq+S_gen
res=dict(D=D,eta=eta,H=H,rows=ntot,nspec=len(spec),ngen=ngen,nsamp=len(gx),cols=C,sumW2=float(w2.sum()),
         diag=float(Dg),time_build=t1-t0,time_kernel=t2-t1)
print(json.dumps(res))
# standard error of generic part (sampling)
for v,nm in enumerate(names):
    se = wgen*math.sqrt(len(gx))*Ag[:,v].std() if wgen>1 else 0.0
    print(f"{nm:11s} S/diag={S[v]/Dg:8.4f} (gen {S_gen[v]/Dg:7.4f} +- {se/Dg:.4f}, sq-rows {S_specq[v]/Dg:7.4f}, 4th-rows {S_spec4[v]/Dg:7.4f})  S/((H+D)*sumW2)={S[v]/((H+D)*w2.sum()):7.4f}  maxrow/D={math.sqrt(max(As[:,v].max(),Ag[:,v].max()))/D:.4f}")
np.savez(f"gout_D{D}_eta{eta:.4f}_n{nsamp}_s{seed}.npz",S=S,S_gen=S_gen,S_specq=S_specq,S_spec4=S_spec4,diag=Dg,names=names,
         gx=gx,gy=gy,Ag=Ag,sx=sx,sy=sy,As=As,stype=stype,res=json.dumps(res))
