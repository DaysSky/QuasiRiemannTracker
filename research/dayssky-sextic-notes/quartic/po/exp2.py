# (a) h0=1 special dual row: Y(1)=sum_c mu(c) gt(c) W(Nc/D)  vs untwisted sum_c gt(c) W  vs sqrt(sum W^2)
# (b) dual mean square T1(Q)=sum_{0<N(h)<=Q} |sum_c b_c conj(chi_c(h))|^2 / diag  for b=mu*gt*W, gt*W, random phase*W
import sys, math, time, json
import numpy as np
from pocore import *
from gsums import *
Dmax=int(float(sys.argv[1])); Q_exp=float(sys.argv[2]) if len(sys.argv)>2 else 2/3
t0=time.time()
kinds,mods,rs,offs,norms,tab=build_primes(Dmax)
G=gens_for(kinds,mods,rs); gp=prime_gauss(kinds,mods,rs,offs,tab,G)
ptr,idx,cn,om=build_columns(norms,Dmax,lo=Dmax/16)
gt,gen=column_gauss(ptr,idx,kinds,mods,rs,offs,tab,G,gp)
mu=(-1.0)**om
print("build",time.time()-t0,"cols",len(cn))
# (a) scan D
for D in [Dmax/8,Dmax/4,Dmax/2,Dmax]:
    W=Wfun(cn/D); s1=np.sum(mu*gt*W); s0=np.sum(gt*W); rw=math.sqrt(np.sum(W*W))
    print(f"D={D:9.0f}  |sum mu gt W|/sqrt(sumW2)={abs(s1)/rw:7.3f}   |sum gt W|/sqrt(sumW2)={abs(s0)/rw:8.3f}   |sum gt W|/D^(3/4)={abs(s0)/D**0.75:.4f}  sumW={np.sum(W):.1f}")
# (b) dual mean square at D=Dmax with rows h, N(h)<=Q
D=Dmax; m=(cn>D/2); W=Wfun(cn/D)
sel=np.nonzero(m)[0]
rng=np.random.default_rng(7)
B=np.stack([mu*gt*W, gt*W, np.exp(2j*np.pi*rng.random(len(cn)))*W],1)[sel]
pp=[ptr[0]]; ii=[]
for c in sel:
    ii+=list(idx[ptr[c]:ptr[c+1]]); pp.append(len(ii))
pp=np.array(pp,np.int64); ii=np.array(ii,np.int64)
Q=int(D**Q_exp); R=math.isqrt(Q)+1
hx=[];hy=[]
for x in range(-R,R+1):
    for y in range(-R,R+1):
        if 0<x*x+y*y<=Q: hx.append(x); hy.append(y)
hx=np.array(hx,np.int64); hy=np.array(hy,np.int64)
cb=np.conj(B)  # sum_c b conj(chi(h)) = conj( sum_c conj(b) chi(h) )
Fre,Fim,dg=kernel(hx,hy,kinds,mods,rs,offs,tab,pp,ii,np.ascontiguousarray(cb.real),np.ascontiguousarray(cb.imag),W[sel]**2)
A=Fre**2+Fim**2
# classify rows: 4th-power class u*lam^e*m^4, square class u*lam^e*m^2
def mul(a,b): return (a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0])
four=set(); sq=set(); Rm=int(Q**0.25)+2
for a in range(-Rm,Rm+1):
    for b in range(-Rm,Rm+1):
        if (a,b)==(0,0): continue
        m2=mul((a,b),(a,b)); m4=mul(m2,m2)
        for u in ((1,0),(0,1),(-1,0),(0,-1)):
            for le in ((1,0),(1,1),(0,2),(-2,2)):
                four.add(mul(mul(u,le),m4))
            for le in ((1,0),(1,1)):
                sq.add(mul(mul(u,le),m2))
cls=np.array([2 if (x,y) in four else (1 if (x,y) in sq else 0) for x,y in zip(hx,hy)])
print(f"(b) D={D} Q={Q} rows={len(hx)} (4th-class {np.sum(cls==2)}, sq-class {np.sum(cls==1)})")
for v,nm in enumerate(['mu*gt (alpha-twisted)','gt (untwisted)','random phase']):
    tot=A[:,v].sum()/dg.sum()
    parts=[A[cls==k,v].sum()/dg.sum() for k in (0,1,2)]
    top=np.argsort(-A[:,v])[:4]
    print(f"   {nm:22s} T1/diag={tot:7.3f}  generic={parts[0]:.3f} sq={parts[1]:.3f} 4th={parts[2]:.3f}  top rows: "+
          ", ".join(f"({hx[i]},{hy[i]}):{A[i,v]/dg[i]:.1f}" for i in top))
