import numpy as np, sys, itertools
from scipy.optimize import linprog
exec(open('tits_grh.py').read().split("def hull_min")[0].replace("L=int(sys.argv[1]) if len(sys.argv)>1 else 10","L=8"))
# keep word labels
words={tuple(np.eye(3).round(6).flatten()):''}; fr=[(np.eye(3),'')]
GL=[(np.eye(3),'')]
for _ in range(8):
    nf=[]
    for g,wd in fr:
        for i in (1,2,3):
            h=S[i]@g; k=tuple(h.round(6).flatten())
            if k not in words: words[k]='s%d'%i+wd; nf.append((h,'s%d'%i+wd)); GL.append((h,'s%d'%i+wd))
    fr=nf
pts=[];lab=[]
for g,wd in GL:
    pts.append(g@rho); lab.append(('abs-apex',wd))
    pts.append(g@p0); lab.append(('grh-apex',wd))
rays=[];rlab=[]
for g,wd in GL:
    for r,nm in zip(R_abs,['e1','e2','e3']): rays.append(g@r); rlab.append(('abs',nm,wd))
    for r,nm in zip(R_grh,['e1','e2','f']): rays.append(g@r); rlab.append(('grh',nm,wd))
pts=np.array(pts); rays=np.array(rays)
t=float(sys.argv[1]) if len(sys.argv)>1 else 0.5
n,m=len(pts),len(rays)
A=np.vstack([np.hstack([pts.T,rays.T]),np.hstack([np.ones(n),np.zeros(m)])])
b=np.array([0,0,t,1.0])
# minimise total ray weight to get a sparse certificate
c=np.concatenate([np.zeros(n),np.ones(m)])
res=linprog(c,A_eq=A,b_eq=b,bounds=(0,None),method='highs')
print("status",res.status)
x=res.x
for i in np.nonzero(x[:n]>1e-9)[0]: print("point",lab[i],np.round(pts[i],4),"weight",round(x[i],6))
for j in np.nonzero(x[n:]>1e-9)[0]: print("ray",rlab[j],np.round(rays[j],4),"weight",round(x[n+j],6))
