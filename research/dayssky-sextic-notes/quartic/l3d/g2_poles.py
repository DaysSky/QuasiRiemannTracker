import numpy as np, itertools
sc=np.array([[-1,0],[3,1]]); sd=np.array([[1,1],[0,-1]])
# enumerate group
G={(1,0,0,1):np.eye(2,dtype=int)}; frontier=[np.eye(2,dtype=int)]; words={(1,0,0,1):""}
while frontier:
    nf=[]
    for M in frontier:
        for nm,g in (("c",sc),("d",sd)):
            N=g@M; k=tuple(N.flatten())
            if k not in G: G[k]=N; words[k]=nm+words[tuple(M.flatten())]; nf.append(N)
    frontier=nf
print("group order",len(G))
twist0=np.array([1,0])
poles=set()
for k,M in G.items():
    t=M@twist0
    for node,idx in (("c",0),("d",1)):
        if t[idx]==0:
            # pole of chamber series at (Mv)_idx = 1/4  -> hyperplane  M[idx]·v = 1/4
            row=tuple(M[idx]); poles.add((node,row,words[k]))
hyp={}
for node,row,w in poles:
    key=row; hyp.setdefault(key,[]).append((node,w))
print("polar hyperplanes  a*x+b*y = 1/4  (x=s-1/2, y=w-1/2):")
for (a,b),v in sorted(hyp.items()):
    # does it separate (1/2,1/2)+ (abs conv corner) from center (0,0)?  value at (0.5,0.5) and (0,0)
    v1=a*0.5+b*0.5; v0=0.0
    sep=(v1-0.25)*(v0-0.25)<0
    # intersection with residue line y=1/4: a x + b/4 = 1/4 -> x=(1/4-b/4)/a
    xi=None if a==0 else (0.25-b*0.25)/a
    print("  (%2d,%2d): from chambers %s ; at abs-conv corner %.2f ; separates corner from centre: %s ; meets y=1/4 at x=%s"%(a,b,v,v1,sep,"%.3f"%xi if xi is not None else "never"))
