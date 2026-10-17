# Root functionals of the (dual) G2 action of S51 §2 and the twist value on each.
# Chamber w has spectral point w.v and twist w.l0; its node functionals are e_x, e_y.
# Root functional attached: f = e_node o w ; twist on that root: f(l0).
import numpy as np
sc=np.array([[-1,0],[3,1]]); sd=np.array([[1,1],[0,-1]])
G={(1,0,0,1):np.eye(2,dtype=int)}; fr=[np.eye(2,dtype=int)]
while fr:
    nf=[]
    for M in fr:
        for g in (sc,sd):
            N=g@M; k=tuple(N.flatten())
            if k not in G: G[k]=N; nf.append(N)
    fr=nf
l0=np.array([1,0])
roots={}
for M in G.values():
    for node in (0,1):
        f=tuple(M[node])                 # functional v -> (M v)[node]
        tw=int(np.dot(M[node],l0))
        if f in roots: assert abs(roots[f][1])==abs(tw), (f,roots[f],tw)
        roots[f]=("short(c)" if node==0 else "long(d)",tw)
pos=sorted([f for f in roots if (f[0],f[1])>(0,0)])
print("positive root functionals a*x+b*y, type, |twist|:")
for f in pos: print("  ",f,roots[f][0],abs(roots[f][1]))
# W-invariance of the normaliser: for each simple reflection, the positive set changes by exactly one sign flip,
# and the flipped root is the simple one whose twist is unchanged in absolute value.
for name,g in (("c",sc),("d",sd)):
    img=set(tuple(np.array(f)@g) for f in pos)
    flipped=[f for f in pos if tuple(-np.array(f)) in img]
    print(name,"-reflection flips:",flipped)
