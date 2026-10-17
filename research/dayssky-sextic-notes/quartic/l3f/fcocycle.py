# F-cocycle for a c-node twist of index-conductor F: s_c at chamber w multiplies by F^{-(M_w v)_x}; s_d by nothing.
import numpy as np, itertools
sc=np.array([[-1,0],[3,1]]); sd=np.array([[1,1],[0,-1]])
def apply(word):
    M=np.eye(2,dtype=int); phi=np.zeros(2); psi=np.zeros(2)
    for g in word:
        if g=='d': phi=phi-M[1]; M=sd@M
        else: psi=psi-M[0]; M=sc@M
    return M,phi,psi
seen={}; bad=0
for L in range(0,13):
    for w in itertools.product('cd',repeat=L):
        if any(w[i]==w[i+1] for i in range(L-1)): continue
        M,phi,psi=apply(w); k=tuple(M.flatten())
        if k in seen:
            if not (np.allclose(seen[k][0],phi) and np.allclose(seen[k][1],psi)): bad+=1
        else: seen[k]=(phi,psi,''.join(w))
print("chambers",len(seen),"inconsistent (h0 or F cocycle):",bad)
xp=1/3
for w in ('c','dcdc'):
    M,phi,psi=apply(w); Mi=np.linalg.inv(M); v=Mi@np.array([xp,0.0])
    print("chamber %-5s v=%s  h0-exp phi.v=%.4f  F-exp psi.v=%.4f"%(w,np.round(v,4),phi@v,psi@v))
w1=0.5+1/(24*xp)
print("weights %.4f/%.4f -> F-cost per unit log F = %.4f (7/12=%.4f)"%(w1,1-w1,w1*(apply('c')[2]@(np.linalg.inv(apply('c')[0])@[xp,0]))+(1-w1)*(apply('dcdc')[2]@(np.linalg.inv(apply('dcdc')[0])@[xp,0])),7/12))
