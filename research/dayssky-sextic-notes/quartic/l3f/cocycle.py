# Path-independence of the h0-factor: Z(v) = N^{phi_w(v)} Z_w(M_w v); s_d at chamber w multiplies by N^{-(M_w v)_y}, s_c by nothing.
import numpy as np, itertools
sc=np.array([[-1,0],[3,1]]); sd=np.array([[1,1],[0,-1]])
def apply(word):
    M=np.eye(2,dtype=int); phi=np.zeros(2)
    for g in word:            # apply generators left to right in time
        if g=='d': phi=phi-M[1]; M=sd@M
        else: M=sc@M
    return M,phi
seen={}; bad=0
for L in range(0,13):
    for word in itertools.product('cd',repeat=L):
        if any(word[i]==word[i+1] for i in range(L-1)): continue
        M,phi=apply(word); k=tuple(M.flatten())
        if k in seen:
            if not np.allclose(seen[k][0],phi): bad+=1; 
        else: seen[k]=(phi,''.join(word))
print("chambers:",len(seen),"inconsistent paths:",bad)
for k,(phi,w) in sorted(seen.items(),key=lambda t:len(t[1][1])):
    M=np.array(k).reshape(2,2); print("%-7s twist=%-8s phi=%s"%(w or 'id',tuple(M@np.array([1,0])),phi))
M,phi=apply('cd'*6); print("(s_c s_d)^6 -> M=",M.flatten(),"phi=",phi)
