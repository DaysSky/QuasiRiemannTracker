# Group of functional equations for the 3-variable series
#   Z_T(x1,x2,y) = sum_m N(m)^{-1/2-y} psi_a(m^3,1/2+x1) conj psi_a(m^3,1/2+conj x2)
# K1, K2: Kubota FE in s1, s2 (conductor N(m)^3); P: Poisson in m (T-form <-> Mobius form);
# J: joint FE of 1/L(1/2+a1, xi_k) 1/L(1/2+a2, conj xi_k) in the Mobius form (conductor N(k)).
import numpy as np, itertools
K1=np.array([[-1,0,0],[0,1,0],[3,0,1]])
K2=np.array([[1,0,0],[0,-1,0],[0,3,1]])
P =np.array([[1,0,1],[0,1,1],[0,0,-1]])
J =np.array([[0,-1,0],[-1,0,0],[-1,-1,1]])
JT=P@J@P
print("JT =",JT.tolist())
for name,M in [("K1",K1),("K2",K2),("P",P),("J",J),("JT",JT)]:
    print(name,"involution:",np.array_equal(M@M,np.eye(3,dtype=int)))
# invariant symmetric form B with M^T B M = B for K1,K2,JT
import sympy as sp
b=sp.symbols('b0:6'); B=sp.Matrix([[b[0],b[1],b[2]],[b[1],b[3],b[4]],[b[2],b[4],b[5]]])
eqs=[]
for M in (K1,K2,JT):
    Ms=sp.Matrix(M); eqs+=list(Ms.T*B*Ms-B)
sol=sp.solve(eqs,b,dict=True); print("invariant forms:",sol)
Bs=B.subs(sol[0]); print("B =",Bs, " eigen:",[sp.nsimplify(e) for e in Bs.subs({s:1 for s in Bs.free_symbols}).eigenvals()])
# orders of pairwise products
def order(M,maxn=50):
    A=np.eye(3,dtype=int)
    for n in range(1,maxn+1):
        A=A@M
        if np.array_equal(A,np.eye(3,dtype=int)): return n
    return None
for (n1,M1),(n2,M2) in itertools.combinations([("K1",K1),("K2",K2),("JT",JT)],2):
    M=M1@M2; print(n1+n2,"order",order(M),"trace",np.trace(M),"eig",np.round(np.linalg.eigvals(M),4))
# group growth
seen={tuple(np.eye(3,dtype=int).flatten())}; frontier=[np.eye(3,dtype=int)]
for depth in range(1,9):
    new=[]
    for A in frontier:
        for G in (K1,K2,JT):
            C=G@A; t=tuple(C.flatten())
            if t not in seen: seen.add(t); new.append(C)
    frontier=new; print("word length",depth,"group elements so far",len(seen))
