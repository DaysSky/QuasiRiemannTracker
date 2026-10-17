import numpy as np, sys, math
exec(open('theta_fit.py').read().split("cases=")[0])
Zs=np.logspace(5,6,9)
e4=((N.astype(np.int64)-1)//4)%4          # chi_c(i) = i^{(N(c)-1)/4}
def K1u(idx,Z,u):   # index i^u * m
    ex=np.zeros(len(N),np.int64); zero=np.zeros(len(N),bool)
    for k in idx: v=E[:,k]; zero|=(v==255); ex+=np.where(v==255,0,v)
    ex=ex+u*e4
    z=np.array([1,1j,-1,-1j])[ex&3]; z[zero]=0; w=Wt(N/Z)
    return (g*z*w).sum()
for name,idx in [("1",[]),("pi17",[5]),("pi41",[11]),("pi5",[0])]:
    vals=[np.array([K1u(idx,Z,u)/Z**0.75 for Z in Zs]).mean() for u in range(4)]
    print("%-6s index i^u*m, u=0..3: "%name+"  ".join("%s"%np.round(v,4) for v in vals))
