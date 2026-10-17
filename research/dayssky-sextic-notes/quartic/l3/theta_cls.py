import numpy as np, sys
exec(open('theta_fit.py').read().split("cases=")[0])
Zs=np.logspace(5,6,9)
cls=(N.astype(np.int64)%8)
def K1c(idx,Z,mask):
    ex=np.zeros(len(N),np.int64); zero=np.zeros(len(N),bool)
    for k in idx: v=E[:,k]; zero|=(v==255); ex+=np.where(v==255,0,v)
    z=np.array([1,1j,-1,-1j])[ex&3]; z[zero]=0; w=Wt(N/Z)
    return (g*z*w*mask).sum()
for name,idx in [("index 1",[]),("pi5 (5 mod 8)",[0]),("pi13 (5 mod 8)",[3]),("pi17 (1 mod 8)",[5]),("pi41 (1 mod 8)",[11])]:
    a=np.array([K1c(idx,Z,cls==1)/Z**0.75 for Z in Zs]).mean(); b=np.array([K1c(idx,Z,cls==5)/Z**0.75 for Z in Zs]).mean()
    print("%-16s c with N≡1(8): %s   c with N≡5(8): %s   total: %s"%(name,np.round(a,4),np.round(b,4),np.round(a+b,4)))
