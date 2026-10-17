import numpy as np, sys
exec(open('theta_fit.py').read().split("cases=")[0])
Zs=np.logspace(5,6,9)
for name,idx in [("pi17",[5]),("pi17^2",[5,5]),("pi17^3",[5,5,5]),("pi41",[11]),("pi41^3",[11,11,11]),("pi5^3",[0,0,0]),("pi17*pi41",[5,11]),("(pi17*pi41)^3",[5,5,5,11,11,11])]:
    v=np.array([K1(idx,Z)[0]/Z**0.75 for Z in Zs]).mean()
    print("%-14s mean K1/Z^{3/4} = %s  |.|=%.4f"%(name,np.round(v,4),abs(v)))
