import numpy as np, sys
exec(open('theta_fit.py').read().split("cases=")[0])   # reuse loader + K1
Zs=np.logspace(5,6,9)
cases=[("pi5^2",[0,0]),("pi13^2",[3,3]),("pi17^2",[5,5]),("3^2",[2,2]),("(pi5 pi13)^2",[0,0,3,3]),("pi5 conj(pi5)=5",[0,1]),
       ("pi5",[0]),("pi13",[3]),("pi17",[5]),("3",[2]),("pi5^3",[0,0,0]),("pi5*pi13",[0,3]),("pi5^2*pi13",[0,0,3])]
for name,idx in cases:
    v=np.array([abs(K1(idx,Z)[0])/Z**0.75 for Z in Zs])
    print("%-16s mean |K1|/Z^{3/4} on [1e5,1e6] = %.4f  (min %.4f max %.4f)   noise ~%.4f"%(name,v.mean(),v.min(),v.max(),0.155*1e5**-0.25*0.75))
