import numpy as np, sys
exec(open('theta_fit.py').read().split("cases=")[0])
Zs=np.logspace(5,6,9)
for name,idx in [("(1,-4)",[5]),("(1,4)",[6]),("(-1,-2)",[0]),("(-1,2)",[1]),("(3,-2)",[3]),("(3,2)",[4]),("(-5,-2) N29",[7]),("(-5,2) N29",[8]),("(-1,-6) N37",[9]),("(-1,6) N37",[10]),("(5,-4) N41",[11]),("(5,4) N41",[12]),("(-7,0) N49",[13]),("(7,-2) N53",[14]),("(7,2) N53",[15])]:
    v=np.array([K1(idx,Z)[0]/Z**0.75 for Z in Zs])
    print("%-12s mean|.|=%.4f  mean value=%s"%(name,np.abs(v).mean(),np.round(v.mean(),4)))
