import numpy as np, sys
for f in sys.argv[1:]:
    raw=open(f,'rb').read(); NR=np.frombuffer(raw[:4],dtype=np.int32)[0]
    a=np.frombuffer(raw[8:],dtype=np.float64).reshape(NR,4)
    N,r1,r0,d=a.T; x1=r1/d; x0=r0/d
    def st(x):
        return "mean %.3f  var %.3f (Exp:1)  max %.2f  P(>3) %.4f (Exp %.4f)  P(>5) %.4f (Exp %.4f)"%(x.mean(),x.var(),x.max(),(x>3).mean(),np.exp(-3),(x>5).mean(),np.exp(-5))
    print(f,NR); print(" twisted  ",st(x1)); print(" untwisted",st(x0))
    i=np.argsort(-x0)[:6]; print(" top untwisted rows N:",N[i].astype(int), np.round(x0[i],2))
    i=np.argsort(-x1)[:6]; print(" top twisted rows N:",N[i].astype(int), np.round(x1[i],2))
    # by norm bands
    for lo,hi in [(0,0.1),(0.1,0.3),(0.3,1.01)]:
        Q=N.max(); s=(N>lo*Q)&(N<=hi*Q); print("  band",lo,hi,"rows",s.sum(),"twisted mean %.3f untw mean %.3f"%(x1[s].mean(),x0[s].mean()))
