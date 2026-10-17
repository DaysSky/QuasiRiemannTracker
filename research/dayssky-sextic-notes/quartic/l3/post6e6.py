import numpy as np, math
raw=open('rows_Z6000000.bin','rb').read(); NR=np.frombuffer(raw[:4],dtype=np.int32)[0]
a=np.frombuffer(raw[8:],dtype=np.float64).reshape(NR,4); N,r1,r0,d=a.T; x=r1/d; y=r0/d; Ni=N.astype(np.int64)
pairs=NR/2
print("rows %d  S_tw/diag %.4f (s.e. %.4f)  S_untw/diag %.4f"%(NR,r1.sum()/d.sum(),math.sqrt(x.var()/pairs),r0.sum()/d.sum()))
print("var %.3f  E[x^2]/E[x]^2 %.3f  E[x^3]/E[x]^3 %.2f  P(>3) %.4f (Exp %.4f)  P(>5) %.4f (Exp %.4f)"%(x.var(),(x**2).mean()/x.mean()**2,(x**3).mean()/x.mean()**3,(x>3).mean(),math.exp(-3),(x>5).mean(),math.exp(-5)))
mx=x.max(); print("max %.2f at N=%d; P(max>=obs | Exp(1), %d pairs) = %.3f; ln(pairs)=%.2f"%(mx,Ni[x.argmax()],pairs,1-(1-math.exp(-mx))**pairs,math.log(pairs)))
for r in (1,5): s=(Ni%8==r); print("  N(m)≡%d mod 8: mean %.3f"%(r,x[s].mean()))
