# Rank-1 check, SYMBOLIC in (r,a): the true Kubota local factor at pi^l || index (Gauss sums g(pi^l,pi^k;eps^t)),
# expanded as in CG (5.12)-(5.21) with the self-symbol sign S(k') = s^((K'+1)k'), satisfies f'(x) = (qx)^(l-K) f'(1/(q^2x)).
from ratf import F, sub, R, A, X1
r=F(R); a=F(A); x=F(X1); q=r*r
for s in (1,-1):
  graw={0:F(-1),1:a,2:r.inv(),3:F(s)/(q*a)}
  for t in (1,3):
    for useS in (True,False):
      res=[]
      for l in range(0,13):
        Nc={}
        for k in range(0,l+2):
            if k<=l: Nc[k]=((q**k-q**(k-1)) if k>0 else F(1)) if k%4==0 else F(0)
            else: Nc[k]=q**l*q*graw[(t*(l+1))%4] if (t*(l+1))%4 else -(q**l)
        Nk={k:sum((v*x**e for e,v in Nc.items() if e%4==k),F(0)) for k in range(4)}
        K=lambda k:(l-2*k)%4
        ok=True
        for k in range(4):
            kp=(l+1-k)%4; Kp=K(kp); f=Nk[k]
            if Kp!=3:
                S=s**(((Kp+1)*kp)%2) if useS else 1
                f=f-F(S)*q**Kp*q*graw[(t*(Kp+1))%4]*x**(Kp+1)*Nk[kp]
            f=f/(1-q**3*x**4)
            if f.iszero(): continue
            if not (f-(q*x)**(l-K(k))*sub(f,1)).iszero(): ok=False
        res.append('T' if ok else 'F')
      print("s=%+d t=%d S=%-5s l=0..12: %s"%(s,t,useS,''.join(res)))
