# Rank-1 sanity check of the slice derivation: the true Kubota local factor at pi | index (pi^l || index),
#   N(x) = sum_k g(pi^l, pi^k; eps^t) x^k,  g = phi(pi^k) [k<=l, 4|k], q^l g(1,pi;eps^{t(l+1)}) [k=l+1], else 0,
# must make f'_k satisfy f'_k(x) = (qx)^(l-K(k)) f'_k(1/(q^2 x)), with the self-symbol sign S(k') = s^((K'+1)(k'+lam-l)), lam=l.
import sympy as sp, sys
x=sp.symbols('x')
r=sp.Integer(3); q=r*r; a=sp.Rational(2,7)
for s in (1,-1):
  graw={0:sp.Integer(-1),1:a,2:1/r,3:sp.Integer(s)/(q*a)}
  for t in (1,3):
    for Sexp in ('derived','none'):
      res=[]
      for l in range(0,9):
        Nc={}
        for k in range(0,l+2):
            if k<=l: Nc[k]=(q**k-q**(k-1) if k>0 else 1) if k%4==0 else 0
            else: Nc[k]=q**l*q*graw[(t*(l+1))%4] if (t*(l+1))%4 else -q**l
        Nk={k:sum(v*x**e for e,v in Nc.items() if e%4==k) for k in range(4)}
        K=lambda k:(l-2*k)%4
        ok=True
        for k in range(4):
            kp=(l+1-k)%4; Kp=K(kp); f=Nk[k]
            if Kp!=3:
                S=sp.Integer(s)**((Kp+1)*kp) if Sexp=='derived' else 1
                f=f-q**Kp*q*graw[(t*(Kp+1))%4]*x**(Kp+1)*S*Nk[kp]
            f=sp.cancel(f/(1-q**3*x**4))
            if f==0: continue
            if sp.cancel(f/((q*x)**(l-K(k))*f.subs(x,1/(q**2*x))))!=1: ok=False
        res.append(ok)
      print("s=%+d t=%d S=%-7s l=0..8:"%(s,t,Sexp),res)
