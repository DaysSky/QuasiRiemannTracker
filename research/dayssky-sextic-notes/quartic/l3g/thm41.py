# Local test of CG Theorem 4.1 for the TRUE slice factor over Q(i) (mu_8 not in F):
#   f'_k(x) = [N^(k)(x) - [K(k')!=3] G(K') x^(K'+1) S(k') N^(k')(x)] / (1 - q^3 x^4),
#   k' = (l+1-k) mod 4, K(k) = (l-2k) mod 4, G(K) = q^K * q*gamma_raw(t(K+1)), S(k') = s^((K'+1)(k'+sum_{j>i} c(j,i) beta_j))  [CG (5.7) twists Psi_i only by j>i],
# where l = lam_i + sum_{j!=i} |c(j,i)| beta_j.  Required (global FE via Kubota): f'_k(x) = (qx)^(l-K(k)) f'_k(1/(q^2 x)).
import sympy as sp, sys, itertools
from ppart import ppart, x1, x2, C, L2
s=int(sys.argv[1]); mode=sys.argv[2]   # mode: raw | bd
r=sp.Integer(3); q=r*r; a=sp.Rational(2,7)
graw={0:sp.Integer(-1),1:a,2:1/r,3:sp.Integer(s)/(q*a)}
gbd={i:graw[i]*sp.Integer(s)**(i*(i-1)//2) for i in graw}
gam=graw if mode=='raw' else gbd
useS = (len(sys.argv)<=3 or sys.argv[3]!='noS')
x=sp.symbols('x')
tot=0; bad=0; ratios={}
for lam in itertools.product(range(2),range(2)):
    Np,den=ppart(q,gam,lam)
    assert den==1, den
    coef=Np.as_dict()
    for i in (1,2):
        j=3-i; t=L2[i]
        betas=sorted(set(e[j-1] for e in coef))
        for bj in betas:
            l=lam[i-1]+(-C[(j,i)])*bj
            Nk={k:sum(v*x**e[i-1] for e,v in coef.items() if e[j-1]==bj and e[i-1]%4==k) for k in range(4)}
            K=lambda k:(l-2*k)%4
            for k in range(4):
                kp=(l+1-k)%4; Kp=K(kp)
                f=Nk[k]
                if Kp!=3:
                    S=sp.Integer(s)**((Kp+1)*(kp+(C[(j,i)]*bj if j>i else 0))) if useS else 1
                    f=f-q**Kp*q*graw[(t*(Kp+1))%4]*x**(Kp+1)*S*Nk[kp]
                f=sp.cancel(f/(1-q**3*x**4))
                if f==0: continue
                tot+=1
                rat=sp.cancel(f/((q*x)**(l-K(k))*f.subs(x,1/(q**2*x))))
                if rat!=1:
                    bad+=1; ratios.setdefault(str(rat),[]).append((lam,i,bj,k))
print("s=%+d mode=%s S=%s: %d slice classes tested, %d fail"%(s,mode,useS,tot,bad))
for rr,v in list(ratios.items())[:8]: print("  ratio",rr,"at",v[:6],"... (%d)"%len(v))
