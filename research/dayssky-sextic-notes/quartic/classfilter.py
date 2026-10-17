# Test: restrict columns n (squarefree primary, N n ~ X) to a class of N(n) mod 8 and/or n mod 4.
# Heuristic (Poisson in row variable + Gauss-sum squares are Hecke): for a single class with
# chi_n(-1)=1 and trivial reciprocity sign, the dual quadratic form is positive and
#   sum_{N h<=H} |sum_n alpha(n)^{-1} g~(n) chi_n(h)|^2 ~ H*X + X^2  (i.e. ratio ~ 1 + X/H),
# while for the full family the sign (-1)^{ab} can cancel the X^2 term.
import numpy as np, math, sys, time
from dualms import build
from gauss4 import chi_prime, gmul

def run(X,H,classes,nrand=2,seed=1):
    hs,hN,P,tb=build(X,H)
    norms=[p[1] for p in P]
    order=sorted(range(len(P)),key=lambda j:norms[j])
    ns=[]
    def rec(start,prod_norm,lst):
        if prod_norm>X/2 and lst: ns.append(list(lst))
        for t in range(start,len(order)):
            j=order[t]
            if prod_norm*norms[j]>X: break
            lst.append(j); rec(t+1,prod_norm*norms[j],lst); lst.pop()
    rec(0,1,[])
    cache={}
    def cpp(i,j):
        if (i,j) not in cache: cache[(i,j)]=chi_prime(P[j][0],P[i][0])
        return cache[(i,j)]
    nH=len(hs)
    acc={c:np.zeros(nH,complex) for c in classes}
    accr={c:np.zeros(nH,complex) for c in classes}
    w2={c:0.0 for c in classes}
    rng=np.random.default_rng(seed)
    for lst in ns:
        n=(1,0); gam=1.0+0j; Nn=1
        for t,j in enumerate(lst):
            gam*=P[j][2]
            for s_ in lst[:t]:
                gam*=cpp(s_,j)*cpp(j,s_)
            n=gmul(n,P[j][0]); Nn*=P[j][1]
        tt=Nn/X
        W=math.sin(math.pi*(2*tt-1))**2 if 0.5<tt<=1 else 0.0
        if W==0: continue
        alpha=complex(*n)/math.sqrt(Nn)
        ex=np.zeros(nH,np.int16); zero=np.zeros(nH,bool)
        for j in lst:
            e=P[j][3]; zero|=(e<0); ex+=np.where(e<0,0,e)
        v=np.where(zero,0,1j**(ex%4))
        key8=Nn%8
        cls4=1 if (n[0]%4,n[1]%4)==(1,0) else 3
        for c in classes:
            if c=='all' or (c=='N1' and key8==1) or (c=='N5' and key8==5) or (c=='n1' and cls4==1) or (c=='n3' and cls4==3):
                acc[c]+=W*np.conj(alpha)*gam*v
                accr[c]+=W*np.exp(2j*np.pi*rng.random())*v
                w2[c]+=W*W
    return hs,hN,acc,accr,w2

if __name__=="__main__":
    X=int(sys.argv[1]); H=int(sys.argv[2])
    t=time.time()
    classes=['all','N1','N5','n1','n3']
    hs,hN,acc,accr,w2=run(X,H,classes)
    print(f"X={X} H={H} rows={len(hs)} time={time.time()-t:.1f}s")
    for Hc in [H//16,H//4,H]:
        m=hN<=Hc; nh=int(m.sum())
        line=f" H'={Hc:6d} #h={nh:6d} X/H'={X/Hc:7.1f} |"
        for c in classes:
            line+=f" {c}:{np.sum(np.abs(acc[c][m])**2)/(nh*w2[c]):7.3f}(r{np.sum(np.abs(accr[c][m])**2)/(nh*w2[c]):.2f})"
        print(line)
