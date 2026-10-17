# Function-field (L3) test: R(m) = sum_{c monic sqfree deg n} T(c) conj(chi_c(m)),  T(c) = mu(c) conj g~(chibar_c),
# compared with the explicit formula R(m) = v_n sum_c chi_c(-1) (c'/c)_4 conj chi_c(m).
import sys, itertools, cmath, math
q=int(sys.argv[1]); n=int(sys.argv[2]); d=int(sys.argv[3])
sys.argv=['x',str(q)]
src=open('ffgauss.py').read().split("random.seed(1)")[0]
exec(src)
cs=[]
for coeffs in itertools.product(range(p),repeat=n):
    c=list(coeffs)+[1]; fac=factor(c)
    if len(fac)!=len(set(map(tuple,fac))): continue
    cs.append((c,fac))
ms=[list(co)+[1] for co in itertools.product(range(p),repeat=d)]
print("q=%d n=%d d=%d: %d sqfree c, %d rows m"%(q,n,d,len(cs),len(ms)))
T=[]; Tf=[]
for c,fac in cs:
    # chibar Gauss sum: sum conj chi_c(x) psi(x/c)
    s=0
    for x in residues(c):
        if not x: continue
        v=chi_c(x,c,fac)
        if v==0: continue
        s+=v.conjugate()*cmath.exp(2j*math.pi*res_inf(x,c)/p)
    gt=s/math.sqrt(q**n)
    mu=(-1)**len(fac)
    T.append(mu*gt.conjugate())
    Tf.append(chi_c([p-1],c,fac)*chi_c(deriv(c),c,fac))
# v_n from ratio
ratios=set(complex(round((a/b).real,6),round((a/b).imag,6)) for a,b in zip(T,Tf))
print("T(c)/(chi_c(-1)(c'/c)_4) takes values:",ratios)
S=0; diag=0; rowvals=[]
for m in ms:
    R=0; dg=0
    for (c,fac),t in zip(cs,T):
        v=chi_c(m,c,fac)
        if v==0: continue
        R+=t*v.conjugate(); dg+=1
    S+=abs(R)**2; diag+=dg; rowvals.append(abs(R)**2/dg if dg else 0)
print("S/diag = %.4f   row |R|^2/diag: min %.3f max %.3f"%(S/diag,min(rowvals),max(rowvals)))
