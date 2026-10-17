# Function field check over F_q[T], q = 5 or 13 (q = 1 mod 4).
# Verify: for monic squarefree c, the normalized quartic Gauss sum g(chi_c) = sum_{x mod c} chi_c(x) psi(x/c)
# equals (unit)^deg(c) * mu(c)-free expression involving (c'/c)_4, i.e. g(chi_c)/sqrt(q^n) = u_n * chi_c(-c') * prod ... ;
# we simply test whether g(chi_c) / (c'/c)_4 depends only on deg(c), omega(c) and chi_c(-1).
import itertools, cmath, math, random, sys
q=int(sys.argv[1]) if len(sys.argv)>1 else 5
p=q  # prime field only
# primitive root and quartic character on F_q^*
def prim_root(p):
    for g in range(2,p):
        if all(pow(g,(p-1)//r,p)!=1 for r in set(f for f in range(2,p) if (p-1)%f==0 and all(f%d for d in range(2,f)))): return g
g0=prim_root(p); I=cmath.exp(2j*math.pi/4)
dlog={pow(g0,k,p):k for k in range(p-1)}
def chi4(a):  # quartic character of F_p^*
    a%=p
    return 0 if a==0 else I**(dlog[a]%4)
# polynomials as lists of coeffs low->high
def trim(a):
    while a and a[-1]%p==0: a.pop()
    return a
def pmod(a,b):
    a=[x%p for x in a]; a=trim(a); inv=pow(b[-1],p-2,p)
    while len(a)>=len(b):
        f=a[-1]*inv%p; s=len(a)-len(b)
        for i,bc in enumerate(b): a[s+i]=(a[s+i]-f*bc)%p
        a=trim(a)
    return a
def pmul(a,b):
    if not a or not b: return []
    r=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b): r[i+j]=(r[i+j]+x*y)%p
    return trim(r)
def ppowmod(a,e,m):
    r=[1]; a=pmod(a,m)
    while e:
        if e&1: r=pmod(pmul(r,a),m)
        a=pmod(pmul(a,a),m); e>>=1
    return r
def deriv(a): return trim([(i*a[i])%p for i in range(1,len(a))])
def is_irred(f):
    n=len(f)-1
    # Rabin test via brute force for small n: no factor of degree <= n/2
    for d in range(1,n//2+1):
        # x^{q^d}-x mod f gcd
        xq=ppowmod([0,1],q**d,f); h=trim([(xq[i] if i<len(xq) else 0)-(1 if i==1 else 0) for i in range(max(len(xq),2))])
        # gcd(h,f)
        a,b=f[:],h
        while b: a,b=b,pmod(a,b)
        if len(a)>1: return False
    return True
def residues(c):
    n=len(c)-1
    for coeffs in itertools.product(range(p),repeat=n):
        yield trim(list(coeffs))
def quartic_symbol_prime(a,pi):
    # (a/pi)_4 = a^{(q^k-1)/4} mod pi, which lies in F_q; map via chi4 composed with identification x -> x (F_q subset)
    k=len(pi)-1
    r=ppowmod(a,(q**k-1)//4,pi)
    if not r: return 0
    assert len(r)==1
    v=r[0]
    # v is a 4th root of unity in F_q; identify with complex via chi4-compatible map: v = g0^{(q-1)/4 * t} -> I^t
    t=dlog[v]//((p-1)//4)
    return I**t
def factor(c):
    # trial division by monic irreducibles of small degree
    fac=[]; rem=c[:]
    n=len(c)-1
    for d in range(1,n+1):
        for coeffs in itertools.product(range(p),repeat=d):
            f=list(coeffs)+[1]
            if len(rem)-1<d: break
            if not is_irred(f): continue
            while len(rem)-1>=d and not pmod(rem,f):
                # divide
                quo=[]; a=rem[:]
                inv=1
                qcoef=[0]*(len(a)-len(f)+1)
                while len(a)>=len(f):
                    co=a[-1]; s=len(a)-len(f); qcoef[s]=co
                    for i,bc in enumerate(f): a[s+i]=(a[s+i]-co*bc)%p
                    a=trim(a)
                rem=trim(qcoef); fac.append(f)
        if len(rem)<=1: break
    return fac
def chi_c(a,c,fac):
    v=1
    for f in fac: v*=quartic_symbol_prime(a,f)
    return v
def res_inf(x,c):
    # coefficient of T^{-1} in x/c (deg x < deg c), c monic: via long division of x*T^{deg c -1}... compute Laurent: x/c = sum a_k T^{-k}
    n=len(c)-1
    # a_1 = leading coeff of x at degree n-1
    return (x[n-1] if len(x)>n-1 else 0)%p
def gauss(c,fac):
    s=0
    for x in residues(c):
        if not x: continue
        v=chi_c(x,c,fac)
        if v==0: continue
        s+=v*cmath.exp(2j*math.pi*res_inf(x,c)/p)
    return s
random.seed(1)
rows=[]
for n in (1,2,3):
    seen=0
    for coeffs in itertools.product(range(p),repeat=n):
        c=list(coeffs)+[1]
        fac=factor(c)
        if len(fac)!=len(set(map(tuple,fac))) : continue   # squarefree only
        G=gauss(c,fac)/math.sqrt(q**n)
        cp=deriv(c)
        sym=chi_c(cp,c,fac)
        mu=(-1)**len(fac)
        chim1=chi_c([p-1],c,fac)
        ratio=G/sym if sym!=0 else None
        rows.append((n,len(fac),round(chim1.real),ratio))
        seen+=1
        if seen>=40: break
# check that ratio depends only on (n, omega, chi_c(-1))
from collections import defaultdict
grp=defaultdict(list)
for n,w,cm,r in rows: grp[(n,w,cm)].append(r)
for k,v in sorted(grp.items()):
    spread=max(abs(x-v[0]) for x in v)
    print("deg %d omega %d chi(-1)=%+d : count %2d  ratio %s  spread %.2e"%(k[0],k[1],k[2],len(v),"%.4f%+.4fi"%(v[0].real,v[0].imag),spread))
