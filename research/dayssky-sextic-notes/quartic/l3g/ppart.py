# CG p-part N(x; lam) for G2, n=4, general twisting parameter lam=(l1,l2), with gamma table supplied.
# Node 1 short (c, |a|^2=1), node 2 long (d, |a|^2=3).  Same conventions as l3e/cg_ppart2.py.
import sympy as sp
x1,x2=sp.symbols('x1 x2'); X=(x1,x2)
n=4; L2={1:1,2:3}; C={(1,1):2,(2,2):2,(1,2):-1,(2,1):-3}; m={1:4,2:4}
pos=[(1,0),(0,1),(1,1),(2,1),(3,1),(3,2)]
words=[[]]
for k in range(1,7):
    for st in (1,2):
        if k==6 and st==2: continue
        words.append([st if j%2==0 else 3-st for j in range(k)])
def ppart(q,gam,l,sg=None):
    if sg is None: sg=lambda b,i,mu: 1
    pair=lambda b,i: b[0]*C[(1,i)]+b[1]*C[(2,i)]
    def dot(b,i): b=list(b); b[i-1]+=-pair(b,i)+(l[i-1]+1); return tuple(b)
    cls=lambda b:(b[0]%n,b[1]%n)
    mu=lambda b,i: l[i-1]+1-pair(b,i)
    P=lambda b,i,x: (q*x)**(l[i-1]+1-(mu(b,i)%m[i]))*(1-1/q)/(1-(q*x)**m[i]/q)
    Q=lambda b,i,x: -sg(b,i,mu(b,i))*gam[(-L2[i]*mu(b,i))%n]*(q*x)**(l[i-1]+1-m[i])*(1-(q*x)**m[i])/(1-(q*x)**m[i]/q)
    def sx(j,x): return tuple(q**(-C[(i,j)])*x[i-1]*x[j-1]**(-C[(i,j)]) for i in (1,2))
    def act(elem,i):
        out={}
        S=sx(i,X)
        for lam,f in elem.items():
            g=f.subs({x1:S[0],x2:S[1]},simultaneous=True)
            lam2=cls(dot(lam,i))
            out[lam]=out.get(lam,0)+P(lam,i,X[i-1])*g
            out[lam2]=out.get(lam2,0)+Q(lam2,i,X[i-1])*g
        return out
    def wx(word):
        y=X
        for i in reversed(word): y=sx(i,y)
        return y
    Delta=lambda y: sp.Mul(*[(1-q**(4*(k1+k2))*y[0]**(4*k1)*y[1]**(4*k2)) for k1,k2 in pos])
    D=sp.Mul(*[(1-q**(4*(k1+k2)-1)*x1**(4*k1)*x2**(4*k2)) for k1,k2 in pos])
    h=0
    for w in words:
        e={(0,0):sp.Integer(1)}
        for i in w: e=act(e,i)
        h+=sum(e.values())/Delta(wx(w))
    N=sp.cancel(sp.together(h*D))
    num,den=sp.fraction(N)
    return sp.Poly(sp.expand(num),x1,x2), sp.expand(den)
