# CG p-part N(x;l) = h D, h = sum_w (1|_l w)(x)/Delta(wx), for G2, n=4, l=(0,0); s = chi_pi(-1) with the binomial fix of S51 §19.
import sympy as sp, sys
from fractions import Fraction as Fr
x1,x2=sp.symbols('x1 x2'); X=(x1,x2)
n=4; L2={1:1,2:3}; C={(1,1):2,(2,2):2,(1,2):-1,(2,1):-3}; m={1:4,2:4}; l=(0,0)
s=int(sys.argv[1]) if len(sys.argv)>1 else 1
q=sp.Integer(25); r=sp.Integer(5); a=sp.Rational(2,7)
gam={0:-1,1:a,2:1/r,3:sp.Integer(s)/(q*a)}
gam={i:gam[i]*sp.Integer(s)**(i*(i-1)//2) for i in gam}
pair=lambda b,i: b[0]*C[(1,i)]+b[1]*C[(2,i)]
def dot(b,i): b=list(b); b[i-1]+=-pair(b,i)+(l[i-1]+1); return tuple(b)
cls=lambda b:(b[0]%n,b[1]%n)
mu=lambda b,i: l[i-1]+1-pair(b,i)
P=lambda b,i,x: (q*x)**(l[i-1]+1-(mu(b,i)%m[i]))*(1-1/q)/(1-(q*x)**m[i]/q)
Q=lambda b,i,x: -gam[(-L2[i]*mu(b,i))%n]*(q*x)**(l[i-1]+1-m[i])*(1-(q*x)**m[i])/(1-(q*x)**m[i]/q)
def sx(j,x): return tuple(q**(-C[(i,j)])*x[i-1]*x[j-1]**(-C[(i,j)]) for i in (1,2))
def act(elem,i):   # elem: dict class -> expr
    out={}
    for lam,f in elem.items():
        g=f.subs({x1:sx(i,X)[0],x2:sx(i,X)[1]},simultaneous=True)
        lam2=cls(dot(lam,i))
        out[lam]=out.get(lam,0)+P(lam,i,X[i-1])*g
        out[lam2]=out.get(lam2,0)+Q(lam2,i,X[i-1])*g
    return out
# words for the 12 elements (right action: 1|w for w = s_{i1}...s_{ik} means apply s_{i1} first)
words=[[]]
for k in range(1,7):
    for st in (1,2):
        w=[st if j%2==0 else 3-st for j in range(k)]
        if k==6 and st==2: continue
        words.append(w)
# positive roots of G2 as (k1,k2) in simple roots: short a1, long a2
pos=[(1,0),(0,1),(1,1),(2,1),(3,1),(3,2)]
def wx(word):  # point w x for the left action composite used in Delta(wx): apply sx along word in order
    y=X
    for i in reversed(word): y=sx(i,y)
    return y
def Delta(y): 
    return sp.Mul(*[(1-q**(4*(k1+k2))*y[0]**(4*k1)*y[1]**(4*k2)) for k1,k2 in pos])
D=sp.Mul(*[(1-q**(4*(k1+k2)-1)*x1**(4*k1)*x2**(4*k2)) for k1,k2 in pos])
h=0
for w in words:
    e={(0,0):sp.Integer(1)}
    for i in w: e=act(e,i)
    one_w=sum(e.values())
    h+=one_w/Delta(wx(w))
N=sp.factor(sp.cancel(sp.together(h*D)))
Np=sp.Poly(sp.expand(sp.numer(sp.together(N))),x1,x2); den=sp.denom(sp.together(N))
print("denominator:",den)
print("support (k1,k2):",sorted(Np.as_dict().keys()))
print("axis check (3.26), x2=0:",sp.simplify(sp.series((N/D).subs(x2,0),x1,0,6).removeO()))
