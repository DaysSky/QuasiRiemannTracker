# Exact rational functions in Z(r,a,x1,x2) on top of FLINT multivariate polynomials, with the G2 monomial substitutions.
import flint
ctx=flint.fmpz_mpoly_ctx.get(('r','a','x1','x2'),'lex')
R,A,X1,X2=ctx.gens()
ONE=ctx.from_dict({(0,0,0,0):1}); ZERO=ctx.from_dict({})
Cm={(1,1):2,(2,2):2,(1,2):-1,(2,1):-3}
class F:
    __slots__=('n','d')
    def __init__(s,n,d=None,norm=True):
        if d is None: d=ONE
        if isinstance(n,int): n=ctx.from_dict({(0,0,0,0):n}) if n else ctx.from_dict({})
        if isinstance(d,int): d=ctx.from_dict({(0,0,0,0):d})
        if norm and not n.is_zero():
            g=n.gcd(d)
            if not g.is_one(): n=n/g; d=d/g
            if d.leading_coefficient()<0: n=-n; d=-d
        if n.is_zero(): d=ONE
        s.n=n; s.d=d
    def __add__(s,o):
        o=o if isinstance(o,F) else F(o)
        if s.d==o.d: return F(s.n+o.n,s.d)
        return F(s.n*o.d+o.n*s.d,s.d*o.d)
    __radd__=__add__
    def __neg__(s): return F(-s.n,s.d,False)
    def __sub__(s,o): return s+(-(o if isinstance(o,F) else F(o)))
    def __rsub__(s,o): return F(o)+(-s)
    def __mul__(s,o):
        o=o if isinstance(o,F) else F(o)
        return F(s.n*o.n,s.d*o.d)
    __rmul__=__mul__
    def inv(s): return F(s.d,s.n)
    def __truediv__(s,o):
        o=o if isinstance(o,F) else F(o)
        return s*o.inv()
    def __pow__(s,e):
        if e>=0: return F(s.n**e,s.d**e,False)
        return F(s.d**(-e),s.n**(-e),False)
    def __eq__(s,o):
        o=o if isinstance(o,F) else F(o)
        return (s.n*o.d-o.n*s.d).is_zero()
    def iszero(s): return s.n.is_zero()
def mono_sub(P,j):
    # P(sigma_j x) as (poly G, monomial M) with P(sigma_j x) = G/M ; sigma_j: x_i -> q^{-c(i,j)} x_i x_j^{-c(i,j)}, q=r^2
    terms=[]
    for e,c in P.to_dict().items():
        er,ea,e1,e2=e
        sh=-(Cm[(1,j)]*e1+Cm[(2,j)]*e2)
        ne=[er+2*sh,ea,e1,e2]; ne[1+j]+=sh
        terms.append((ne,c))
    mins=[min(0,min(t[0][k] for t in terms)) for k in range(4)] if terms else [0]*4
    G=ctx.from_dict({tuple(t[0][k]-mins[k] for k in range(4)):t[1] for t in terms}) if terms else ctx.from_dict({})
    M=ctx.from_dict({tuple(-m for m in mins):1})
    return G,M
def sub(f,j):
    Gn,Mn=mono_sub(f.n,j); Gd,Md=mono_sub(f.d,j)
    return F(Gn*Md,Gd*Mn)
