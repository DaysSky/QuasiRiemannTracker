# Compare the slice-signed p-part with the BD-signed one (S51 §19) monomial by monomial, lam=(0,0).
import sympy as sp
from ppart import ppart, L2, C
s=-1; r=sp.Integer(3); q=r*r; a=sp.Rational(2,7)
graw={0:sp.Integer(-1),1:a,2:1/r,3:sp.Integer(s)/(q*a)}
gbd={i:graw[i]*sp.Integer(s)**(i*(i-1)//2) for i in graw}
def B(b,i): return b[i-1]+sum(C[(j,i)]*b[j-1] for j in (1,2) if j>i)
sl=lambda b,i,mu: s**((mu*(B(b,i)+1))%2)
for lam in [(0,0),(1,0),(0,1),(1,1)]:
    Ns,_=ppart(q,graw,lam,sl); Nb,_=ppart(q,gbd,lam)
    ds,db=Ns.as_dict(),Nb.as_dict()
    keys=sorted(set(ds)|set(db))
    rat={k:(sp.nsimplify(ds.get(k,0)/db[k]) if db.get(k,0)!=0 else ('only slice' if ds.get(k,0)!=0 else 0)) for k in keys}
    print(lam,"support equal:",set(ds)==set(db)," ratios slice/BD:",rat)
