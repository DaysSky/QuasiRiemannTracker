# q-exponents of the slice-signed CG p-part coefficients for twisting parameter l, and the worst normalised monomial size at a point
# (x', y') in chamber coordinates: e(k) - k1(1+x') - k2(1+y').  Exponent e(k) measured as log_r |coef| at r = 10007, 10009
# (gamma(1) = 1/r, gamma(2) = 1/r, gamma(3) = s/(q gamma(1))), i.e. Gauss sums of exact size sqrt(q).
import sys, math
from fractions import Fraction as Fr
from ratf import F, ctx
import thm41sym as T
def exps(s,l):
    N,_=T.ppart(s,l,'slice')
    out=None
    for r0 in (10007,10009):
        res={}
        num=N.n.to_dict(); den=N.d.to_dict()
        # evaluate coefficient polys in (r,a) at r=r0, a=1/r0
        def ev(d):
            acc={}
            for e,c in d.items():
                key=(int(e[2]),int(e[3]))
                acc[key]=acc.get(key,Fr(0))+Fr(int(c))*Fr(r0)**int(e[0])*Fr(1,r0)**int(e[1])
            return acc
        nv=ev(num); dv=ev(den); dval=dv[(0,0)]
        for k,v in nv.items():
            if v!=0:
                w=abs(v/dval); res[k]=(math.log(w.numerator)-math.log(w.denominator))/math.log(r0)/2   # q = r^2: exponent in q
        out=res if out is None else {k:(out[k],res[k]) for k in res}
    return out
s=int(sys.argv[1]); pts=[(Fr(1,3),Fr(0)),(Fr(1,2),Fr(0)),(Fr(1,3),Fr(1,4))]
for l in [tuple(map(int,z.split(','))) for z in sys.argv[2].split(';')]:
    E=exps(s,l)
    ee={k:round((a+b)/2*2)/2 for k,(a,b) in E.items()}   # half-integers
    worst={}; worst2={}
    for (xp,yp) in pts:
        w=max(((ee[k]-k[0]*(1+xp)-k[1]*(1+yp)),k) for k in ee if k not in ((0,0),))
        worst[(float(xp),float(yp))]=(round(float(w[0]),3),w[1])
        cand=[((ee[k]-k[0]*(1+xp)-k[1]*(1+yp)),k) for k in ee if k not in ((0,0),(1,0),(0,1))]
        if cand:
            w2=max(cand); worst2[(round(float(xp),3),float(yp))]=(round(float(w2[0]),3),w2[1])
    nontriv=sorted(k for k in ee if k!=(0,0))
    print("l=%s: %d mon.; exps consistent=%s\n   worst incl. first order: %s\n   worst excl. first order: %s"%(l,len(ee),all(abs(a-b)<1e-3 for a,b in E.values()),worst,worst2),flush=True)
