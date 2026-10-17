# EXACT (rational) certificate of the cone lemma:
#   for every block (rho,tC,tD) with tK=1-rho, tE=rho, 0<=tD<=(1+tC)/2  (theta=1/2; y=0 so dual blocks coincide),
#   min(B1,B2) <= h := c1*Delta + rho/2,   at x = 1/3 + dp, for any table S <= L + Delta*1[q' in [2/3,5/6)] (conservative sigma).
# f is piecewise linear; on each linearity cell the max of min(B1,B2) is at a vertex of the cell split by B1=B2.
import itertools, sys
from fractions import Fraction as F
def run(dp,Delta,c1=F(3),M=F(10**5)):
    x=F(1,3)+dp
    # linear forms a0 + a1*rho + a2*tC + a3*tD as tuples
    def L(a0=0,ar=0,ac=0,ad=0): return (F(a0),F(ar),F(ac),F(ad))
    def add(*fs): return tuple(sum(c) for c in zip(*fs))
    def sc(k,f): return tuple(k*c for c in f)
    tK=L(1,-1); tC=L(0,0,1); tD=L(0,0,0,1); one=L(1)
    P1=[add(tK,tC), tD, sc(F(2,3),add(tK,tC,tD))]
    P2=[tK, tD, sc(F(2,3),add(tK,tD))]
    sig=[ (add(sc(F(2,3),tC),sc(2,tD)),            [add(sc(3,tD),sc(-1,tC)), add(sc(5,tC),sc(-9,tD))]),          # 1/3<=q<=5/9
          (add(sc(F(2,3),tC),sc(2,tD),sc(Delta,add(sc(3,tD),sc(-1,tC)))), [add(sc(9,tD),sc(-5,tC)), add(sc(2,tC),sc(-3,tD))]),  # 5/9..2/3
          (add(sc(F(4,3)+Delta,tC),tD),               [add(sc(3,tD),sc(-2,tC)), add(sc(5,tC),sc(-6,tD))]),          # 2/3..5/6
          (add(sc(F(4,3),tC),tD),                     [add(sc(6,tD),sc(-5,tC))]) ]                                   # q>=5/6 (incl. tC=0)
    h=add(sc(c1*Delta,one),L(0,F(1,2)))
    dom=[L(0,1), L(1,-1), tC, add(sc(M,one),sc(-1,tC)), tD, add(L(F(1,2)),sc(F(1,2),tC),sc(-1,tD))]   # each >= 0
    def ev(f,v): return f[0]+f[1]*v[0]+f[2]*v[1]+f[3]*v[2]
    def solve3(rows):
        A=[list(r[1:]) for r in rows]; b=[-r[0] for r in rows]
        # Cramer / Gaussian elimination with Fractions
        m=[A[i]+[b[i]] for i in range(3)]
        for col in range(3):
            piv=next((r for r in range(col,3) if m[r][col]!=0),None)
            if piv is None: return None
            m[col],m[piv]=m[piv],m[col]
            for r in range(3):
                if r!=col and m[r][col]!=0:
                    k=m[r][col]/m[col][col]; m[r]=[a-k*c for a,c in zip(m[r],m[col])]
        return tuple(m[i][3]/m[i][i] for i in range(3))
    best=(None,None); nv=0
    for (sg,secs) in sig:
        for i in range(3):
            for j in range(3):
                B1=add(sc(-x,tC),sc(F(1,2),add(P1[i],sc(-1,tK))))
                B2=add(sc(-(F(1,2)+x),tC),sc(F(1,2),sg),sc(F(1,2),add(P2[j],sc(-1,tK),sc(-1,tD))))
                cons=dom+secs+[add(P1[i],sc(-1,P1[k])) for k in range(3) if k!=i]+[add(P2[j],sc(-1,P2[k])) for k in range(3) if k!=j]
                for side in (0,1):
                    extra=add(B2,sc(-1,B1)) if side==0 else add(B1,sc(-1,B2))   # side0: B1<=B2 -> f=B1
                    cc=cons+[extra]; obj=B1 if side==0 else B2
                    for tri in itertools.combinations(cc,3):
                        v=solve3(tri)
                        if v is None: continue
                        if all(ev(c,v)>=0 for c in cc):
                            nv+=1; val=ev(obj,v)-ev(h,v)
                            if best[0] is None or val>best[0]: best=(val,v)
    return best,nv
for dp in (F(1,100),F(1,1000)):
    for Delta in (F(0),F(1,1000),F(1,50),F(1,12)):
        (val,v),nv=run(dp,Delta)
        print("dp=%s Delta=%s: max(min(B1,B2) - (3*Delta+rho/2)) = %s (=%.3e) at (rho,tC,tD)=%s  [%d feasible vertices]"%(dp,Delta,val,float(val),tuple(str(a) for a in v),nv),flush=True)
