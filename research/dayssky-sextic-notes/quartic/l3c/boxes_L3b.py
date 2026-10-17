# (L3) boxes, log_Z units. Rows m = m1 m2^2 m3^3 m4^4 (N(m)<=Q=Z^q), columns c ~ Z (T-family, ||T||^2 = Z).
# Per-box bounds for Delta(box) = sup_{lambda} sum_{m in box}|sum_c lambda_c chi_c(m)^3|^2/||lambda||^2  (normalised so diag = Z -> 1).
import math
q=2/3
def B2(k,n):   # all rows ~K=Z^k, sqfree quartic columns ~N=Z^n: min(Lemma A, Poisson-improved)
    a=max(k, k/4+n, 2*(k+n)/3)
    b=max(k, n, 0.75*k+0.5*n, k/3+n)
    return min(a,b)
def B1(x,n):   # sqfree rows ~X, columns ~N, improved + duality + BGL
    def one(M,N): return max(M,N,0.75*M+0.5*N, min(M/3+N, 2*(M+N)/3))
    return min(one(x,n), one(n,x))
def B2q(k,n): return max(k,n,k/2+n)   # quadratic, all rows
def bounds(x1,x2,x3,y,use=('row','cls','conv','pois')):
    P=x1+x2+x3; c={}
    if 'row' in use:
        c['c1']=x2+x3+y+B1(x1,1.0); c['c2']=x1+x3+y+max(x2,1.0); c['c3']=x1+x2+y+B1(x3,1.0)
    if 'cls' in use: c['c4']=y+max(1.0,2*P)
    if 'conv' in use: c['c5']=P+y+(3*x1+2*x2+x3)/2
    if 'pois' in use:
        h=2*P-1
        if h<0: c['p0']=y+1.0
        else:
            c['p1']=y+max(1.0,1-P+x2+x3+B2(h,x1))
            c['p2']=y+max(1.0,1-P+x1+x3+B2q(h,x2))
            c['p3']=y+max(1.0,1-P+x1+x2+B2(h,x3))
            c["p4"]=y+max(1.0,1-P+max(2*h,h/4+x1+2*x2+3*x3))
    return c
def worst(use,n=40):
    best=(-1,None,None)
    for i1 in range(n+1):
      x1=q*i1/n
      for i2 in range(n+1):
        x2=q/2*i2/n
        if x1+2*x2>q+1e-12: break
        for i3 in range(n+1):
          x3=q/3*i3/n
          if x1+2*x2+3*x3>q+1e-12: break
          y=(q-x1-2*x2-3*x3)/4
          for yy in (y,0.0,y/2):
            c=bounds(x1,x2,x3,yy,use); k=min(c,key=c.get); v=c[k]
            if v>best[0]: best=(v,(round(x1,4),round(x2,4),round(x3,4),round(yy,4)),k)
    return best
for use in [('row','cls','conv'),('row','cls','conv','pois'),('row','cls','pois'),('cls','conv','pois')]:
    v,b,k=worst(use)
    print(f"{'+'.join(use):22s} worst={v:.5f} Q-exp={v/q:.5f} at {b} (binding {k}); target Q^1.625 -> {1.625*q:.5f}")
