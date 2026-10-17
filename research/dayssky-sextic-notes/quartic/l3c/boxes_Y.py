# Q-units. Rows m=m1 m2^2 m3^3 m4^4, N(m)<=Q; columns c~Y=Q^eta (T-family). Bound Delta(Q,Y)/1 in Q-exponent; target 13/8.
def B2(k,n):
    return min(max(k, k/4+n, 2*(k+n)/3), max(k, n, 0.75*k+0.5*n, k/3+n))
def bounds(x1,x2,x3,y,eta):
    P=x1+x2+x3; c={}
    c['cls']=y+max(eta,2*P)
    c['conv']=P+y+(3*x1+2*x2+x3)/2
    h=2*P-eta
    if h<0: c['p0']=y+eta
    else: c['p1']=y+max(eta, eta-P+x2+x3+B2(h,x1))
    return c
def worst(eta,n=36):
    best=(-1,None,None)
    for i1 in range(n+1):
      x1=i1/n
      for i2 in range(n+1):
        x2=0.5*i2/n
        if x1+2*x2>1+1e-12: break
        for i3 in range(n+1):
          x3=i3/(3*n)
          if x1+2*x2+3*x3>1+1e-12: break
          ymax=(1-x1-2*x2-3*x3)/4
          for yy in (ymax,0.0,ymax/2,ymax/4,3*ymax/4):
            c=bounds(x1,x2,x3,yy,eta); k=min(c,key=c.get); v=c[k]
            if v>best[0]: best=(v,(round(x1,3),round(x2,3),round(x3,3),round(yy,3)),k)
    return best
for eta in [0.25,0.5,0.75,1.0,1.1,1.2,1.25,1.3,1.4,1.45,1.5]:
    v,b,k=worst(eta); print(f"eta={eta:5.2f} worst={v:.5f} {'OK' if v<=13/8+1e-9 else 'FAIL'} at {b} via {k}")
