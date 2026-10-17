exec(open('boxes_L3b.py').read().split("def worst")[0])
use=('row','cls','conv','pois')
res=[]
n=60
for i1 in range(n+1):
  x1=q*i1/n
  for i2 in range(n+1):
    x2=q/2*i2/n
    if x1+2*x2>q+1e-12: break
    for i3 in range(n+1):
      x3=q/3*i3/n
      if x1+2*x2+3*x3>q+1e-12: break
      ymax=(q-x1-2*x2-3*x3)/4
      for yy in (ymax,0.0,ymax/2):
        c=bounds(x1,x2,x3,yy,use); k=min(c,key=c.get); res.append((c[k],(round(x1,3),round(x2,3),round(x3,3),round(yy,3)),k))
res.sort(reverse=True)
print("top boxes (Z-exponent of S, (L3) Q-exponent = 1.5*(v-1)... v in log_Z of S/||T||^2 * ... ):")
seen=0
for v,b,k in res:
    if b[0]>=q-1e-9 and b[1]==0 and b[2]==0: tag="MAIN"
    else: tag=""
    print("  %.5f  Q-exp(L3)=%.4f  box=%s binding=%s %s"%(v,1.5*v, b,k,tag)); seen+=1
    if seen>12: break
# best non-main
nm=[r for r in res if not (r[1][0]>=q-0.02)]
print("worst box with x1 <= q-0.02:",nm[0], " -> (L3) Q-exp %.4f"%(1.5*nm[0][0]))
