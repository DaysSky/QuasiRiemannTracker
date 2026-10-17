import sys
sys.argv=['x','0.0104166666',  '3']
exec(open('recur.py').read().split("for r in (4/3")[0])
def Tdbg(q,c,d,mu):
    q=rd(q); n=int(round(q/STEP)); rows=[]
    for i4 in range(0,n//4+1):
        for i3 in range(0,(n-4*i4)//3+1):
            for i2 in range(0,(n-4*i4-3*i3)//2+1):
                i1=n-4*i4-3*i3-2*i2
                x1,x2,x3,x4=i1*STEP,i2*STEP,i3*STEP,i4*STEP
                opts=[x2+x3+Bsf(x1,c,d), x1+x3+max(x2,c), x1+x2+Bsf(x3,c,d), max(2*(x1+x2+x3),c)]
                if mu: opts.append((x1+x2+x3)+1.5*(x1+2*x2+3*x3))
                rows.append((x4+min(opts),(round(x1,3),round(x2,3),round(x3,3),round(x4,3)),[round(o,3) for o in opts]))
    rows.sort(reverse=True); return rows[:6]
for r in (4/3,7/6):
    q=2-r
    for d in (1,2,3):
        print(r,d,"Ball",Ball(r,1.0,d,True)-r, Ball(r,1.0,d,False)-r)
        for row in Tdbg(q,1.0,d,True): print("   ",row)
