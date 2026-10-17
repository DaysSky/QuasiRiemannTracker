# Mean-value Bochner hull for Z(s,w) = sum_c T(c) N(c)^{-s} psi(h0 c, w), coordinates x = Re s - 1/2, y = Re w - 1/2.
# Regions where the series converge (mean-value inputs), with h0-exponent b:
#  R_d : x > 1/2 + max(0,-y)          [LoA for psi(h0 c, .) over c; b = 1/6 + max(0,-y)]
#  R_c : y > g(x), g = lam/2 - 3x (x<=0), lam/2 (0<x<=1/2), 1/2 (x>1/2)   [(L3)-type mean square, exponent lam; b = 0]
#  Sd(R_c) : (x+y, -y) in R_c         [d-FE then R_c; b = -y]
#  Sc(R_d) : (-x, y+3x) in R_d        [c-FE then R_d for the l=3 family; b = 1/6 + max(0,-(y+3x))]
import numpy as np, itertools, sys
lam=float(sys.argv[1]) if len(sys.argv)>1 else 13/8
use=sys.argv[2] if len(sys.argv)>2 else "dc"
def g(x): return lam/2-3*x if x<=0 else (lam/2 if x<=0.5 else 0.5)
def h(y): return 0.5+max(0,-y)
pts=[]  # sample boundary points (x,y,b) of each region (slightly inside)
e=1e-6
for y in np.linspace(-6,6,2401):
    pts.append((h(y)+e,y,1/6+max(0,-y),'d'))
for x in np.linspace(-3,3,2401):
    pts.append((x,g(x)+e,0.0,'c'))
if 'D' in use:   # Sd(R_c): points (x,y) with -y = g(x+y)+e  -> param u = x+y, y = -(g(u)+e), x = u - y
    for u in np.linspace(-3,3,2401):
        y=-(g(u)+e); x=u-y; pts.append((x,y,-y,'Dc'))
if 'C' in use:   # Sc(R_d): (-x, y+3x) = (h(v)+e, v) -> x = -(h(v)+e), y = v - 3x
    for v in np.linspace(-6,6,2401):
        x=-(h(v)+e); y=v-3*x; pts.append((x,y,1/6+max(0,-v),'Cd'))
P=np.array([(p[0],p[1],p[2]) for p in pts])
# residue line y=1/4: minimal x over convex combinations of two sample points (hull in 2D: two points suffice for a line section)
def line_min(y0, weight_b=0.0):
    best=(9,None)
    A=P[P[:,1]<=y0]; B=P[P[:,1]>=y0]
    # subsample for speed
    A=A[::3]; B=B[::3]
    for a in A:
        t=(B[:,1]-y0)/(B[:,1]-a[1]+1e-300)   # weight on a
        x=t*a[0]+(1-t)*B[:,0]; b=t*a[2]+(1-t)*B[:,2]
        val=x+weight_b*b
        i=np.argmin(val)
        if val[i]<best[0]: best=(val[i],(x[i],b[i],a,B[i]))
    return best
r=line_min(0.25)
print("lam=%.4f use=%s : residue line y=1/4 reachable for x > %.4f  (b=%.3f)  pi = %.4f"%(lam,use,r[1][0],r[1][1],0.5+r[1][0]))
print("   from", np.round(r[1][2],4), "and", np.round(r[1][3],4))
