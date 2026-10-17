# (L3) box check in log_Z units: rows m=m1 m2^2 m3^3 m4^4, N(m)<=Q=Z^q, columns c~Z.
# Bounds (normalised by ||T||^2 = Z):
#  c1 = X2X3X4 * B1(X1,Z), c2 = X1X3X4*(X2+Z), c3 = X1X2X4*B1(X3,Z), c4 = X4*(Z+P^2), c5 = X1X2X3X4*(X1^3X2^2X3)^{1/2}
#  B1(X,Z) (X<=Z, squarefree rows) = Z + Z^{3/4}X^{1/2} + min(Z^{1/3}X,(XZ)^{2/3}) + X   [improved LS via duality]
import itertools, sys
def g_new(x): return max(1.0, 0.75+x/2, min(1/3+x, 2*(x+1)/3), x)
def g_bgl(x): return max(1.0, 2*(x+1)/3, x)
def worst(q, g, use5=True, n=48):
    best=(-1,None)
    for i2 in range(n+1):
        x2=q/2*i2/n
        for i3 in range(n+1):
            x3=q/3*i3/n
            for i4 in range(n+1):
                y=q/4*i4/n
                rest=q-2*x2-3*x3-4*y
                if rest< -1e-12: continue
                for x1 in (rest, rest/2, rest/4, 0.0):
                    c=[x2+x3+y+g(x1), x1+x3+y+max(x2,1.0), x1+x2+y+g(x3), y+max(1.0,2*(x1+x2+x3))]
                    if use5: c.append(x1+x2+x3+y+(3*x1+2*x2+x3)/2)
                    v=min(c)
                    if v>best[0]: best=(v,(round(x1,4),round(x2,4),round(x3,4),round(y,4)))
    return best
q=2/3
for name,g in (("new",g_new),("bgl",g_bgl)):
    for u5 in (True,False):
        v,b=worst(q,g,u5)
        print(f"{name:4s} conv={u5!s:5s} worst={v:.5f} (target 13/12={13/12:.5f}; in Q-units {v/q:.5f}) box(x1,x2,x3,y)={b}")
