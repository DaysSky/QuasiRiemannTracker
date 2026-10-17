import sys; sys.argv=['x','all','0.1667']
src=open('hull6.py').read().split("which=sys.argv")[0]
exec(src)
# evaluation constraints with fix (a): drop x>=0, keep 3x+y>=0 (long root (3,1)), 2x+y>=0, x+y>=0, y>=0
def lpmin2(P,obj,ylo,yhi):
    n=len(P); A_eq=np.ones((1,n))
    A_ub=np.vstack([P[:,1],-P[:,1],-(3*P[:,0]+P[:,1]),-(2*P[:,0]+P[:,1]),-(P[:,0]+P[:,1])])
    r=linprog(obj,A_ub=A_ub,b_ub=[yhi,-ylo,0,0,0],A_eq=A_eq,b_eq=[1],bounds=(0,None),method="highs")
    return r.fun if r.status==0 else 9e9
def Xexp2(P,eta):
    main=1.0+lpmin2(P,P[:,0]+eta*P[:,2],0.25,0.25)
    dual=5/6+lpmin2(P,P[:,0]+(2/3)*P[:,1]+eta*P[:,2],0.0,0.2499)
    nos=5/6+lpmin2(P,P[:,0]+(2/3)*P[:,1]+eta*P[:,2],0.2501,5.0)
    return min(nos,max(main,dual)),main,dual,nos
def book2(P,H=1/3,n=25):
    etas=np.linspace(0,H,n); ex=[Xexp2(P,e)[0] for e in etas]
    best=(9,None)
    for i,k0 in enumerate(etas):
        offT=max(1/3+etas[j]+(H-etas[j])/4+2*ex[j] for j in range(i+1))
        s=max(2+(H-k0)/4,(1+offT)/2,2)
        if s<best[0]: best=(s,k0)
    return best,ex
lam=13/8
for it in range(6):
    P=sample(lam,"all"); (s,k0),ex=book2(P)
    print("fix(a): lam=%.4f X(eta=0)=%.4f X(eta=1/3)=%.4f  S<<Z^%.4f (L3)<<Q^%.4f k0=%.3f"%(lam,ex[0],ex[-1],s,1.5*(s-1),k0),flush=True)
    if lam-1.5*(s-1)<1e-4: break
    lam=1.5*(s-1)
