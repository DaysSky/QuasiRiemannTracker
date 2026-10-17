import numpy as np, json, sys
G={}; exec(open('recur9.py').read().split('if __name__')[0],G)
globals().update(G)
d=json.load(open("recur9_tab.json")); tab=np.array(d["tab"]); lam=d["lam"]
loc=local_points(lam,tab); bc={}
def tagged(tK,tE):
    pts=[];tags=[]
    lb={pol:[local_b(it,tK,tE,tab,pol,bc) for it in loc[pol]] for pol in (False,True)}
    for k,(M,c) in ch.items():
        pol=POLAR[k]; Mi=np.linalg.inv(M); tw=tuple((M@np.array([1,0])).astype(int))
        for it,b in zip(loc[pol],lb[pol]):
            v=Mi@np.array(it[0]); pts.append((v[0],v[1],b+(c@v)*(tK+tE))); tags.append((tw,it[1],pol,tuple(np.round(it[0],3)),b))
    return np.array(pts),tags
def lp(P,obj,ylo,yhi):
    n=len(P); A_ub=np.vstack([P[:,1],-P[:,1],-(3*P[:,0]+P[:,1]),-(2*P[:,0]+P[:,1]),-(P[:,0]+P[:,1])])
    return linprog(obj,A_ub=A_ub,b_ub=[yhi,-ylo,0,0,0],A_eq=np.ones((1,n)),b_eq=[1],bounds=(0,None),method="highs")
q=2/3
for eta in [float(a) for a in sys.argv[1:]] or (0.0,0.1667,0.3056):
    P,tags=tagged(eta,0.0); o=P[:,0]+P[:,2]
    rm=lp(P,o,0.25,0.25); rd=lp(P,o+q*P[:,1],0.0,0.2499); rn=lp(P,o+q*P[:,1],0.2501,5.0)
    m=0.75*q+0.5+rm.fun; d=0.5+0.5*q+rd.fun; nn=0.5+0.5*q+rn.fun
    print("eta=%.4f  main %.4f dual %.4f nos %.4f -> X=%.4f ; MDS needs X <= %.4f"%(eta,m,d,nn,min(nn,max(m,d)),1.2986-0.375*eta))
    for nm,r in (("main",rm),("dual",rd),("nos",rn)):
        for i in np.where(r.x>1e-7)[0]: print("   %s w=%.3f pt(x,y,b)=%s twist=%s kind=%s polar=%s local=%s b_loc=%.4f"%(nm,r.x[i],np.round(P[i],4),*tags[i]))
