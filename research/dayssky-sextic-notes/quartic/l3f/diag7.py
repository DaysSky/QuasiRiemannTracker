import numpy as np, json
src=open('recur7.py').read().split('if __name__=="__main__":')[0]
exec(src)
d=json.load(open('recur7_tab.json')); tab=np.interp(QS,np.array(d['QS']),np.array(d['tab'])); lam=d['lam']
loc=local_points(lam,tab)
def points_tagged(loc,tK,tE,tab):
    pts=[];tags=[]
    lb={pol:[local_b(it,tK,tE,tab,pol) for it in loc[pol]] for pol in (False,True)}
    names={}
    for k,(M,c) in ch.items():
        pol=POLAR[k]; Mi=np.linalg.inv(M); tw=tuple((M@np.array([1,0])).astype(int))
        for it,b in zip(loc[pol],lb[pol]):
            v=Mi@np.array(it[0]); pts.append((v[0],v[1],b+(c@v)*(tK+tE))); tags.append((tw,it[1],pol))
    return np.array(pts),tags
q=2/3
for eta in (0.25,0.3056,1/3):
    P,tags=points_tagged(loc,eta,0.0,tab); o=P[:,0]+P[:,2]
    n=len(P); A_ub=np.vstack([P[:,1],-P[:,1],-(3*P[:,0]+P[:,1]),-(2*P[:,0]+P[:,1]),-(P[:,0]+P[:,1])])
    r=linprog(o,A_ub=A_ub,b_ub=[0.25,-0.25,0,0,0],A_eq=np.ones((1,n)),b_eq=[1],bounds=(0,None),method="highs")
    print("eta=%.4f main X=%.4f lam=%.4f"%(eta,0.75*q+0.5+r.fun,lam))
    for i in np.where(r.x>1e-9)[0]: print("   w=%.3f pt=%s twist=%s kind=%s polar=%s"%(r.x[i],np.round(P[i],3),tags[i][0],tags[i][1],tags[i][2]))
