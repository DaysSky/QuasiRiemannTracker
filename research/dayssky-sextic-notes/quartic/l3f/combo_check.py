import numpy as np, json
src=open('recur6.py').read().split('if __name__=="__main__":')[0]
exec(src)
d=json.load(open('recur6_tab.json')); tab=np.interp(QS,np.array(d['QS']),np.array(d['tab']))
for lam in (1.5253,1.5104,1.5053,1.5):
    loc=local_points(lam,tab)
    P=points(loc,1/3,0.0,tab); o=P[:,0]+P[:,2]
    n=len(P); A_ub=np.vstack([P[:,1],-P[:,1],-(3*P[:,0]+P[:,1]),-(2*P[:,0]+P[:,1]),-(P[:,0]+P[:,1])])
    r=linprog(o,A_ub=A_ub,b_ub=[0.25,-0.25,0,0,0],A_eq=np.ones((1,n)),b_eq=[1],bounds=(0,None),method="highs")
    idx=np.where(r.x>1e-9)[0]
    print("lam=%.4f main X(1/3,0)=%.4f  combo:"%(lam,1+r.fun))
    for i in idx: print("    w=%.3f  (x,y,b)=%s"%(r.x[i],np.round(P[i],4)))
    dual=lpmin(P,o+(2/3)*P[:,1],0.0,0.2499); print("    dual exp=%.4f"%(5/6+dual))
