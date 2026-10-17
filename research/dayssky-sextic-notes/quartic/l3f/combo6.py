exec(open('fixa6.py').read().split("lam=13/8")[0])
P=sample(1.5365,"all")
def combo(P,obj,ylo,yhi):
    n=len(P); A_eq=np.ones((1,n))
    A_ub=np.vstack([P[:,1],-P[:,1],-(3*P[:,0]+P[:,1]),-(2*P[:,0]+P[:,1]),-(P[:,0]+P[:,1])])
    r=linprog(obj,A_ub=A_ub,b_ub=[yhi,-ylo,0,0,0],A_eq=A_eq,b_eq=[1],bounds=(0,None),method="highs")
    idx=np.where(r.x>1e-9)[0]; return r.fun,[(round(r.x[i],3),np.round(P[i],3)) for i in idx]
for eta in (0.0,0.236,0.333):
    print("eta",eta)
    print("  main:",combo(P,P[:,0]+eta*P[:,2],0.25,0.25))
    print("  dual:",combo(P,P[:,0]+(2/3)*P[:,1]+eta*P[:,2],0.0,0.2499))
