import numpy as np, json, sys
G={}; exec(open('recur9.py').read().split('if __name__')[0],G); globals().update(G)
d=json.load(open('recur9_tab.json')); tab=np.array(d['tab']); lam=d['lam']
print("lam",lam,"s(2/3)",np.interp(2/3,QS,tab))
loc=local_points(lam,tab); bc={}
q=2/3; H=2*q-1; s=float(sys.argv[1]) if len(sys.argv)>1 else 2.001
etas=np.linspace(0,H,13)
for g in etas:
    emax=(H-g)/4
    ee=emax
    for t in np.linspace(0,emax,41):
        if all(2-q+(H-g)/4-u+lXi(g,q-u)<=s+1e-9 for u in np.linspace(t,emax,9) if u>0): ee=t; break
    xs=[Xq(loc,q,g,t,tab,bc) for t in np.linspace(0,ee,3)]
    mds=(1+(1-q)+g+(H-g)/4+2*max(xs))/2
    print("eta=%.4f e0=%.4f X(eps)=%s  MDS s=%.5f  gen=%.5f"%(g,ee,np.round(xs,5),mds,2-q+(H-g)/4+lXi(g,q)))
