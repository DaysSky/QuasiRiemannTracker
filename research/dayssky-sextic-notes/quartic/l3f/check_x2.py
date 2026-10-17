import numpy as np, json
src=open('recur6.py').read().split('if __name__=="__main__":')[0]
exec(src)
d=json.load(open('recur6_tab.json')); tab=np.interp(QS,np.array(d['QS']),np.array(d['tab']))
lam=1.5104
loc=local_points(lam,tab)
q=2/3; s=2.0035; H=1/3
for eta in (0.3056,0.3194,0.3333):
    X0=Xq(loc,q,eta,0.0,tab)
    budget=(2*s-4/3-eta-(H-eta)/4)/2
    print("eta=%.4f X=%.4f budget=%.4f generic(k0=eta)=%.4f"%(eta,X0,budget,2+(H-eta)/4),flush=True)
