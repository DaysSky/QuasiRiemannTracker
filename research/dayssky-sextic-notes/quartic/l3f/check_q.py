import numpy as np, json, sys
src=open('recur6.py').read().split('if __name__=="__main__":')[0]
exec(src)
d=json.load(open('recur6_tab.json')); tab=np.interp(QS,np.array(d['QS']),np.array(d['tab'])); lam=d['lam']
q=float(sys.argv[1]); s=float(sys.argv[2])
loc=local_points(lam,tab)
H=2*q-1
print("q=%.4f claimed s=%.4f lam=%.4f H=%.4f"%(q,s,lam,H))
for eta in np.linspace(0,H,9):
    emax=(H-eta)/4
    # e-tail threshold
    e0=emax
    for t in np.linspace(0,emax,41):
        if all(2-q+(H-eta)/4-u+lXi(eta,q-u)<=s+1e-9 for u in np.linspace(t,emax,9) if u>0): e0=t; break
    xs=[Xq(loc,q,eta,t,tab) for t in (0.0,e0)]
    budget=(2*s-1-(1-q)-eta-(H-eta)/4)/2
    print("  eta=%.4f e0=%.4f X(eta,0)=%.4f X(eta,e0)=%.4f budget=%.4f  ok=%s  generic=%.4f"%(eta,e0,xs[0],xs[1],budget,max(xs)<=budget+1e-9,2-q+(H-eta)/4+lXi(eta,q)),flush=True)
