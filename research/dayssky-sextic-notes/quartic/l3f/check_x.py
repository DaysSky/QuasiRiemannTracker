import numpy as np
src=open('recur6.py').read().split('if __name__=="__main__":')[0]
exec(src)
tab=None; lam=1.5253
loc=local_points(lam,tab)
q=2/3
for eta in (0.0,0.2222,0.2778,0.3056,0.3333):
    print("eta=%.4f  X(eta,0)=%.4f   X(eta,0.01)=%.4f"%(eta,Xq(loc,q,eta,0.0,tab),Xq(loc,q,eta,0.01,tab)),flush=True)
