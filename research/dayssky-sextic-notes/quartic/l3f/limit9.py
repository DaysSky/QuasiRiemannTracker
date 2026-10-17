import numpy as np, json
G={}; exec(open('recur9.py').read().split('if __name__')[0],G); globals().update(G)
d=json.load(open('recur9_tab.json')); tabF=np.array(d['tab'])
def argmax_block(x,y,tK,tE,tab,th):
    C=GC[:,None]*np.ones_like(U)[None,:]; tD=th*(tK+tE+C)*U[None,:]
    v=blocks_val(x,y,C,tD,tK,sigma_vec(C,tD,tab)); i=np.unravel_index(np.argmax(v),v.shape)
    return float(v[i]),float(C[i]),float(tD[i])
for name,tab in (("S49 table",None),("recur9 final",tabF)):
    print("==",name, "s(2/3)=",S(2/3,tab))
    for dlt in (1e-2,1e-3,1e-4):
        x=1/3+dlt
        for tK in (0.1,0.2,1/3):
            bs=[(th,b_twist(x,0.0,tK,0.0,tab,{},th=th,polar=False)) for th in THS if okd_theta(x,0.0,th,tab,False)]
            th,b=min(bs,key=lambda z:z[1]) if bs else (None,None)
            print("  x=1/3+%g tK=%.3f: feasible th=%s  best b=%s at th=%s  argmax block=%s"%(dlt,tK,[round(t,1) for t,_ in bs],None if b is None else round(b,5),th,None if th is None else np.round(argmax_block(x,0.0,tK,0.0,tab,th),4)))
