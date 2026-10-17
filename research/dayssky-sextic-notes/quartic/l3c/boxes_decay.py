exec(open('boxes_rest.py').read().split("res.sort")[0])
import numpy as np
for d in (0.0,0.01,0.02,0.04,0.06,0.08,0.1,0.15,0.2):
    w=max([r for r in res if r[1][0]<=q-d+1e-9], key=lambda r:r[0])
    print("x1<=q-%.2f: worst (L3) exponent %.4f at box %s (%s)"%(d,1.5*w[0],w[1],w[2]))
