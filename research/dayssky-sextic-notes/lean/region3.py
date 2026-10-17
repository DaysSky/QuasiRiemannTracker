exec(open('region2.py').read().split("k_=k")[0])
k_=k
import numpy as np
for dl in np.linspace(.40,.52,13):
    row=[]
    for p in (.48,.49,.495,.4985,.5):
        t=tN(dl,p); r=rN(p,t); m=t-r
        row.append(f'{m-m_lo:+.2e}')
    print(f'{dl:.3f}',row)
for (d1,d2,p1,p2) in [(.445,.485,.495,.5),(.445,.485,.497,.5),(.44,.49,.498,.5),(.446,.484,.498,.5)]:
    print((d1,d2,p1,p2)); check(d1,d2,p1,p2,40)
