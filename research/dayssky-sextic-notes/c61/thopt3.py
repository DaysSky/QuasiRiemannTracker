import sys
sys.argv=['x','new']
src=open('thopt2.py').read().split("m, h = best_h(0.874861)")[0]
exec(src)
from scipy.optimize import minimize_scalar
def hmax(th):
    ell=11/3-4*th
    # bmax=2y/dm-1<=4/25 with dm=h-1e-3 ; y=h-2ell
    h1=(4*ell-1.16e-3)/0.84
    # zmax<=37/1000 : nmin=1-y/dm >= (B-0.037Lam)/2
    nm=(B-0.037*Lam)/2
    # y <= (1-nm) dm ; h-2ell <= (1-nm)(h-1e-3) ; h*nm <= 2ell-(1-nm)1e-3
    h2=(2*ell-(1-nm)*1e-3)/nm
    return min(h1,h2)
def best_h_c(th):
    hm=hmax(th)
    res = minimize_scalar(lambda h: -margin(th, h), bounds=(0.785, hm), method='bounded', options={'xatol':1e-9})
    return -res.fun, res.x
for th in (0.874861,0.87485,0.874846):
    print(th, hmax(th), best_h_c(th))
lo,hi=0.87483,0.874861
for _ in range(18):
    mid=(lo+hi)/2; m,h=best_h_c(mid)
    if m>=0: hi=mid; hb=h
    else: lo=mid
print('theta_min constrained',hi,'h',hb,'hmax',hmax(hi))
