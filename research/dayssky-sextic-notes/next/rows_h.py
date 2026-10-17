import numpy as np
exec(open('/home/user/QuasiRiemannTracker/proofs/detector-affine-c0125154/scripts/fastopt.py').read().split("for cap in")[0])
th=0.8748446
ell=11/3-4*th; print('hcap',hcap(th))
for h in np.linspace(0.790,0.800,11):
    y=h-2*ell; k=2*th-1
    R=Rold(k,dl,P); R=np.where(OK,np.minimum(R,RN),R)
    s=[]
    for d,C,v in ((h,(2*y-(2-Bc)*h)/L,1e-6),(h-1e-3,min(h-1e-3,2*y-h+1e-3)/(6*k),1e-5),(.5,min(.5,2*y-.5)/(6*k),1e-3)):
        G=(1+dl)/2+13*h/75-y/2-ell/2+ell*dl*P+d*(1/4+3*R/4-17/50)-dl*P*C/2
        M=th-v-G; i=np.unravel_index(np.argmin(M),M.shape)
        s.append('%+.2e@(%.3f,%.3f)'%(M[i],dl[i[0],0],P[0,i[1]]))
    print('%.4f'%h,' '.join(s))
