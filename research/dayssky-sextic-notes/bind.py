import numpy as np, sys
sys.path.insert(0,'.')
exec(open('fastopt.py').read().split("for cap in")[0])
th=0.8748446; h=0.7952291
ell=11/3-4*th; y=h-2*ell; k=2*th-1
R=Rold(k,dl,P); R=np.where(OK,np.minimum(R,RN),R)
names=['new h','junction h-1e-3','half 1/2']
for (d,C,v),nm in zip(((h,(2*y-(2-Bc)*h)/L,1e-6),(h-1e-3,min(h-1e-3,2*y-h+1e-3)/(6*k),1e-5),(.5,min(.5,2*y-.5)/(6*k),1e-3)),names):
    G=(1+dl)/2+13*h/75-y/2-ell/2+ell*dl*P+d*(1/4+3*R/4-17/50)-dl*P*C/2
    M=th-v-G; i=np.unravel_index(np.argmin(M),M.shape)
    print(f'{nm:18s} margin {M[i]:+.3e} at delta={dl[i[0],0]:.4f} p={P[0,i[1]]:.4f} R={R[i]:.5f} usesNew={bool(OK[i] and RN[i]<Rold(k,dl,P)[i])}')
    # second-best region (exclude neighborhood)
    Mm=M.copy(); Mm[max(0,i[0]-60):i[0]+60,:]=9
    j=np.unravel_index(np.argmin(Mm),M.shape); print(f'{"":18s} next    {Mm[j]:+.3e} at delta={dl[j[0],0]:.4f} p={P[0,j[1]]:.4f}')
