exec(open('region.py').read().split("ds=np.linspace")[0])
m_lo=max(.41,(Bf-.04*L)/2); m_hi=min(.44,(Bf-.034*L)/2,(1-.034*4.5)/2)
print('window',m_lo,m_hi)
def cross_old(p,k,t):
    b=p/(3*k); return ((2-2*b)*t-p+b)/(3-p-2*b)
def Pold(dl,p,k,m): return 1-2*dl*m-2*dl*p*(1-2*m)/(6*k)
def shortN(dl,p,t): return 1-dl*(p+(1-p)*rN(p,t))
import numpy as np
def check(d1,d2,p1,p2,n=60):
    worst={}
    def upd(k,v):
        worst[k]=min(worst.get(k,9),v)
    for dl in np.linspace(d1,d2,n):
        for p in np.linspace(p1,p2,n):
            t=tN(dl,p); r=rN(p,t); m=t-r
            upd('t-1',t-1); upd('1.5-t',1.5-t)
            upd('r-23/37',r-23/37); upd('1-r',1-r)
            upd('m-mlo',m-m_lo); upd('mhi-m',m_hi-m)
            for k in (k_,0.75):
                upd('cross-(t-mhi)',cross_old(p,k,t)-(t-m_hi))
                upd('short-Pold(mhi)',shortN(dl,p,t)-Pold(dl,p,k,m_hi))
            upd('Rold-Rnew',Rold(k_,dl,p)-Rnew(dl,p))
            upd('G new h',T-1e-6-G(hf,newC(hf),dl,p,Rnew(dl,p)))
            for name,d,C,v in rows[1:]:
                upd('G '+name, T-v-G(d,C,dl,p,Rnew(dl,p)))
    for k,v in worst.items(): print(f'{k:20s} {v:.3e}')
k_=k
check(.44,.49,.49,.5)
