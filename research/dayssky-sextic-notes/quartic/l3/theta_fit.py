import numpy as np, sys
raw=open(sys.argv[1],'rb').read(); K=int(np.frombuffer(raw[:4],dtype=np.int32)[0])
pr=np.frombuffer(raw[4:4+16*K],dtype=np.float64).reshape(K,2)
rec=np.dtype([('N','<f8'),('tr','<f8'),('ti','<f8'),('a','<f8'),('b','<f8'),('e','u1',(K,))])
d=np.frombuffer(raw[4+16*K:],dtype=rec)
N=d['N']; T=d['tr']+1j*d['ti']; al=(d['a']+1j*d['b'])/np.sqrt(N); E=d['e'].astype(np.int64)
g=T*np.conj(al)                     # gtilde(conj chi_c)
mu=np.sign(np.real(T**2/al))        # T^2 = mu*alpha
print("K",K,[(int(x),int(y)) for x,y in pr],"#c",len(N), "mu check", np.abs(T**2-mu*al).max())
Wt=lambda x: np.where((x>0.5)&(x<1), np.exp(-1/np.clip((x-0.5)*(1-x),1e-300,None)+16), 0.0)
Zs=np.logspace(4,6,13)
def K1(idx, Z, mob=False):
    ex=np.zeros(len(N),np.int64); zero=np.zeros(len(N),bool)
    for k in idx: v=E[:,k]; zero|=(v==255); ex+=np.where(v==255,0,v)
    z=np.array([1,1j,-1,-1j])[ex&3]; z[zero]=0      # chi_c(m) (index m in the chi-bar convention)
    w=Wt(N/Z); coef=g*(mu if mob else 1)
    return (coef*z*w).sum(), np.sqrt(((w**2)*(~zero)).sum())
cases=[("1",[]),("pi5",[0]),("pi5^2",[0,0]),("pi5^3",[0,0,0]),("pi13",[3]),("pi5*pi13",[0,3]),("3 (inert)",[2])]
for name,idx in cases:
    out=[]
    for Z in Zs:
        k1,sd=K1(idx,Z); out.append(k1/Z**0.75)
    out=np.array(out)
    print("%-10s K1/Z^{3/4}: "%name+" ".join("%6.3f"%abs(v) for v in out)+"   noise level ~%.3f"%(0.155*Zs[-1]**-0.25))
# S13 follow-up: Moebius rows at small index (= conj R(m) up to index convention) on a fine Z grid; decade-binned mean squares
Zf=np.logspace(3.5,6,121)
for name,idx in [("1",[]),("pi5",[0]),("pi13",[3])]:
    vals=np.array([abs(K1(idx,Z,mob=True)[0])**2 for Z in Zf]); sds=np.array([K1(idx,Z,mob=True)[1]**2 for Z in Zf])
    for lo,hi in [(10**3.5,10**4.25),(10**4.25,10**5),(10**5,10**5.5),(10**5.5,10**6)]:
        s=(Zf>=lo)&(Zf<hi)
        print("Moebius index %-5s Z in [%.0e,%.0e): mean |K|^2/Z^{3/4} = %.3f   mean |K|^2/diag = %.3f"%(name,lo,hi,(vals[s]/Zf[s]**0.75).mean(),(vals[s]/sds[s]).mean()))
