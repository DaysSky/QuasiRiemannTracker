import numpy as np, sys
raw=open(sys.argv[1],'rb').read(); K=np.frombuffer(raw[:4],dtype=np.int32)[0]
pr=np.frombuffer(raw[4:4+16*K],dtype=np.float64).reshape(K,2)
rec=np.dtype([('N','<f8'),('tr','<f8'),('ti','<f8'),('e','u1',(K,))])
d=np.frombuffer(raw[4+16*K:],dtype=rec)
N=d['N']; T=d['tr']+1j*d['ti']; E=d['e']
print("K=%d test primes:"%K, [(int(a),int(b)) for a,b in pr], " #c=%d"%len(N))
def W(x):
    out=np.zeros_like(x); s=(x>0.5)&(x<1); out[s]=np.exp(-1/((x[s]-0.5)*(1-x[s]))+16); return out
def row(m_idx, Z):
    # m given as list of test-prime indices (with multiplicity); chi_c(m)-bar weight
    w=W(N/Z); ex=np.zeros(len(N),dtype=np.int64); zero=np.zeros(len(N),bool)
    for k in m_idx:
        v=E[:,k]; zero|=(v==255); ex+=np.where(v==255,0,v)
    z=np.array([1,-1j,-1,1j])[ex&3]; z[zero]=0
    R=(T*z*w).sum(); D=((w**2)*(~zero)).sum(); return R,D
Zs=np.unique(np.round(np.logspace(3.5,6,26))).astype(float)
for name,m in [("m=1",[]),("m=pi5",[0]),("m=pi13",[2]),("m=pi5*pi13",[0,2])]:
    print(name)
    for Z in Zs:
        R,D=row(m,Z)
        print("  Z=%9.0f |R|=%8.2f sqrtDiag=%7.2f |R|/sqrtD=%5.2f |R|/Z^{3/8}=%5.2f"%(Z,abs(R),np.sqrt(D),abs(R)/np.sqrt(D),abs(R)/Z**0.375))
