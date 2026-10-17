import numpy as np
raw=open('T15e4.bin','rb').read(); K=np.frombuffer(raw[:4],dtype=np.int32)[0]
pr=np.frombuffer(raw[4:4+16*K],dtype=np.float64).reshape(K,2)
rec=np.dtype([('N','<f8'),('tr','<f8'),('ti','<f8'),('e','u1',(K,))])
d=np.frombuffer(raw[4+16*K:],dtype=rec); N=d['N']; T=d['tr']+1j*d['ti']; E=d['e']
Z=150000.; Q=Z**(2/3)
s=(N>Z/2)&(N<=Z); N=N[s]; T=T[s]; E=E[s].astype(np.int16)
x=N/Z; W=np.exp(-1/((x-0.5)*(1-x))+16)
pN=np.round(pr[:,0]**2+pr[:,1]**2)
rows=[]
def rdfs(start,n,fl):
    rows.append(list(fl))
    for k in range(start,K):
        if n*pN[k]>Q: break
        rdfs(k,n*pN[k],fl+[k])
order=np.argsort(pN); assert (order==np.arange(K)).all()
rdfs(0,1,[])
print("c count",len(N),"rows",len(rows))
ph=np.array([1,-1j,-1,1j])
for t in [0,3,30,300,3000]:
    w=W*np.exp(-1j*t*np.log(N))
    S=0; D=0
    for fl in rows:
        ex=np.zeros(len(N),np.int64); zero=np.zeros(len(N),bool)
        for k in fl:
            v=E[:,k]; zero|=(v==255); ex+=np.where(v==255,0,v)
        z=ph[ex&3]; z[zero]=0
        R=(T*z*w).sum(); S+=abs(R)**2; D+=((W**2)*(~zero)).sum()
    print("t=%5d  S/diag = %.4f"%(t,S/D))
