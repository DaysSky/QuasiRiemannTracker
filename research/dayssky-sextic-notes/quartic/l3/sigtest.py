import numpy as np
d=np.fromfile('pg1e6.bin').reshape(-1,5); a,b,p,Gr,Gi=d.T
G=Gr+1j*Gi; pi=a+1j*b; al=pi/np.abs(pi); T=al*G/np.sqrt(p)
print("check T^2=-alpha: max err", np.abs(T**2+al).max())
print("check G(conj)=conj G: ", end='')
key={ (int(x),int(y)):g for x,y,g in zip(a,b,G)}
print(max(abs(key[(int(x),-int(y))]-np.conj(g)) for x,y,g in zip(a[:2000],b[:2000],G[:2000])))
s=b>0; a,b,p,T,pi=a[s],b[s],p[s],T[s],pi[s]
th=np.angle(pi); u=T*np.exp(-0.5j*th); print("Re part max", np.abs(u.real).max())
sg=np.sign(u.imag); n=len(sg); print("n",n,"mean sign",sg.mean(), "sqrt n", np.sqrt(n))
# correlations with angular characters and residue classes
best=[]
for j in range(-12,13):
    for mod in [1]:
        v=(sg*np.exp(1j*j*th)).sum()/np.sqrt(n); best.append((abs(v),'ang j=%d'%j))
# residue classes mod small rational integers m (a,b mod m)
for m in [3,4,5,7,8,9,12,13,16]:
    cls={}
    for x,y,s_ in zip(a.astype(int)%m,b.astype(int)%m,sg): cls.setdefault((x,y),[]).append(s_)
    for k,v in cls.items():
        if len(v)>50: best.append((abs(np.sum(v))/np.sqrt(len(v)),'mod %d class %s n=%d'%(m,k,len(v))))
# quadratic/quartic symbols (pi/q) for small rational primes q=3,7,11 (inert) via norm residues and split q
def leg(x,q): x%=q; return 0 if x==0 else (1 if pow(int(x),(q-1)//2,q)==1 else -1)
for q in [3,7,11,19,23]:
    v=np.array([leg(int(pp),q) for pp in p]); best.append((abs((sg*v).sum())/np.sqrt(n),'leg(N,%d)'%q))
# correlate with N mod 16 and arg-sector
for m in [8,16,32]:
    for r in range(1,m,4):
        sel=(p%m==r)
        if sel.sum()>50: best.append((abs(sg[sel].sum())/np.sqrt(sel.sum()),'N mod %d = %d'%(m,r)))
best.sort(reverse=True)
for x in best[:15]: print("%.2f  %s"%x)
print("total tests",len(best))
# sector test: sign vs angle bins
bins=np.linspace(0,np.pi,13); idx=np.digitize(th,bins)
print("by angle sector:",[round(sg[idx==k].mean(),3) for k in range(1,13)])
