# RMS over Z-windows of |S_m(Z)|/sqrt(sum W^2) for m = 1, small squarefree m, and squares; Zmax large.
import sys, math, time, numpy as np
from autoA import prime_list, build_c, W
from gauss4 import gmul, gnorm, primary_assoc, primary_prime_over
Zmax=int(sys.argv[1]); t0=time.time()
P=prime_list(Zmax)
pr={p: primary_prime_over(p)[0] for p in (5,13,17,29)}
ms={'1':(1,0)}
for p in (5,13,17,29): ms[f'pi{p}']=pr[p]
for p in (5,13): ms[f'pi{p}^2']=gmul(pr[p],pr[p])
ms['q3']=primary_assoc((3,0)); ms['q3^2']=gmul(ms['q3'],ms['q3'])
Ns,A,CH,cs=build_c(Zmax,P,list(set(ms.values())))
print(f"#c={len(Ns)} build {time.time()-t0:.0f}s",flush=True)
edges=[30,100,300,1000,3000,10000,30000,Zmax]
print("window-RMS of |S_m(Z)|/||W||   windows:",edges)
for name,m in ms.items():
    ch=CH[m]; row=f"{name:>7} N={gnorm(m):4d}:"
    for lo,hi in zip(edges[:-1],edges[1:]):
        vals=[]
        for Z in np.geomspace(lo,hi,12):
            w=W(Ns/Z); nrm=math.sqrt(np.sum(w*w*(ch!=0)))
            vals.append(abs(np.sum(A*np.conj(ch)*w))**2/nrm**2)
        row+=f" {math.sqrt(np.mean(vals)):6.3f}"
    print(row,flush=True)
