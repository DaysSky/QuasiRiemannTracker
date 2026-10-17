import sys, math, numpy as np, time
from gauss6 import build, chi_of
sys.path.insert(0, '../sextic/mmnum')
from eis import mul, norm, ZERO, primary_assoc
Ymax = int(float(sys.argv[1]))
t0=time.time()
prs, tabs, gs = build(Ymax)
np_ = len(prs)
# small test rows m (any nonzero elements), use u = -m
ms = []
B=8
for a in range(-B,B+1):
    for b in range(-B,B+1):
        if (a,b)!=(0,0) and norm((a,b))<=40: ms.append((a,b))
U = [(-a,-b) for (a,b) in ms]
# chi_p(u) table: E[p, j]
E = np.array([[chi_of(prs[i], tabs[i], u) for u in U] for i in range(np_)], dtype=np.int16)
# cross symbols between primes: C[i,j] = chi_{p_i}(p_j)
print('primes', np_, 'prep', round(time.time()-t0,1), flush=True)
# enumerate squarefree products (DFS) with norm <= Ymax
Ns=[]; Gs=[]; Es=[]
def dfs(start, s_g, s_N, gval, eu, plist):
    if s_N > 1:
        Ns.append(s_N); Gs.append(gval); Es.append(eu.copy())
    for i in range(start, np_):
        P = prs[i]['P']
        if s_N*P > Ymax: break
        if any(prs[k]['P']==P for k in plist):   # avoid both primes above same rational prime? allowed actually; keep distinct ideals
            pass
        # cross symbol: chi_s(p_i) * chi_{p_i}(s)
        cs = 0; ok=True
        for k in plist:
            a = chi_of(prs[k], tabs[k], prs[i]['g']); b = chi_of(prs[i], tabs[i], prs[k]['g'])
            if a==ZERO or b==ZERO: ok=False; break
            cs += a + b
        if not ok: continue
        newg = gval*gs[i]*np.exp(1j*np.pi*(cs%6)/3)
        ne = eu.copy()
        bad = (E[i]==ZERO) | (ne==ZERO)
        ne = np.where(bad, ZERO, (ne + E[i]) % 6)
        dfs(i+1, mul(s_g, prs[i]['g']), s_N*P, newg, ne, plist+[i])
dfs(0, (1,0), 1, 1.0+0j, np.zeros(len(U), np.int16), [])
Ns=np.array(Ns,float); Gs=np.array(Gs); Es=np.array(Es)
print('squarefree s:', len(Ns), 'time', round(time.time()-t0,1), flush=True)
gt = Gs/np.sqrt(Ns)
# chi_s(-m) conj
chiv = np.where(Es==ZERO, 0, np.exp(-1j*np.pi*Es/3))   # conj chi
def bump(u):
    u=np.asarray(u,float); s=4*(u-1)-1  # support [1,1.5]? map [1,1.5]->[-1,1]
    s = (u-1.25)/0.25
    out=np.zeros_like(u); k=np.abs(s)<1; out[k]=np.exp(1-1/(1-s[k]**2)); return out
def mom(Wf, a):
    uu=np.linspace(0.5,4,200001); return np.trapezoid(Wf(uu)*uu**a, uu)
W0 = bump
lam = mom(lambda u: bump(u),-1/3)/mom(lambda u: bump(u/1.5),-1/3)
Wk = lambda u: bump(u) - lam*bump(u/1.5)
print('check kill moment', mom(Wk,-1/3), 'other moment', mom(Wk,0.0))
for Y in [Ymax/16, Ymax/8, Ymax/4, Ymax/2.5]:
    res=[]
    for name,Wf in (('W0',W0),('Wkill',Wk)):
        w = Wf(Ns/Y)
        F = (w*gt) @ chiv            # F(m) for each test m
        rnd = np.sqrt(np.sum(w**2))
        res.append((name, np.mean(np.abs(F))/rnd, np.max(np.abs(F))/rnd))
    print(f"Y={Y:9.0f} #s~{int(np.sum(bump(Ns/Y)>0))}", '  '.join(f"{n}: mean|F|/rand={a:.2f} max={b:.2f}" for n,a,b in res), flush=True)
