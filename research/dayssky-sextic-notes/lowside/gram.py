# Gram test: G = sum_m Phi(N(m)/Q) |F(m)|^2,  F(m) = sum_s W(N(s)/Y) gt(chi_s) conj(chi_s(-m)),
# compared with diagonal D = sum_m Phi * sum_s |W gt|^2.   Paper's Gram bound predicts excess ~ P_a^{1/6}, P_a = Y^2/Q.
import sys, math, numpy as np, time
sys.path.insert(0, '../sextic/mmnum')
from gauss6 import chi_table_split, gauss_prime, chi_of
from eis import prime_ideals, mul, norm, ZERO
from mmdata import Ideals, rows, run_kernel
Y = float(sys.argv[1]); expo = float(sys.argv[2]); Q = Y**expo
t0=time.time()
Ymax = 2.3*Y
prs = [p for p in prime_ideals(int(Ymax)+1) if p['kind']=='split']
tabs = [chi_table_split(p) for p in prs]
gs = [gauss_prime(p, t) for p,t in zip(prs,tabs)]
sq = {}
def dfs(start, s_g, s_N, gval, plist):
    if s_N > 1: sq[s_g] = gval/np.sqrt(s_N)
    for i in range(start, len(prs)):
        P = prs[i]['P']
        if s_N*P > Ymax: break
        cs=0
        for k in plist:
            cs += chi_of(prs[k], tabs[k], prs[i]['g']) + chi_of(prs[i], tabs[i], prs[k]['g'])
        dfs(i+1, mul(s_g, prs[i]['g']), s_N*P, gval*gs[i]*np.exp(1j*np.pi*(cs%6)/3), plist+[i])
dfs(0,(1,0),1,1.0+0j,[])
I = Ideals(int(Ymax)+2)
lo, hi = I.rng(Y, Ymax)
def bump(u):
    u=np.asarray(u,float); s=(u-1.25)/0.25; out=np.zeros_like(u); k=np.abs(s)<1; out[k]=np.exp(1-1/(1-s[k]**2)); return out
def mom(Wf,a):
    uu=np.linspace(0.5,4,200001); return np.trapezoid(Wf(uu)*uu**a, uu)
lam = mom(bump,-1/3)/mom(lambda u: bump(u/1.5),-1/3)
WK = lambda u: bump(u) - lam*bump(u/1.5)
specs=[]; diag={}
for name, Wf, rng in (('W0', bump, (Y, 1.5*Y)), ('Wk', WK, (Y, 2.25*Y))):
    l, h = I.rng(*rng) if name=='W0' else (lo, hi)
    if name=='Wk':
        l, h = I.rng(Y, min(2.25*Y, Ymax))
    w = np.zeros(h-l, complex)
    for j in range(l, h):
        g = I.g[j]
        if g in sq: w[j-l] = Wf(I.N[j]/Y)*np.conj(sq[g])
    diag[name] = float(np.sum(np.abs(w)**2))
    specs.append((l, h, w.real.copy())); specs.append((l, h, w.imag.copy()))
print('squarefree s', len(sq), 'ideals in range', hi-lo, 'prep', round(time.time()-t0,1), flush=True)
rx, ry = rows(int(Q), sixth_power_free=False)
Nm = rx*rx - rx*ry + ry*ry
phi = bump(Nm/(Q/1.5))            # smooth row ball weight on N(m) in [Q/1.5*1, Q/1.5*1.5]
keep = phi > 0
rx, ry, phi = rx[keep], ry[keep], phi[keep]
print('rows', len(rx), flush=True)
V = run_kernel(I, rx, ry, specs, f"gram{int(Y)}")
for idx,name in enumerate(('W0','Wk')):
    K = V[:,2*idx] + 1j*V[:,2*idx+1]
    G = np.sum(phi*np.abs(K)**2); D = np.sum(phi)*diag[name]
    print(f"{name}: Gram/diag = {G/D:.4f}", flush=True)
Pa = Y**2/Q
print(f"Y={Y:.0f} Q={Q:.0f} P_a={Pa:.2f} P_a^(1/6)={Pa**(1/6):.3f} P_a^2/Y={Pa**2/Y:.3f} total {time.time()-t0:.0f}s")
