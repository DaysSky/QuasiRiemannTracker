# Sextic Gauss sums over O=Z[w] for squarefree s built from split primes; twisted-multiplicativity check.
import sys, math, numpy as np
sys.path.insert(0, '../sextic/mmnum')
from eis import prime_ideals, mul, norm, Z6, ZERO
def chi_table_split(pr):
    P, r = pr['P'], pr['r']
    # exponent table e(t) for t in Z/P (t = image of x in O/p), via primitive root
    m = P - 1; fac=[]; t=m; d=2
    while d*d<=t:
        if t%d==0:
            fac.append(d)
            while t%d==0: t//=d
        d+=1
    if t>1: fac.append(t)
    g=2
    while any(pow(g,m//f,P)==1 for f in fac): g+=1
    h = pow(g, m//6, P)
    z6 = (1 + r) % P          # image of zeta6 = 1+w
    k0 = next(k for k in range(6) if pow(z6, k, P) == h)   # g^{m/6} = zeta6^{k0}
    tab = np.full(P, ZERO, np.int8)
    cur = 1
    for j in range(m):
        tab[cur] = (j*k0) % 6
        cur = cur*g % P
    return tab
def gauss_prime(pr, tab):
    P = pr['P']; a, b = pr['g']
    t = np.arange(P)
    z = np.exp(1j*np.pi*tab/3.0); z[tab==ZERO] = 0
    return np.sum(z*np.exp(-2j*np.pi*t*b/P))
def build(Ymax):
    prs = [p for p in prime_ideals(Ymax) if p['kind']=='split']
    tabs = {i: chi_table_split(p) for i,p in enumerate(prs)}
    gs = {i: gauss_prime(p, tabs[i]) for i,p in enumerate(prs)}
    return prs, tabs, gs
def chi_of(pr, tab, u):
    P, r = pr['P'], pr['r']; t = (u[0] + u[1]*r) % P
    return int(tab[t])
if __name__ == "__main__":
    prs, tabs, gs = build(400)
    # check |g|^2 = P and twisted multiplicativity on a few pairs by direct summation mod s=p1p2 (O/s = Z/N(s))
    for i in range(3):
        print(prs[i]['g'], prs[i]['P'], abs(gs[i])**2)
    for (i,j) in [(0,1),(1,4),(2,7)]:
        p1, p2 = prs[i], prs[j]
        if p1['P']==p2['P']: continue
        s = mul(p1['g'], p2['g']); N = norm(s); bs = s[1]
        t = np.arange(N)
        e1 = tabs[i][t % p1['P'] if False else ((t) % p1['P'])]
        # chi_s(t) = chi_p1(t) chi_p2(t); t is rational integer so image in O/p is t mod P
        e2 = tabs[j][t % p2['P']]
        ok = (e1!=ZERO)&(e2!=ZERO)
        val = np.where(ok, np.exp(1j*np.pi*((e1.astype(int)+e2.astype(int))%6)/3), 0)
        gdir = np.sum(val*np.exp(-2j*np.pi*t*bs/N))
        c12 = chi_of(p1, tabs[i], p2['g']); c21 = chi_of(p2, tabs[j], p1['g'])
        pred = gs[i]*gs[j]*np.exp(1j*np.pi*(c12+c21)/3)
        print('pair', p1['g'], p2['g'], 'direct', np.round(gdir,6), 'pred', np.round(pred,6), 'c12,c21', c12, c21)
