# Sextic Gauss sums for all squarefree primary s (split and inert primes) via twisted multiplicativity.
import sys, numpy as np
sys.path.insert(0, '../sextic/mmnum')
from eis import prime_ideals, mul, norm, ZERO, chi_prime, _fq2_pow, Z6
from gauss6 import chi_table_split, gauss_prime
def gauss_inert(pr):
    q = pr['q']; tot = 0j
    for x0 in range(q):
        for y0 in range(q):
            if x0==0 and y0==0: continue
            e = chi_prime(pr, (x0, y0))
            tot += np.exp(1j*np.pi*e/3)*np.exp(-2j*np.pi*y0/q)
    return tot
def all_primes(Ymax):
    prs = prime_ideals(int(Ymax)+1)
    gs = []
    for p in prs:
        if p['kind']=='split': gs.append(gauss_prime(p, chi_table_split(p)))
        else: gs.append(gauss_inert(p))
    return prs, gs
def squarefree_gauss(Ymax):
    prs, gs = all_primes(Ymax)
    cache = {}
    def cx(i,k):
        key=(i,k)
        if key not in cache: cache[key]=chi_prime(prs[i], prs[k]['g'])
        return cache[key]
    out = {}
    def dfs(start, s_g, s_N, gval, plist):
        if s_N > 1: out[s_g] = (s_N, gval/np.sqrt(s_N), tuple(plist))
        for i in range(start, len(prs)):
            P = prs[i]['P']
            if s_N*P > Ymax: break
            cs = 0
            for k in plist: cs += cx(k,i) + cx(i,k)
            dfs(i+1, mul(s_g, prs[i]['g']), s_N*P, gval*gs[i]*np.exp(1j*np.pi*(cs%6)/3), plist+[i])
    dfs(0,(1,0),1,1.0+0j,[])
    return prs, gs, out
if __name__=="__main__":
    prs, gs, out = squarefree_gauss(400)
    inert=[(p['g'],abs(g)**2,p['P']) for p,g in zip(prs,gs) if p['kind']=='inert']
    print('inert |g|^2 vs N:', inert[:4]); print('count sqfree', len(out))
