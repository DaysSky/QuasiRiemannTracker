# Sublattice mean squares of sextic Kubota sums:  M(m) = mean_{N(h)<=H} |sum_a w(a) gt(h m^4, a)|^2  over moduli a coprime to m,
# against the untwisted M = mean_h |sum_a w(a) gt(h,a)|^2 over the same moduli.  w = single smooth weight, or the convolution W1*W2.
import sys, math, numpy as np
sys.path.insert(0, '.')
from kub import build, coeffs, W, g1_prime
from g6 import ZETA6, mul, nrm
Y = int(float(sys.argv[1])); H = int(float(sys.argv[2])); mode = sys.argv[3] if len(sys.argv) > 3 else 'single'
ps, mods = build(Y)
for P in ps: P.g1 = g1_prime(P)
# weights on moduli in [Y/2, Y] (single) or all moduli <= Y with convolution weight (product)
if mode == 'single':
    sel = [(n, idx) for (n, idx) in mods if Y/2 <= n <= Y]
    wt = np.array([W(n/Y) for (n, idx) in sel])
else:
    Y1 = math.sqrt(Y)  # two factors of norm ~ sqrt(Y) each
    import itertools
    sel = [(n, idx) for (n, idx) in mods if Y/4 <= n <= Y]
    wt = []
    for (n, idx) in sel:
        s = 0.0
        # sum over ordered factorizations a = a1 a2 (subsets of the prime list)
        for r in range(len(idx)+1):
            for sub in itertools.combinations(idx, r):
                n1 = 1
                for i in sub: n1 *= ps[i].p
                n2 = n//n1
                s += float(W(n1/Y1))*float(W(n2/Y1))
        wt.append(s)
    wt = np.array(wt); keep = wt > 0; sel = [s for s, k in zip(sel, keep) if k]; wt = wt[keep]
# gt(1,a) for selected moduli
c = np.zeros(len(sel), complex)
for j, (n, idx) in enumerate(sel):
    v = complex(1); k = 0
    for i in idx: v *= ps[i].g1
    for a in range(len(idx)):
        for b in range(a+1, len(idx)):
            k += ps[idx[a]].sym(ps[idx[b]].g) + ps[idx[b]].sym(ps[idx[a]].g)
    c[j] = wt[j]*v*ZETA6[k % 6]
# rows h: all nonzero elements with N(h) <= H
B = int(math.isqrt(4*H//3)) + 2
hs = [(x, y) for x in range(-B, B+1) for y in range(-B, B+1) if 0 < x*x - x*y + y*y <= H]
hx = np.array([h[0] for h in hs]); hy = np.array([h[1] for h in hs])
used = sorted(set(i for (n, idx) in sel for i in idx))
# exponent of (h/pi)_6 for all h, per used prime; None (pi | h) -> mark
E = {}; Z0 = {}
for i in used:
    P = ps[i]; v = (hx + hy*P.r) % P.p
    tab = np.full(P.p, -1); 
    for x in range(1, P.p): tab[x] = P.k[pow(x, P.e, P.p)]
    e = tab[v]; Z0[i] = (e < 0); E[i] = np.where(e < 0, 0, e)
def chi_exp_rows(idx):
    tot = np.zeros(len(hs), int); zero = np.zeros(len(hs), bool)
    for i in idx: tot += E[i]; zero |= Z0[i]
    return tot, zero
CH = np.zeros((len(sel), len(hs)), complex)
for j, (n, idx) in enumerate(sel):
    tot, zero = chi_exp_rows(idx)
    CH[j] = np.where(zero, 0, np.exp(-2j*np.pi*tot/6))   # conj(chi_a(h))
print(f"Y={Y} H={H} mode={mode}: moduli {len(sel)}, rows {len(hs)}, sum|w|^2 = {np.sum(np.abs(c)**2):.2f}", flush=True)
for mi in range(0, 14):
    Pm = ps[mi]; m = Pm.g
    if Pm.p > 200: break
    cop = np.array([mi not in idx for (n, idx) in sel])
    # chi_a(m) exponent for coprime moduli
    km = np.array([sum(ps[i].sym(m) for i in idx) % 6 if cop[j] else 0 for j, (n, idx) in enumerate(sel)])
    cc = c*cop
    K0 = cc @ CH                                   # untwisted, moduli coprime to m
    Km = (cc*np.exp(-2j*np.pi*4*km/6)) @ CH        # gt(h m^4, a) = conj(chi_a(m^4)) gt(h,a)
    M0 = np.mean(np.abs(K0)**2); Mm = np.mean(np.abs(Km)**2)
    # controls: twist by a random character pattern of the same order structure (random sixth roots per prime)
    rng = np.random.default_rng(mi)
    rk = rng.integers(0, 6, len(ps)); kr = np.array([sum(rk[i] for i in idx) % 6 for (n, idx) in sel])
    Mr = np.mean(np.abs((cc*np.exp(2j*np.pi*kr/6)) @ CH)**2)
    print(f"  N(m)={Pm.p:4d}: M(m)/M0 = {Mm/M0:.3f}   random-twist control {Mr/M0:.3f}   M0/sum|w|^2 = {M0/np.sum(np.abs(cc)**2):.3f}", flush=True)
