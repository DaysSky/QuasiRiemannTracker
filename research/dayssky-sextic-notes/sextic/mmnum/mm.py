# Task 2: mixed moment MM(r,m) = U^{-1} sum_{0<N(u)<=U, u 6pf} |M_u(U^r) S_u(U^m)|^2 for the sextic family.
#   M_u(D) = D^{-1/2} sum_{n primary odd} c(n) chi_n(u) W(Nn/D), c = mu (or random-sign replacements)
#   S_u(N) = N^{-1/2} sum_n chi_n(u) W(Nn/N), W = bump on [1/2,1].
# usage: python3 mm.py U [seed]
import sys, json, time, math
import numpy as np
from collections import defaultdict
from mmdata import Ideals, rows, run_kernel, W
from eis import mul

U = int(float(sys.argv[1])); seed = int(sys.argv[2]) if len(sys.argv) > 2 else 12345
t0 = time.time()
I = Ideals(U)
rx, ry = rows(U)
nrows = len(rx)
rng = np.random.default_rng(seed)
sqf = (I.mu != 0)
eps_ind = rng.choice([-1, 1], size=len(I.N)) * sqf               # independent random signs on squarefree n
fp = rng.choice([-1, 1], size=len(I.primes))                        # random multiplicative (Rademacher) at primes
eps_mul = np.zeros(len(I.N), np.int64); eps_mul[0] = 1
for i in range(1, len(I.N)):
    if sqf[i]: eps_mul[i] = -1 * 0 + fp[I.spf[i]] * eps_mul[I.cof[i]]
eps_mul *= sqf
COEF = {'mu': I.mu.astype(float), 'rnd': eps_ind.astype(float), 'rmul': eps_mul.astype(float)}

RCONF = {1.0: [0.0, 0.15, 0.25, 1/3, 0.4, 0.5], 0.9: [0.0, 0.15, 0.25, 0.4, 0.5],
         1.1: [0.15, 0.25, 0.4, 0.5], 1.25: [0.25, 0.4, 0.5]}
configs = []
for tot, ms in RCONF.items():
    for m in ms:
        configs.append((round(tot - m, 6), m))
specs = []; spec_id = {}
def add_spec(key, length, coef):
    if key in spec_id: return
    lo, hi = I.rng(length/2, length)
    w = coef[lo:hi] * W(I.N[lo:hi] / length)
    spec_id[key] = len(specs); specs.append((lo, hi, w))
for r, m in configs:
    D = U**r
    for v in COEF: add_spec(('M', v, r), D, COEF[v])
    if m > 0: add_spec(('S', m), U**m, np.ones(len(I.N)))
for m in (0.25, 0.5, 0.75): add_spec(('S', m), U**m, np.ones(len(I.N)))
print(f"U={U} rows={nrows} ideals={len(I.N)} primes={len(I.primes)} specs={len(specs)} prep {time.time()-t0:.1f}s", flush=True)
t1 = time.time()
V = run_kernel(I, rx, ry, specs, f"U{U}")
print(f"kernel {time.time()-t1:.1f}s", flush=True)

# --- diagonal prediction ---
gp = np.array([(1 - 1/P) / (1 - P**-6.0) for P in I.Pnorm])
FAC = {}
def fac(i):
    if i not in FAC: FAC[i] = I.factor(i)
    return FAC[i]
def diag(specA, specB, D, N):
    """(|R|/U) (DN)^{-1} sum_{k ~6 k'} A_k A_k' prod_{p | kk'} g_p ; specB None -> k = n1."""
    A = defaultdict(float); F = {}
    loA, hiA, wA = specA
    if specB is None:
        it = [((i,), wA[i - loA]) for i in range(loA, hiA) if wA[i - loA] != 0]
    else:
        loB, hiB, wB = specB
        it = [((i, j), wA[i - loA] * wB[j - loB]) for i in range(loA, hiA) if wA[i - loA] != 0
              for j in range(loB, hiB) if wB[j - loB] != 0]
    for key, wt in it:
        g = (1, 0); fl = []
        for i in key:
            g = mul(g, I.g[i]); fl += fac(i)
        A[g] += wt
        F[g] = fl
    classes = defaultdict(list)
    for g, a in A.items():
        if a == 0: continue
        cnt = defaultdict(int)
        for p in F[g]: cnt[p] += 1
        ck = tuple(sorted((p, v % 6) for p, v in cnt.items() if v % 6))
        classes[ck].append((a, frozenset(cnt)))
    tot = 0.0
    for ck, L in classes.items():
        for a1, s1 in L:
            for a2, s2 in L:
                tot += a1 * a2 * np.prod(gp[list(s1 | s2)]) if (s1 | s2) else a1 * a2
    return nrows / U * tot / (D * N)

res = {'U': U, 'rows': nrows, 'seed': seed, 'configs': []}
def val(key, length):
    return V[:, spec_id[key]] / math.sqrt(length)
# bounded rows: u supported on S (units * lambda^i * 2^j): flag them
def s_supported(x, y):
    a, b = int(x), int(y)
    while (a*a - a*b + b*b) % 3 == 0:
        # divide by lambda = 1 + 2w : u * conj(lambda) / 3, conj(lambda) = (1-2) + (-2)w = -1 - 2w
        a, b = (a*(-1) - b*(-2)) // 3, (a*(-2) + b*(-1) - b*(-2)) // 3
    while a % 2 == 0 and b % 2 == 0:
        a, b = a // 2, b // 2
    return a*a - a*b + b*b == 1
bounded = np.array([s_supported(x, y) for x, y in zip(rx, ry)])
for r, m in configs:
    D, N = U**r, U**m
    for v in COEF:
        Mv = val(('M', v, r), D)
        Sv = val(('S', m), N) if m > 0 else np.ones(nrows)
        prod2 = np.abs(Mv * Sv)**2
        MM = prod2.sum() / U
        dg = diag(specs[spec_id[('M', v, r)]], specs[spec_id[('S', m)]] if m > 0 else None, D, N)
        rec = dict(r=r, m=m, tot=round(r + m, 4), var=v, MM=MM, diag=dg, ratio=MM / dg,
                   maxrow_over_U=prod2.max() / U, bounded_rows_share=prod2[bounded].sum() / U / MM,
                   se=float(prod2.std()*math.sqrt(nrows)/U), z=float((MM-dg)/(prod2.std()*math.sqrt(nrows)/U)) if prod2.std()>0 else 0.0,
                   nM=int(np.count_nonzero(specs[spec_id[('M', v, r)]][2])),
                   nS=int(np.count_nonzero(specs[spec_id[('S', m)]][2])) if m > 0 else 0,
                   EM2=float(np.mean(np.abs(Mv)**2)), ES2=float(np.mean(np.abs(Sv)**2)))
        res['configs'].append(rec)
        print(f"r={r:.3f} m={m:.3f} {v:4s} MM={MM:9.4f} diag={dg:9.4f} ratio={MM/dg:7.4f} z={rec['z']:+.2f} "
              f"max|MS|^2/U={rec['maxrow_over_U']:.2e} nM={rec['nM']} nS={rec['nS']}", flush=True)
# calibration: plain x plain at total length 1 (fourth moment at m=1/2, and (1/4,3/4))
for m1, m2 in ((0.5, 0.5), (0.25, 0.75)):
    S1 = val(('S', m1), U**m1); S2 = val(('S', m2), U**m2)
    prod2 = np.abs(S1 * S2)**2
    MM = prod2.sum() / U
    dg = diag(specs[spec_id[('S', m1)]], specs[spec_id[('S', m2)]], U**m1, U**m2)
    se = float(prod2.std()*math.sqrt(nrows)/U)
    res['configs'].append(dict(r=None, m=m1, m2=m2, var='plainxplain', MM=MM, diag=dg, ratio=MM / dg, se=se, z=(MM-dg)/se,
                               maxrow_over_U=prod2.max() / U, bounded_rows_share=prod2[bounded].sum() / U / MM))
    print(f"plain x plain m1={m1} m2={m2} MM={MM:.4f} diag={dg:.4f} ratio={MM/dg:.4f} z={(MM-dg)/se:+.2f} max/U={prod2.max()/U:.2e}", flush=True)
json.dump(res, open(f"res_U{U}_s{seed}.json", 'w'), indent=1)
print(f"total {time.time()-t0:.1f}s")
