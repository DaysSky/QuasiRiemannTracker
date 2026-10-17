# brute-force check of ff_lm.py for small q, d: per-row L-polynomial roots and M_u(k)
import sys, itertools, numpy as np
q, d = int(sys.argv[1]), int(sys.argv[2])
def pmod(a, b):  # coef lists low->high, b monic
    a = a[:]
    while len(a) >= len(b):
        c = a[-1]
        if c:
            for i in range(len(b)):
                a[len(a)-len(b)+i] = (a[len(a)-len(b)+i] - c*b[i]) % q
        a.pop()
    while a and a[-1] == 0: a.pop()
    return a
def pmul(a, b):
    if not a or not b: return []
    c = [0]*(len(a)+len(b)-1)
    for i, x in enumerate(a):
        for j, y in enumerate(b): c[i+j] = (c[i+j] + x*y) % q
    return c
def ppow(a, e, m):
    r = [1]; a = pmod(a, m)
    while e:
        if e & 1: r = pmod(pmul(r, a), m)
        a = pmod(pmul(a, a), m); e >>= 1
    return r
def monics(k):
    for t in itertools.product(range(q), repeat=k): yield list(t) + [1]
def is_irr(p):
    k = len(p)-1
    for a in range(1, k//2+1):
        for f in monics(a):
            if not pmod(p, f): return False
    return True
irr = {k: [p for p in monics(k) if is_irr(p)] for k in range(1, 7)}
def factor(n):
    fs = []; n = n[:]
    for k in range(1, len(n)):
        for P in irr.get(k, []):
            while len(n) > 1 and not pmod(n, P):
                # divide
                quo = []; a = n[:]
                out = [0]*(len(a)-len(P)+1)
                for i in range(len(a)-len(P), -1, -1):
                    c = a[i+len(P)-1]; out[i] = c
                    for j in range(len(P)): a[i+j] = (a[i+j] - c*P[j]) % q
                n = out; fs.append(tuple(P))
    return fs
g = next(g for g in range(2, q) if len({pow(g, i, q) for i in range(q-1)}) == q-1)
w0 = pow(g, (q-1)//6, q); logw = {pow(w0, i, q): i for i in range(6)}
def sym(a, P):  # (a/P)_6 exponent or None
    r = ppow(a, (q**(len(P)-1)-1)//6, P)
    if not r: return None
    assert len(r) == 1
    return logw[r[0]]
import random
random.seed(1)
us = [list(random.choice(list(itertools.product(range(q), repeat=d)))) + [1] for _ in range(6)]
for u in us:
    fu = factor(u)
    def chi(n):
        e = 0
        for P in fu:
            s = sym(n, list(P))
            if s is None: return 0
            e += s
        return np.exp(2j*np.pi*e/6)
    # L coefficients
    a = [sum(chi(n) for n in monics(k)) for k in range(d+1)]
    M = []
    for k in range(0, d+2):
        tot = 0
        for n in monics(k):
            fn = factor(n)
            if len(set(fn)) < len(fn): continue
            tot += (-1)**len(fn) * chi(n)
        M.append(tot)
    roots = np.roots(a[:d][::-1]) if d > 1 else []
    print("u", u, "factors", fu, "a_d≈", abs(a[d]).round(6), "|roots|*sqrt(q):", np.round(np.abs(roots)*np.sqrt(q), 4), "|M_k|^2/q^k:", [round(abs(m)**2/q**k, 3) for k, m in enumerate(M)])
