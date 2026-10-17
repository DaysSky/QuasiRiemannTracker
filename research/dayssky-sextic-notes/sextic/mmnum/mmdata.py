# Data preparation for the sextic mixed-moment numerics.
import numpy as np, os, subprocess, struct, math
from eis import prime_ideals, mul, norm, primary_assoc, divides

HERE = os.path.dirname(os.path.abspath(__file__))

def W(t):
    """smooth bump supported on [1/2,1], max 1 at t=3/4."""
    t = np.asarray(t, float); s = 4*t - 3
    out = np.zeros_like(t); m = np.abs(s) < 1
    out[m] = np.exp(1 - 1/(1 - s[m]**2))
    return out

class Ideals:
    """odd primary ideals (generators = 1 mod 3, coprime to 2) with norm <= X, sorted by norm; index 0 = (1)."""
    def __init__(self, X):
        self.X = X
        self.primes = prime_ideals(X)
        B = int(math.isqrt(4*X//3)) + 3
        L = []
        for b in range(-B, B + 1):
            if b % 3 != 0: continue
            for a in range(-B, B + 1):
                if a % 3 != 1: continue
                n = a*a - a*b + b*b
                if 0 < n <= X and not (a % 2 == 0 and b % 2 == 0):
                    L.append((n, a, b))
        L.sort()
        self.N = np.array([t[0] for t in L], np.int64)
        self.g = [(t[1], t[2]) for t in L]
        assert self.g[0] == (1, 0)
        self.idx = {g: i for i, g in enumerate(self.g)}
        n = len(L)
        spf = -np.ones(n, np.int64); cof = np.zeros(n, np.int64)
        for j, pr in enumerate(self.primes):
            P = pr['P']; gp = pr['g']
            for i in range(n):
                if self.N[i] * P > X: break
                k = self.idx[mul(gp, self.g[i])]
                if spf[k] < 0:
                    spf[k] = j; cof[k] = i
        assert (spf[1:] >= 0).all()
        self.spf, self.cof = spf, cof
        mu = np.zeros(n, np.int64); mu[0] = 1
        for i in range(1, n):
            c = cof[i]
            mu[i] = 0 if (c > 0 and spf[c] == spf[i]) else -mu[c]
        self.mu = mu
        self.Pnorm = np.array([p['P'] for p in self.primes], np.int64)

    def factor(self, i):
        """list of prime indices with multiplicity."""
        out = []
        while i > 0:
            out.append(int(self.spf[i])); i = int(self.cof[i])
        return out

    def rng(self, lo, hi):
        """index range of ideals with lo <= N < hi... we use N in [lo, hi]."""
        a = int(np.searchsorted(self.N, lo, 'left')); b = int(np.searchsorted(self.N, hi, 'right'))
        return a, b

def rows(U, sixth_power_free=True):
    B = int(math.isqrt(4*U//3)) + 3
    y, x = np.mgrid[-B:B+1, -B:B+1]
    x = x.ravel().astype(np.int64); y = y.ravel().astype(np.int64)
    n = x*x - x*y + y*y
    m = (n > 0) & (n <= U)
    x, y = x[m], y[m]
    if sixth_power_free:
        bad = np.zeros(len(x), bool)
        # prime ideals pi with N(pi)^6 <= U: lambda (norm 3), 2 (norm 4), norm 7, 13 (13^6 = 4.8e6), ...
        divs = []
        for p in [3, 4, 7, 13, 19, 25, 31, 37]:
            if p**6 <= U:
                divs.append(p)
        gens = {3: [(1, 2)], 4: [(2, 0)]}
        from eis import prime_ideals as _pi
        for pr in _pi(40):
            gens.setdefault(pr['P'], []).append(pr['g'])
        for p in divs:
            for g in gens[p]:
                g6 = (1, 0)
                for _ in range(6): g6 = mul(g6, g)
                # u divisible by g6 <=> u * conj(g6) == 0 mod N(g6)
                Ng = norm(g6); cg = (g6[0] - g6[1], -g6[1])
                pa = x*cg[0] - y*cg[1]; pb = x*cg[1] + y*cg[0] - y*cg[1]
                bad |= (pa % Ng == 0) & (pb % Ng == 0)
        x, y = x[~bad], y[~bad]
    return x, y

def run_kernel(I, rx, ry, specs, tag):
    """specs: list of (lo_idx, hi_idx, weights array). returns complex array [nrows, nspecs]."""
    fin = os.path.join(HERE, f"in_{tag}.bin"); fout = os.path.join(HERE, f"out_{tag}.bin")
    np_ = len(I.primes); ni = len(I.N); nr = len(rx); ns = len(specs)
    kind = np.array([0 if p['kind'] == 'split' else 1 for p in I.primes], np.int64)
    P = np.array([p['P'] for p in I.primes], np.int64)
    par = np.array([p['r'] if p['kind'] == 'split' else p['q'] for p in I.primes], np.int64)
    with open(fin, 'wb') as f:
        np.array([np_, ni, nr, ns], np.int64).tofile(f)
        kind.tofile(f); P.tofile(f); par.tofile(f)
        I.spf.astype(np.int64).tofile(f); I.cof.astype(np.int64).tofile(f)
        rx.astype(np.int64).tofile(f); ry.astype(np.int64).tofile(f)
        np.array([s[0] for s in specs], np.int64).tofile(f)
        np.array([s[1] for s in specs], np.int64).tofile(f)
        for s in specs:
            np.asarray(s[2], np.float64).tofile(f)
    subprocess.run([os.path.join(HERE, 'mmk'), fin, fout], check=True)
    out = np.fromfile(fout, np.float32).reshape(nr, ns, 2)
    os.remove(fin); os.remove(fout)
    return out[..., 0].astype(np.float64) + 1j*out[..., 1].astype(np.float64)
