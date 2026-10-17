# LM numerics at larger D with table-free kernel mmk2 and row subsampling.
import sys, os, math, numpy as np, subprocess, time
import mmdata
from mmdata import Ideals, rows, W
U = int(float(sys.argv[1])); smax = float(sys.argv[2]); nsamp = int(sys.argv[3]); kern = sys.argv[4] if len(sys.argv) > 4 else 'mmk2'
t0 = time.time()
I = Ideals(int(U**smax) + 2); rx, ry = rows(U)
print(f"ideals {len(I.N)} rows {len(rx)} build {time.time()-t0:.0f}s", flush=True)
rng = np.random.default_rng(1); sel = rng.choice(len(rx), min(nsamp, len(rx)), replace=False); rx, ry = rx[sel], ry[sel]
sqf = I.mu != 0; fp = rng.choice([-1, 1], size=len(I.primes)); em = np.zeros(len(I.N), np.int64); em[0] = 1
for i in range(1, len(I.N)):
    if sqf[i]: em[i] = fp[I.spf[i]]*em[I.cof[i]]
SS = [s for s in (1.0, 1.1, 1.2, 1.3, 1.4, 1.5) if s <= smax + 1e-9]
specs = []; key = []
for s in SS:
    D = U**s; lo, hi = I.rng(D/2, D)
    for lab, c in (('mu', I.mu.astype(float)), ('rmul', em.astype(float))):
        w = c[lo:hi]*W(I.N[lo:hi]/D); specs.append((lo, hi, w)); key.append((s, lab, float((w**2).sum())))
orig = os.path.join(mmdata.HERE, 'mmk')
# swap kernel binary
def run(I, rx, ry, specs, tag):
    fin = os.path.join(mmdata.HERE, f"in_{tag}.bin"); fout = os.path.join(mmdata.HERE, f"out_{tag}.bin")
    np_ = len(I.primes); ni = len(I.N); nr = len(rx); ns = len(specs)
    kind = np.array([0 if p['kind'] == 'split' else 1 for p in I.primes], np.int64)
    P = np.array([p['P'] for p in I.primes], np.int64)
    par = np.array([p['r'] if p['kind'] == 'split' else p['q'] for p in I.primes], np.int64)
    with open(fin, 'wb') as f:
        np.array([np_, ni, nr, ns], np.int64).tofile(f); kind.tofile(f); P.tofile(f); par.tofile(f)
        I.spf.astype(np.int64).tofile(f); I.cof.astype(np.int64).tofile(f); rx.astype(np.int64).tofile(f); ry.astype(np.int64).tofile(f)
        np.array([s[0] for s in specs], np.int64).tofile(f); np.array([s[1] for s in specs], np.int64).tofile(f)
        for s in specs: np.asarray(s[2], np.float64).tofile(f)
    subprocess.run([os.path.join(mmdata.HERE, kern), fin, fout], check=True)
    out = np.fromfile(fout, np.float32).reshape(nr, ns, 2); os.remove(fin); os.remove(fout)
    return out[..., 0].astype(np.float64) + 1j*out[..., 1].astype(np.float64)
V = run(I, rx, ry, specs, f"lmb{U}")
for j, (s, lab, d) in enumerate(key):
    a2 = np.abs(V[:, j])**2
    se = a2.std()/math.sqrt(len(a2))/d
    print(f"U={U:.0e} s={s:.2f} {lab:4s} mean|M|^2/diag = {a2.mean()/d:.4f} +- {se:.4f}  max/mean = {a2.max()/a2.mean():.1f}", flush=True)
print(f"total {time.time()-t0:.0f}s")
