# Marked-moment ratio at the binding Mobius length r* = 0.745 with slots beyond the paper's r+2z<=1, sampled rows at large U.
import sys, os, math, numpy as np, subprocess, time
import mmdata
from mmdata import Ideals, W
from eis import mul, norm
U = float(sys.argv[1]); r = float(sys.argv[2]); zs = [float(x) for x in sys.argv[3].split(',')]; nsamp = int(sys.argv[4])
t0 = time.time()
N = U**r; Zmax = U**max(zs)
I = Ideals(int(max(N, Zmax)) + 2)
print(f"ideals {len(I.N)} build {time.time()-t0:.0f}s", flush=True)
# sample rows u = x + y w with N(u) <= U, sixth-power-free w.r.t. small primes (checked against prime ideals with norm^6 <= U)
rng = np.random.default_rng(7)
B = int(math.isqrt(int(4*U/3))) + 3
xs, ys = [], []
small = [p for p in I.primes if p['P']**6 <= U]
while len(xs) < nsamp:
    x = int(rng.integers(-B, B+1)); y = int(rng.integers(-B, B+1))
    n = x*x - x*y + y*y
    if n == 0 or n > U: continue
    ok = True
    for p in small:
        g = p['g']; g6 = (1, 0)
        for _ in range(6): g6 = mul(g6, g)
        Ng = norm(g6); cg = (g6[0] - g6[1], -g6[1])
        pa = x*cg[0] - y*cg[1]; pb = x*cg[1] + y*cg[0] - y*cg[1]
        if pa % Ng == 0 and pb % Ng == 0: ok = False; break
    if ok: xs.append(x); ys.append(y)
rx = np.array(xs, np.int64); ry = np.array(ys, np.int64)
isprime = np.zeros(len(I.N), bool); isprime[1:] = (I.cof[1:] == 0)
specs = []; key = []
lo, hi = I.rng(N/2, N); specs.append((lo, hi, I.mu[lo:hi]*W(I.N[lo:hi]/N)))
for z in zs:
    Z = U**z; lo, hi = I.rng(Z/2, Z); base = isprime[lo:hi]*W(I.N[lo:hi]/Z)
    ph = np.exp(2j*np.pi*rng.random(hi-lo))
    for lab, a in (('one', np.ones(hi-lo)), ('rand', ph)):
        specs.append((lo, hi, (base*a).real)); specs.append((lo, hi, (base*a).imag)); key.append((z, lab, int(isprime[lo:hi].sum())))
def run(I, rx, ry, specs, tag):
    fin = os.path.join(mmdata.HERE, f"in_{tag}.bin"); fout = os.path.join(mmdata.HERE, f"out_{tag}.bin")
    np_ = len(I.primes); ni = len(I.N); nr = len(rx); ns = len(specs)
    kind = np.array([0 if p['kind'] == 'split' else 1 for p in I.primes], np.int64)
    P = np.array([p['P'] for p in I.primes], np.int64)
    par = np.array([p['r'] if p['kind'] == 'split' else p['q'] for p in I.primes], np.int64)
    with open(fin, 'wb') as f:
        np.array([np_, ni, nr, ns], np.int64).tofile(f); kind.tofile(f); P.tofile(f); par.tofile(f)
        I.spf.astype(np.int64).tofile(f); I.cof.astype(np.int64).tofile(f); rx.tofile(f); ry.tofile(f)
        np.array([s[0] for s in specs], np.int64).tofile(f); np.array([s[1] for s in specs], np.int64).tofile(f)
        for s in specs: np.asarray(s[2], np.float64).tofile(f)
    subprocess.run([os.path.join(mmdata.HERE, 'mmk2'), fin, fout], check=True)
    out = np.fromfile(fout, np.float32).reshape(nr, ns, 2); os.remove(fin); os.remove(fout)
    return out[..., 0].astype(np.float64) + 1j*out[..., 1].astype(np.float64)
V = run(I, rx, ry, specs, f"sb{int(U)}")
M = V[:, 0]; m2 = np.mean(np.abs(M)**2); c = 1
for (z, lab, npr) in key:
    Pp = V[:, c] + 1j*V[:, c+1]; c += 2
    p2 = np.mean(np.abs(Pp)**2); x = np.abs(M*Pp)**2
    rho = x.mean()/(m2*p2); se = x.std()/math.sqrt(len(x))/(m2*p2)
    print(f"U={U:.0e} r={r} z={z} ({npr} prime ideals) r+z={r+z:.3f} r+2z={r+2*z:.3f} a={lab}: rho = {rho:.4f} +- {se:.4f}", flush=True)
print(f"total {time.time()-t0:.0f}s")
