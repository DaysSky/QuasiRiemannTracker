# Function-field analogue of LM. Over F_q[t], q = 1 mod 6, rows u monic sixth-power-free of degree d,
# M_u(k) = sum_{deg n = k, n monic} mu(n) (n/u)_6  = T^k coefficient of 1/L(T, chi_u).
# L(T, chi_u) is a polynomial of degree <= d-1, obtained from chi_u(P), deg P < d, by Newton's identities.
# Output: MS(k)/(#u q^k) (U(N)-Haar value E|h_k|^2 = 1) and MS(k)/TD(k) (exact coprime squarefree diagonal).
import sys, numpy as np, time

def prim_root(q):
    for g in range(2, q):
        if all(pow(g, (q-1)//p, q) != 1 for p in set(f for f in range(2, q) if (q-1) % f == 0 and all(f % r for r in range(2, f)))):
            return g

def idx_to_coefs(N, m, q):
    a = np.arange(N, dtype=np.int64); out = np.zeros((N, m), dtype=np.int64)
    for j in range(m):
        out[:, j] = a % q; a //= q
    return out

def mulmod(A, B, P, q):
    # A, B: (N, m) residues mod monic P (len m+1, P[m] = 1)
    N, m = A.shape
    C = np.zeros((N, 2*m-1), dtype=np.int64)
    for i in range(m):
        C[:, i:i+m] += A[:, i:i+1] * B
    C %= q
    for deg in range(2*m-2, m-1, -1):
        c = C[:, deg].copy()
        C[:, deg-m:deg+1] -= c[:, None] * P[None, :]
        C %= q
    return C[:, :m]

def powmod(A, e, P, q):
    N, m = A.shape
    R = np.zeros_like(A); R[:, 0] = 1
    B = A.copy()
    while e:
        if e & 1: R = mulmod(R, B, P, q)
        e >>= 1
        if e: B = mulmod(B, B, P, q)
    return R

def monic_irreducibles(q, m):
    # all monic of degree m as coefficient arrays (len m+1); sieve out products
    N = q**m
    allc = idx_to_coefs(N, m, q)
    red = np.zeros(N, dtype=bool)
    for a in range(1, m//2 + 1):
        A = np.concatenate([idx_to_coefs(q**a, a, q), np.ones((q**a, 1), dtype=np.int64)], axis=1)
        b = m - a
        B = np.concatenate([idx_to_coefs(q**b, b, q), np.ones((q**b, 1), dtype=np.int64)], axis=1)
        # products: (q^a * q^b, m+1)
        prod = np.zeros((q**a, q**b, m+1), dtype=np.int64)
        for i in range(a+1):
            prod[:, :, i:i+b+1] += A[:, None, i:i+1] * B[None, :, :]
        prod %= q
        idx = np.zeros(prod.shape[:2], dtype=np.int64)
        for j in range(m-1, -1, -1):
            idx = idx * q + prod[:, :, j]
        red[idx.ravel()] = True
    irr = np.nonzero(~red)[0]
    P = np.concatenate([allc[irr], np.ones((len(irr), 1), dtype=np.int64)], axis=1)
    return P

def run(q, d, K):
    t0 = time.time()
    g = prim_root(q); w0 = pow(g, (q-1)//6, q)
    logw = {pow(w0, i, q): i for i in range(6)}
    lut = np.full(q, -1, dtype=np.int64)
    for v, i in logw.items(): lut[v] = i
    NU = q**d
    U = np.concatenate([idx_to_coefs(NU, d, q), np.ones((NU, 1), dtype=np.int64)], axis=1)  # monic deg d
    keep = np.ones(NU, dtype=bool)
    if d >= 6:
        # exclude sixth powers (t+a)^6 (d = 6 only)
        assert d == 6
        from math import comb
        for a in range(q):
            co = [comb(6, i) * pow(a, 6-i, q) % q for i in range(7)]
            ii = 0
            for j in range(5, -1, -1): ii = ii * q + co[j]
            keep[ii] = False
    zeta = np.exp(2j*np.pi*np.arange(6)/6)
    c = np.zeros((d, NU), dtype=complex)          # c_j, j = 1..d-1
    td = np.zeros((NU, K+1))
    # base series (1 - q T^2)/(1 - q T): s_0 = 1, s_1 = q, s_k = q^k - q^{k-1} (k >= 2)
    base = np.array([1.0, q] + [q**k - q**(k-1) for k in range(2, K+1)])[:K+1]
    td[:] = base
    hasfac = np.zeros(NU, dtype=bool)
    sgnexp = 3 if (((q-1)//6) * d) % 2 else 0      # reciprocity (-1)^{((q-1)/6) deg P deg u}, per unit of deg P
    nP = 0
    for m in range(1, d):
        Ps = monic_irreducibles(q, m)
        e = (q**m - 1) // 6
        res = idx_to_coefs(q**m, m, q)
        for P in Ps:
            nP += 1
            pw = powmod(res, e, P, q)
            assert (pw[1:, 1:] == 0).all()
            sym = lut[pw[:, 0]]; sym[0] = -1
            # reduce U mod P
            R = U.copy()
            for deg in range(d, m-1, -1):
                cc = R[:, deg].copy()
                R[:, deg-m:deg+1] -= cc[:, None] * P[None, :]
                R %= q
            ii = np.zeros(NU, dtype=np.int64)
            for j in range(m-1, -1, -1): ii = ii * q + R[:, j]
            s = sym[ii]
            div = s < 0
            ex = (s + sgnexp * m) % 6
            for r in range(1, (d-1)//m + 1):
                val = np.where(div, 0, zeta[(r * ex) % 6])
                c[m*r] += m * val
            if div.any():
                hasfac |= div
                rows = np.nonzero(div)[0]
                for k in range(m, K+1):
                    td[rows, k] -= td[rows, k-m]
    # irreducible u (no factor of degree < d): divide by (1 + T^d)
    rows = np.nonzero(~hasfac)[0]
    for k in range(d, K+1):
        td[rows, k] -= td[rows, k-d]
    # Newton: k a_k = sum_{j=1}^k c_j a_{k-j}
    a = np.zeros((d, NU), dtype=complex); a[0] = 1
    for k in range(1, d):
        a[k] = sum(c[j] * a[k-j] for j in range(1, k+1)) / k
    b = np.zeros((K+1, NU), dtype=complex); b[0] = 1
    for k in range(1, K+1):
        b[k] = -sum(a[j] * b[k-j] for j in range(1, min(k, d-1) + 1))
    nu = float(keep.sum())
    # rows whose L-polynomial has a repeated root (CM-type members with repeated Frobenius eigenvalues)
    rep = np.zeros(NU, dtype=bool)
    for i in np.nonzero(keep)[0]:
        co = a[:, i]
        while len(co) > 1 and abs(co[-1]) < 1e-9: co = co[:-1]
        if len(co) <= 2: continue
        r = np.roots(co[::-1]) * np.sqrt(q)
        dd = np.abs(r[:, None] - r[None, :]) + np.eye(len(r)) * 9
        rep[i] = dd.min() < 1e-4
    gen = keep & ~rep
    print(f"  rows with repeated L-root: {int(rep[keep].sum())} of {int(nu)}")
    out = []
    for k in range(K+1):
        MS = (np.abs(b[k][keep])**2).sum()
        TD = td[keep, k].sum()
        MSg = (np.abs(b[k][gen])**2).sum(); TDg = td[gen, k].sum()
        out.append((k, k/d, MS/(nu * float(q)**k), MS/TD, (np.abs(b[k][keep])**2 / float(q)**k).max(), MSg/TDg, (np.abs(b[k][gen])**2 / float(q)**k).max()))
    print(f"q={q} d={d} rows={nu} primes={nP} time={time.time()-t0:.1f}s")
    for k, s, r1, r2, mx, rg, mxg in out:
        print(f"  k={k:2d} s=k/d={s:5.2f}  MS/(#u q^k)={r1:.4f}  MS/truediag={r2:.4f}  max|M|^2/q^k={mx:.2f}  | distinct-root rows: MS/truediag={rg:.4f} max={mxg:.2f}")
    sys.stdout.flush()

if __name__ == "__main__":
    q, d, K = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
    run(q, d, K)
