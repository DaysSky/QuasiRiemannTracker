# Check g6(h, n m) = g6(h, m) * sigma(n,m) * g6(h m^4, n), sigma a 6th root of unity independent of h  (n, m coprime squarefree, split primes)
import sys, math, numpy as np, itertools, random
exec(open('g2per.py').read().split('rows = []')[0].replace("P = int(sys.argv[1]) if len(sys.argv) > 1 else 400", "P = 200"))
def gh(hval, n, N, facs):
    # g6(h, n) = sum_{d mod N} chi_n(d) e(Tr(h d / n)/...) using the same additive convention as gauss(): e(d*T/N), T from n; parameter h multiplies d via its image mod N
    d = np.arange(N); ch = np.ones(N, complex)
    for (pi, p, wr, j0, ind) in facs:
        r = d % p
        ch *= np.where(r == 0, 0, W6[(j0*(ind[r] % 6)) % 6])
    T = (2*n[0] - n[1]) % N
    # image of h in Z/N: h = a + b*w, w -> root of x^2+x+1 mod N consistent with n; use per-prime CRT of w-images
    a, b = hval
    hw = 0; M = 1
    img = []
    for (pi, p, wr, j0, ind) in facs:
        img.append(((a + b*wr) % p, p))
    # CRT
    x = 0; mod = 1
    for (r, p) in img:
        t = ((r - x) * pow(mod, -1, p)) % p
        x += mod*t; mod *= p
    return np.sum(ch*np.exp(2j*np.pi*(((x*d) % N)*T % N)/N))/math.sqrt(N)
random.seed(5)
primes = PR
ok = 0; bad = 0; worst = 0
for trial in range(60):
    i, j = random.sample(range(len(primes)), 2)
    if primes[i][1] == primes[j][1]: continue
    m = primes[i]; nf = [primes[j]]
    n_el = nf[0][0]; N_n = nf[0][1]; N_m = m[1]
    nm = mul(n_el, m[0]); Nnm = N_n*N_m
    ratios = []
    for t in range(4):
        h = (random.randint(-30, 30), random.randint(-30, 30))
        if (h[0]*h[0] - h[0]*h[1] + h[1]*h[1]) % N_n == 0 or (h[0]*h[0] - h[0]*h[1] + h[1]*h[1]) % N_m == 0: continue
        m4 = mul(mul(m[0], m[0]), mul(m[0], m[0])); hm4 = mul(h, m4)
        lhs = gh(h, nm, Nnm, [nf[0], m]); r = lhs/(gh(h, m[0], N_m, [m]) * gh(hm4, n_el, N_n, nf))
        ratios.append(r)
    if len(ratios) < 2: continue
    ph = [round((np.angle(r)/(2*np.pi))*6, 6) % 6 for r in ratios]; mods = [abs(r) for r in ratios]
    const = max(ph) - min(ph) < 1e-6 and all(abs(x - round(x)) < 1e-6 for x in ph) and all(abs(mm - 1) < 1e-6 for mm in mods)
    ok += const; bad += (not const); worst = max(worst, max(abs(mm-1) for mm in mods))
print("pairs tested", ok + bad, "ratio a sixth root of unity independent of h:", ok, " failures:", bad, " max||ratio|-1|:", worst)
