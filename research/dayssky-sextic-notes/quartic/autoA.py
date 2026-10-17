# Approach (A) numerics.
# a(c) = alpha(c) * gtilde(conj chi_c) = mu(c) eps(c) gtilde(chi_c)   (identity I1)
# Tests:
#  (P1) alpha(c) gtilde(conj chi_c)^2 == mu(c) eps(c) chi_c(-1)
#  (P2) conj chi_{c1}(c2) chi_{c2}(c1) == (-1)^{((Nc1-1)/4)((Nc2-1)/4)} (quartic reciprocity, composite)
#  (N3) conductor detection: S_m(Z) = sum_{c sqfree primary,(c,m)=1} a(c) conj chi_c(m) W(Nc/Z)
#       normalised by sqrt(sum W^2); find Z beyond which it collapses, compare with N(m)^2 vs N(m)^3.
import sys, math, time, cmath
import numpy as np
sys.path.insert(0, '.')
from gauss4 import (gmul, gnorm, gconj, gmod, chi_prime, chi, gfactor, gauss,
                    mobius, is_primary, primary_assoc)
from dualms import primes_upto, prim_root

def eps_of(n):
    return 1 if (n[0] % 4, n[1] % 4) == (1, 0) else -1

def prime_list(Zmax):
    """primary primes with norm <= Zmax and their normalised Gauss sums gtilde(chi_pi)"""
    P = []
    for p in primes_upto(Zmax):
        p = int(p)
        if p == 2: continue
        if p % 4 == 1:
            a = next(a for a in range(1, int(p**.5)+1) if int(round((p-a*a)**.5))**2 == p-a*a)
            b = int(round((p-a*a)**.5))
            g = prim_root(p)
            ind = np.zeros(p, np.int64); x = 1
            for t in range(p-1): ind[x] = t; x = x*g % p
            for pi in (primary_assoc((a, b)), primary_assoc((a, -b))):
                A, B = pi
                r = (-A*pow(B, -1, p)) % p      # image of i in F_p (pi | A + B i)
                z = pow(g, (p-1)//4, p)
                e = 1 if z == r else (3 if z == (p-r) % p else None)
                xs = np.arange(1, p)
                k = (e*ind[xs]) % 4
                # psi(x/pi) for rational x: Im(x conj(pi))/p = -B x/p
                gs = np.sum((1j**k)*np.exp(-2j*np.pi*B*xs/p))
                P.append((pi, p, gs/math.sqrt(p)))
        else:
            if p*p > Zmax: continue
            pi = primary_assoc((p, 0))
            gs = gauss(pi, 1)
            P.append((pi, p*p, gs/p))
    P.sort(key=lambda t: t[1])
    return P

def build_c(Zmax, P, ms):
    """enumerate squarefree primary c with N(c)<=Zmax; return norms, a(c), chi_c(m) for m in ms"""
    out_N = []; out_a = []; out_chi = {m: [] for m in ms}; out_c = []
    norms = [t[1] for t in P]
    chim_p = {m: [chi_prime(m, t[0]) for t in P] for m in ms}   # chi_pi(m)
    def rec(start, c, N, gt, mu, chims, plist):
        if N > 1:
            out_N.append(N); out_c.append(c)
            out_a.append(mu*eps_of(c)*gt)
            for m in ms: out_chi[m].append(chims[m])
        for j in range(start, len(P)):
            pi, Np, gp = P[j]
            if N*Np > Zmax: break
            # cross symbols chi_c(pi) * chi_pi(c)
            cross = 1
            for q in plist: cross *= chi_prime(pi, P[q][0])
            cross *= chi_prime(c, pi) if N > 1 else 1
            newc = gmul(c, pi)
            rec(j+1, newc, N*Np, gt*gp*cross, -mu,
                {m: chims[m]*chim_p[m][j] for m in ms}, plist+[j])
    rec(0, (1, 0), 1, 1+0j, 1, {m: 1 for m in ms}, [])
    return (np.array(out_N, float), np.array(out_a), {m: np.array(v) for m, v in out_chi.items()}, out_c)

def W(t):
    out = np.zeros_like(t)
    m = (t > 0.5) & (t < 1.0)
    u = t[m]
    out[m] = np.exp(-1.0/((u-0.5)*(1.0-u))/4.0)
    return out

if __name__ == "__main__":
    Zmax = int(sys.argv[1]) if len(sys.argv) > 1 else 20000
    t0 = time.time()
    P = prime_list(Zmax)
    print(f"primes: {len(P)}  ({time.time()-t0:.1f}s)", flush=True)
    # test m's (primary): primes of norm 5,13,17,29,37,41, inert 3 (=-3 primary), squares, composites
    def pp(p, k=0): return primary_prime_over_k(p, k)
    from gauss4 import primary_prime_over
    pr = {p: primary_prime_over(p)[0] for p in (5, 13, 17, 29, 37, 41, 53)}
    m3 = primary_assoc((3, 0))
    ms = {}
    ms['1'] = (1, 0)
    for p, pi in pr.items(): ms[f'pi{p}'] = pi
    ms['q3'] = m3
    for p in (5, 13, 17):
        ms[f'pi{p}^2'] = gmul(pr[p], pr[p])
    ms['q3^2'] = gmul(m3, m3)
    ms['pi5*pi5b'] = gmul(primary_prime_over(5)[0], primary_prime_over(5)[1])
    ms['pi5*pi13'] = gmul(pr[5], pr[13])
    ms['pi5^3'] = gmul(gmul(pr[5], pr[5]), pr[5])
    mlist = list(set(ms.values()))
    Ns, A, CH, cs = build_c(Zmax, P, mlist)
    print(f"squarefree c: {len(Ns)}  ({time.time()-t0:.1f}s)", flush=True)
    np.savez('/tmp/claude-0/-home-user-QuasiRiemannTracker/726b9ed3-631a-533e-b4f5-2a9b6d3bc0e0/scratchpad/quartic/autoA_data.npz',
             Ns=Ns, A=A, **{f'chi_{k}': CH[v] for k, v in ms.items()})
    # (P1)/(I1) spot check against direct Gauss sums for small c
    bad = 0; cnt = 0; badP1 = 0
    for c, Nc, a in zip(cs, Ns, A):
        if Nc > 600: continue
        fac = gfactor(c)
        al = complex(*c)/math.sqrt(Nc)
        gconjt = gauss(c, -1, fac=fac)/math.sqrt(Nc)
        cnt += 1
        if abs(al*gconjt - a) > 1e-8: bad += 1
        chim1 = chi((-1, 0), c, fac)
        if abs(al*gconjt**2 - mobius(fac)*eps_of(c)*chim1) > 1e-8: badP1 += 1
    print(f"a(c) via multiplicativity vs direct alpha*gtilde(conj chi): {cnt} checked, {bad} mismatches; (P1) mismatches {badP1}")
    # conductor detection
    Zs = [int(x) for x in np.geomspace(20, Zmax, 26)]
    print("\n|S_m(Z)|/sqrt(sum W^2):  rows = m, columns = Z")
    print(f"{'m':>10} {'N(m)':>6} " + " ".join(f"{z:>7d}" for z in Zs))
    for name, m in ms.items():
        ch = CH[m]
        line = f"{name:>10} {gnorm(m):>6} "
        for Z in Zs:
            w = W(Ns/Z)
            S = np.sum(A*np.conj(ch)*w)
            nrm = math.sqrt(np.sum(w*w*(ch != 0)))
            line += f" {abs(S)/nrm:7.3f}" if nrm > 0 else "    nan"
        print(line, flush=True)
    print(f"total {time.time()-t0:.1f}s")
