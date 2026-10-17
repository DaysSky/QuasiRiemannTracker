# Numerical test of the "no theta pole" claim for the alpha-twisted Kubota series.
#   S_j(r;Z) = sum_{c sqfree primary} alpha(c)^j * gtilde(conj chi_c, r) * W(N c / Z)
# with r = 1, m^2, m^3 (m = test primes in the dataset).  gtilde(conj chi_c, r) = chi_c(r) gtilde(conj chi_c).
# Pole at s=3/4 (normalised)  <=>  S ~ c * What(3/4) * Z^{3/4};  no pole  <=>  S = O(Z^{1/2+eps}).
# Blocks Z_i = Z0 * 2^i have DISJOINT supports (W supported in [1/2,1]) -> independent under a random model.
# GLS estimate of c, standard error from empirical residual variance, z^2 ~ chi^2_2 under H0 (no pole).
import sys, math, numpy as np
sys.path.insert(0, '../condB'); sys.path.insert(0, '..')
from ser import Data, W

s = Data(sys.argv[1] if len(sys.argv) > 1 else '../condB/t2e6.bin')
Nf = s.Nf
xs = np.linspace(0.5, 1, 200001); Wx = W(xs)
What = lambda u: np.trapezoid(Wx * xs ** (u - 1), xs)
order = np.argsort(Nf); Nf_s = Nf[order]

def blocks(Z0):
    Zs = []; Z = Z0
    while Z <= s.Zmax: Zs.append(Z); Z *= 2
    return np.array(Zs, float)

def smooth(cf, Zs):
    out = np.zeros(len(Zs), complex); sig = np.zeros(len(Zs))
    cfs = cf[order]
    for i, Z in enumerate(Zs):
        lo = np.searchsorted(Nf_s, Z / 2); hi = np.searchsorted(Nf_s, Z, 'right')
        w = W(Nf_s[lo:hi] / Z)
        out[i] = np.sum(cfs[lo:hi] * w); sig[i] = math.sqrt(np.sum(np.abs(cfs[lo:hi]) ** 2 * w ** 2))
    return out, sig

def fit(cf, Z0s=(1000., 1414.)):
    Sall = []; Xall = []; sgall = []
    for Z0 in Z0s:
        Zs = blocks(Z0); S, sg = smooth(cf, Zs)
        Sall += list(S); Xall += list(What(0.75) * Zs ** 0.75); sgall += list(sg)
    S = np.array(Sall); X = np.array(Xall); sg = np.array(sgall)
    wgt = X / sg ** 2
    c = np.sum(wgt * S) / np.sum(wgt * X)
    res = (S - c * X) / sg
    kappa = np.mean(np.abs(res) ** 2)            # empirical noise variance relative to random-phase model
    se = math.sqrt(kappa / np.sum(X ** 2 / sg ** 2))
    ratio_top = np.abs(S[-1]) / sg[-1]
    return c, se, abs(c) / se, kappa, np.max(np.abs(S) / sg)

print("data: n=%d  Zmax=%d  K=%d test primes" % (s.n, s.Zmax, s.K))
print("%-28s %10s %9s %7s %7s %9s" % ("series", "|c|", "SE", "z", "kappa", "max|S|/sig"))
G = s.Gt
for j in (0, 1, -1, 2):
    a = s.alpha ** j
    c, se, z, ka, mx = fit(a * G)
    print("%-28s %10.5f %9.5f %7.2f %7.2f %9.2f" % (f"alpha^{j}, r=1", abs(c), se, z, ka, mx))
# indices m^2 and m^3 for the test primes
for (lab, pw) in (("m^3", 3), ("m^2", 2), ("m", 1)):
    for j in (0, 1):
        zs = []; cs = []
        for k in range(s.K):
            cf = (s.alpha ** j) * G * s.chipow(k, pw)
            c, se, z, ka, mx = fit(cf)
            zs.append(z); cs.append(abs(c))
        zs = np.array(zs)
        print("%-28s mean z^2=%6.2f  median z=%5.2f  max z=%6.2f  #(z>3)=%2d of %d   mean|c|=%.4f" %
              (f"alpha^{j}, r={lab}", np.mean(zs ** 2), np.median(zs), zs.max(), int(np.sum(zs > 3)), s.K, np.mean(cs)))
