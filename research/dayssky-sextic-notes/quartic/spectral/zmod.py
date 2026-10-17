# Arithmetic of Z[i]/c in a vectorised model, for the metaplectic Kloosterman-sum experiments.
#  * c = a+bi with gcd(a,b)=1 (split primes, products of distinct split primes, squares of split primes):
#    Z[i]/c = Z/N, N=N(c), via i -> iota = -a*b^{-1} mod N.  psi(z/c)=e(Im(z conj c)/N) = e(-b*x/N) for x = image of z.
#  * inert p = 3 mod 4: Z[i]/p = F_{p^2}, elements (u,v) <-> u+v i, index u+p*v, psi(z/(-p)) = e(-v/p) (pi=-p primary).
# chi_c(x) = prod_k (x/pi_k)_4 (quartic residue symbol, depends only on the ideal), values in {0,1,i,-1,-i}.
import numpy as np, math, sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from gauss4 import gfactor, gnorm, primary_assoc, is_primary

I4 = np.array([1, 1j, -1, -1j])

class ZmodSplit:
    """Z[i]/c with gcd(Re c, Im c)=1, c odd primary (squarefree or square of a split prime)."""
    def __init__(self, c):
        self.c = c; a, b = c; N = a*a+b*b; self.N = N; self.b = b
        assert math.gcd(a, b) == 1 and N % 2 == 1
        self.iota = (-a*pow(b, -1, N)) % N
        x = np.arange(N, dtype=np.int64)
        self.fac = gfactor(c)
        chi = np.ones(N, complex)
        for pi in self.fac:
            p = gnorm(pi); io = (-pi[0]*pow(pi[1], -1, p)) % p
            assert (pi[0]+pi[1]*io) % p == 0
            pw = [pow(io, k, p) for k in range(4)]   # images of 1,i,-1,-i mod pi
            v = np.array([pow(int(t), (p-1)//4, p) for t in x % p])
            lut = np.zeros(N, complex)
            for k in range(4): lut[v == pw[k]] = I4[k]
            lut[(x % p) == 0] = 0
            chi *= lut
        self.chi = chi
        self.units = np.array([t for t in range(N) if math.gcd(t, N) == 1], dtype=np.int64)
        self.inv = np.array([pow(int(t), -1, N) for t in self.units], dtype=np.int64)
    def e(self, x):            # psi(x/c) for residues x (int array)
        return np.exp(-2j*np.pi*((self.b*(x % self.N)) % self.N)/self.N)
    def kloost(self, r, rp):   # S(r,r';c) = sum_{d unit} chi_c(d) psi((r d + r' dbar)/c)
        d = self.units; di = self.inv; N = self.N
        return np.sum(self.chi[d]*self.e((r*d + rp*di) % N))
    def kloost_cubes(self):    # vector over all x mod N of S(x^3,x^3;c)
        N = self.N; d = self.units; di = self.inv; s = (d+di) % N; ch = self.chi[d]
        x = np.arange(N, dtype=np.int64); x3 = (x*x % N)*x % N
        out = np.empty(N, complex)
        for k in range(N):
            out[k] = np.sum(ch*self.e((x3[k]*s) % N))
        return out
    def dual(self, Sx):        # Ssum(j) = sum_x Sx[x] psi(j x / c), all j mod N (DFT)
        N = self.N
        # psi(jx/c) = exp(-2 pi i b j x / N); j -> b j is a bijection mod N since gcd(b,N)=1
        F = np.fft.fft(Sx)        # F[k] = sum_x Sx[x] exp(-2 pi i k x/N)
        j = np.arange(N); return F[(self.b*j) % N]

class FieldInert:
    """F_{p^2} = Z[i]/p, p = 3 mod 4, pi = -p primary."""
    def __init__(self, p):
        assert p % 4 == 3
        self.p = p; q = p*p; self.N = q
        u = np.arange(q) % p; v = np.arange(q)//p
        self.u, self.v = u, v
        # chi(x) = x^{(q-1)/4} in F_{p^2}; compute by repeated squaring on (u,v)
        e = (q-1)//4
        ru = np.ones(q, np.int64); rv = np.zeros(q, np.int64); bu = u.astype(np.int64); bv = v.astype(np.int64)
        while e:
            if e & 1:
                ru, rv = (ru*bu - rv*bv) % p, (ru*bv + rv*bu) % p
            bu, bv = (bu*bu - bv*bv) % p, (2*bu*bv) % p
            e >>= 1
        chi = np.zeros(q, complex)
        chi[(ru == 1) & (rv == 0)] = 1; chi[(ru == 0) & (rv == 1)] = 1j
        chi[(ru == p-1) & (rv == 0)] = -1; chi[(ru == 0) & (rv == p-1)] = -1j
        chi[0] = 0
        self.chi = chi
        self.units = np.arange(1, q)
        # inverse: x^{-1} = conj(x)/N(x)
        nn = (u*u + v*v) % p
        ninv = np.array([pow(int(t), -1, p) if t else 0 for t in nn])
        self.invu = (u*ninv) % p; self.invv = ((-v)*ninv) % p
    def mul(self, x, y):
        p = self.p; xu, xv = x % p, x//p; yu, yv = y % p, y//p
        return ((xu*yu - xv*yv) % p) + p*((xu*yv + xv*yu) % p)
    def e(self, x):  # psi(x/pi), pi=-p: Im(x*conj(pi))/p^2 = -v/p
        return np.exp(-2j*np.pi*(x//self.p)/self.p)
    def kloost(self, r, rp):
        d = self.units; di = self.invu[d] + self.p*self.invv[d]
        a = (self.mul(np.full_like(d, r), d)); bb = self.mul(np.full_like(d, rp), di)
        s = ((a % self.p + bb % self.p) % self.p) + self.p*((a//self.p + bb//self.p) % self.p)
        return np.sum(self.chi[d]*self.e(s))
    def add(self, x, y):
        p = self.p; return ((x % p + y % p) % p) + p*((x//p + y//p) % p)
    def kloost_cubes(self):
        q = self.N; d = self.units; di = self.invu[d] + self.p*self.invv[d]
        s = self.add(d, di); ch = self.chi[d]
        x = np.arange(q); x3 = self.mul(self.mul(x, x), x)
        out = np.empty(q, complex)
        for k in range(q):
            out[k] = np.sum(ch*self.e(self.mul(np.full_like(s, x3[k]), s)))
        return out
    def dual(self, Sx):
        q = self.N; x = np.arange(q); out = np.empty(q, complex)
        for j in range(q):
            out[j] = np.sum(Sx*self.e(self.mul(np.full_like(x, j), x)))
        return out

def split_primes(maxN, minN=5):
    out = []
    for p in range(minN, maxN+1):
        if p % 4 != 1 or any(p % t == 0 for t in range(2, int(p**0.5)+1)): continue
        for a in range(1, int(p**0.5)+1):
            b2 = p-a*a; b = int(round(b2**0.5))
            if b*b == b2:
                out.append(primary_assoc((a, b))); break
    return out
