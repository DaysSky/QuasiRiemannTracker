# Shared problem statement: the quartic lemma (L^3) / (G) / (PO)

## Setting
K = Q(i), O = Z[i], lambda = 1+i. An element n is *primary* if n ≡ 1 mod (2+2i) (i.e. (Re n, Im n) ≡ (1,0) or (3,2) mod 4).
For primary squarefree n coprime to 2: chi_n(x) = (x/n)_4 quartic residue symbol (zero on non-units), extended multiplicatively in n.
Additive character psi(z) = exp(2πi Im z) (Z[i] is self-dual for this character).
Gauss sums: g(chi, r; n) = sum_{x mod n} chi(x) psi(r x / n); normalised gtilde = g / sqrt(N n).
Note g(chi_n, r) = conj(chi_n(r)) g(chi_n) for (r,n)=1; g(conj chi_n, r) = chi_n(r) g(conj chi_n).
alpha(n) = n/|n| (Hecke character of infinity type).

## Verified identities (numerically, all squarefree primary n with N(n)<=1300, and 513 random (n,h))
(I1)  mu(n) * gtilde(conj chi_n) = eps(n) * conj(alpha(n)) * gtilde(chi_n),   eps(n)=+1 if n≡1 mod 4, -1 if n≡3+2i mod 4
(I2)  mu(n) * g(chi_n, h) = eps(n) * alpha(n) * (h/n)_2 * g(conj chi_n, h)
Proof of (I1) for primes: J(chi,chi) = -chi(-1) pi for primary pi (Ireland–Rosen 9.9.4), g(chi)^2 = J(chi,chi) g(chi^2),
normalised quadratic Gauss sum = (b/p) = ±1 determined by pi mod 4; composite n by twisted multiplicativity
g(mn) = chi_m(n) chi_n(m) g(m) g(n) (ratios cancel).

## Why this matters
OpenAI's quasi-Riemann proof (Re s > 7/8 for all Dirichlet L-functions; github.com/openai/math, preprint
"The Quasi-Riemann Hypothesis", Sept 30 2026) uses the SEXTIC family over Q(sqrt -3): Poisson in the row variable turns
mu(n) * (sextic Gauss sum) into conj(alpha(n)) * (cubic Gauss sum) * (fixed ray factor), and cubic Gauss sums are exactly the
Fourier coefficients of Patterson's cubic theta function, whose automorphy gives a "reflection" shortening twisted sums
to length N(h)^2/X with a QUADRATIC dual twist, handled by the quadratic large sieve.
Exponent bookkeeping (derived, checked against the paper at k=6): for a family of order k the method gives
sigma0 = 1 - 1/(2k) - 3 ell/(2k) (ell = compensation length, ell < 1/5). k=6: 11/12 - ell/4 (current best ~0.874861).
k=4 (quartic family over Q(i)): 7/8 - 3 ell/8, i.e. 13/16 at ell=1/6; more generally sigma0 = 1 - M/4 where M>=2/3 is the
row-range exponent for which the lemma below holds. Any M>1/2 with everything else working beats 7/8; M=2/3 gives 5/6.

## THE MISSING LEMMA (three equivalent-ish forms)
(G)  For Q = Z^M with 2/3 <= M < 1 (smooth W supported in [1/2,1]):
     sum_{0<N(m)<=Q} | sum_{c sqfree primary} alpha(c) gtilde(conj chi_c) conj(chi_c(m)) W(N c / Z) |^2  <<  Q Z^{1+eps}.
     Note: the row sum has length ~Z > Q, so this needs a genuine "reflection"; the large-sieve shape (Q+Z)Z is NOT enough.
(L^3) Equivalently (via Mellin, no poles right of 1/2 for the alpha-twisted series):
     sum_{N(m)<=Q} |psi_alpha(m^3, 1/2+it)|^2 << Q^{1+eps} (1+|t|)^A, where psi_alpha(r,s) = sum_c alpha(c) gtilde(conj chi_c, r) N(c)^{-s}
     is (up to normalisation) the r-th Fourier-coefficient Dirichlet series of Kubota's quartic Eisenstein series (with K-type for alpha).
     Each psi(m^3, .) has conductor ~N(m)^3 by Kubota's functional equation (index m^3 is 4th-power-free): a THIN family
     (Q members, conductor Q^3), so approximate functional equation + any large sieve gives only Q^{3/2}.
(PO) Poisson-dual form: G at M  <=>  perfect orthogonality for Mobius coefficients:
     sum_{0<N(k)<=H} | sum_{N(c)~D} mu(c) chi_c(k) W(N c/D) |^2 << (H + D) D^{1+eps}  for H >= D^{2-M} (so H >= D^{4/3} at M=2/3).
     Heuristically equivalent to zero density N(sigma) << H^{1 - 3 delta/4 + eps} (delta = 2 sigma - 1) for the quartic Hecke family
     {L(s, chi_k): N(k) <= H} over Q(i).

## What is known / proved so far
- Quartic large sieve (Gao–Zhao): sum_{m<=M} |sum_{n<=N} a_n (n/m)_4|^2 << (M + N + (MN)^{3/4})(MN)^eps ||a||^2. Gives only
  sum |psi(m^3,1/2)|^2 << Q^{15/8}.
- de Faveri–Dunn–Hoffstein (arXiv 2607.07911, 2026): quartic and cubic large sieves are NOT perfectly orthogonal; obstruction = Gauss-sum
  bias; they prove Lindelöf on average over ALL indices r for Kubota's series (gives only Q^3 for cube indices).
- David–Dunn–Hamieh–Lin (arXiv 2306.11875): quartic Gauss sums over primes; use twisted multiplicativity g(mn) = (m/n)_2 * sign * g(m)g(n)
  to reduce bilinear forms to the QUADRATIC large sieve over Q(i); Suzuki evaluated quartic theta coefficients only at squares.
- Quartic theta coefficients at squarefree indices are unknown (Suzuki, Eckhardt–Patterson). Friedberg–Ginzburg descent identification conjectural.
- Numerics (X up to 80000, Q up to 3000): alpha-twisted rows behave exactly like random phases (no bias rows); untwisted show bias at h = 1, -4, 16, ...
- OpenAI's zero-free region covers Hecke L-functions over Q(sqrt -3) and Dirichlet L-functions, NOT general Hecke L-functions over Q(i).
- arXiv and most paper mirrors are blocked from this container (only github, raw.githubusercontent, pypi etc. reachable); WebSearch gives snippets.

## Code
gauss4.py (Z[i] arithmetic, quartic symbols, Gauss sums), dualms.py/run2.py (mean-square experiment), test_identity.py, test_dual.py.
