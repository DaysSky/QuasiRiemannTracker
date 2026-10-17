# Quartic route for the zero-free half-plane: status notes

## 1. Identity (verified numerically, 452 squarefree primary n, N(n) <= 1300)
For squarefree primary n in Z[i] (n ≡ 1 mod (2+2i)), quartic symbol chi_n = (./n)_4,
psi(z) = exp(2πi Im z), gamma_j(n) = N(n)^{-1/2} sum_x chi_n(x)^j psi(x/n), alpha(n)=n/|n|:

    mu(n) gamma_{-1}(n) = eps(n) conj(alpha(n)) gamma_1(n),   eps(n) = +1 (n≡1 mod 4), -1 (n≡3+2i mod 4)

Proof sketch: J(chi,chi) = -chi(-1) pi for primary pi; g(chi)^2 = J g(chi^2); quadratic Gauss sum
normalised = (b/p) = ±1 by pi mod 4; twisted multiplicativity factors cancel in the ratio.

## 2. Dual coefficient after Poisson
mu(n) g(chi_n, h) = eps alpha(n) chi_n(h)^{-1}... = eps alpha(n) * g(conj chi_n, h) * (h/n)_2
=> quartic Gauss sum at index h times a QUADRATIC twist (h/n)_2  (sextic analogue: cubic Gauss sum
at index h times quadratic chi^3).

## 3. Exponent bookkeeping (order-k family, OpenAI Part II architecture)
C(s) = s - 1 + x/2 + h/k ; low = x/2 + b/(2k) (Gram: exceptional frequencies = k-th powers, count K^{1/k});
h = 1 - x + ell, x + y + ell = 1, b = y - x
=> sigma0 = 1 - 1/(2k) - 3 ell/(2k).     k=6: 11/12 - ell/4 (matches OpenAI/Cycle61). k=4: 7/8 - 3 ell/8.
delta=0 generic-row constraint: x(1/k - 1/2) + ell(1 + 1/(2k)) - 1/(2k) <= 0, with x <= (1-ell)/2  => ell <= 1/5 for k=4 and k=6.
Ceilings: sextic 13/15 ; quartic 4/5.  At OpenAI's ell = 1/6: quartic 13/16 = 0.8125.
Without compensation (ell=0) quartic gives exactly 7/8 (no gain).
High side at ell=1/6, x=17/48: same inequality shape as sextic Part II, C0 = -5/192 (sextic -4/192).

## 4. Missing input (unproved)
Quartic completed-row energy: for M in the needed range,
  sum_{N(m) <= Z^M} | Z^{-1/2} sum_c eta(c) alpha(c) gtilde(conj chi_c, m) (m/c)_2 W(N c / Z) |^2  << Z^{M+eps}.
Sextic version is proved via Patterson's cubic theta (Gauss sums = theta coefficients). No GL2 form has
quartic Gauss sums as coefficients (quartic theta coefficients known only at squares, Suzuki).
Plain quartic large sieve loses ~Z^{13/24}; budget for a lossy version is ~Z^{1/8}.
Also needed: quartic analogues of marked inverse moment, plain moments with prime slots, density counts.

## 5. Numerics (dualms.py, run2.py)
Mean square / (#rows * sum W^2): alpha-twisted (k=1) ≈ random (0.90 vs 0.90) at X=20000, 80000, H<=3000.
Untwisted (k=0) shows bias rows at h = 1, -4, 16, -(4+3i)^2, ... growing ~X^{1/2}.

## 6. Structural results from the (L^3) attack (my own + subagents A, B)
- Conductor of B_m is N(m)^3 (B): (m/c)_2 = chi_c(m)^2 is a Kummer twist for the 4-fold cover; Gamma_0(m)-nebentypus route returns index m^3.
- Reflection is the involution M <-> M/(3M-1), fixed point M = 2/3 = (L^3). Reflection + optimal dual large sieve proves (G) only for M <= 1/2
  (Gao–Zhao: M <= 3/8). Needed: M >= 2/3 (delta=0 constraint 5 ell + 3 g <= 1 with M = 1 - ell - g).
- Without beating the large-sieve shape (Q+Z) there is no gain: sigma0 = 7/8 + ell/8 + 3g/8 >= 7/8.
- Uniqueness of the sextic: for an order-k family, mu(n) g(chi^a, h) ∝ lambda(n) g(chi^{a+k/2}, h) (h/n)_2. The quadratic twist is Kummer
  (thin rows) iff ord(chi^{a+k/2}) is even. Usable needs that order k' odd with known theta: k'=3 -> k_eff=6 (sextic). Quartic: k'=4 even => thin.
- (L^3) => |psi(m^3,1/2)| << N(m)^{1/2+eps} = cond^{1/6}: Weyl-strength subconvexity for each member (convexity cond^{1/4}).
  Lossy version sufficient for < 0.8749 still implies cond-exponent ~0.19. Heuristic: rows ~ L(s, sym^2(theta_4) ⊗ chi_m), GL3 x quartic twist.
- Sextic architecture ceiling (with slack g): sigma0 = 11/12 - ell/4 - g/12 s.t. 5 ell + 3g <= 1 -> 13/15.

## 7. Subagent C (direct large sieve) and D (redesign) results
- Exact identity (D, 69/69): gtilde(chi_c)^2 = mu(c) alpha(c) for squarefree primary c.
- Correction (C, D): (G) <=> off-diagonal of the Mobius mean square at H = Z^{2-M} is O(D^{2+eps}) (asymptotic with relative saving),
  STRONGER than the large-sieve-shape (PO). Loss for (G) from any (H+D)D^{1+lambda} bound is Z^{1-M+lambda} >= Z^{1/3} at M=2/3.
- C: recursion (Poisson + quartic/quadratic/multiplicative large sieves by square class) proves the PO-shape with lambda = 1/6 at H=D^{4/3}
  (any |a_c|<=1), 1/9 conditional on uniform convexity for Kubota sums, floor 1/12 for these methods.
  Barrier: at M=5/6, lambda<1/8 would force zero-free strips Re s > 15/16 + lambda/2 for weight-1 dihedral Hecke L-functions over Q(i).
  Numerics (D up to 1e5, H=D^{4/3}, D^{7/6}): mu shows no off-diagonal excess (ratio 1.000 ± 0.004); g-tilde calibration shows predicted bias.
- D: no design in the natural class works (quadratic twists circular, Stickelberger no-go for importing sextic data, norm family loses the
  w=1 pole, WMDS has no mu, triple series holomorphy equivalent to zero-free region). Side remark: sextic argument over Q(zeta_12) would give
  7/8 for Hecke L-functions over Q(i) (needs cubic theta over that field).

## 8. Second wave (E spectral, F moments, G barrier) + sextic mixed-moment lead
- E: thin Kuznetsov is a reformulation; Weil gives Q^{5/2}; natural barrier Q^2 (sqrt cancellation in m only). Exact identity
  S(x^3,x^3;pi) = chi(x) Kl_chi(x^6); complete sums are rank-6 hypergeometric (moments 1,3,15 Gaussian). Needed Lemma N equivalent to (L^3)+cusp analogue.
- F: second-moment MDS group is hyperbolic Coxeter (2,6,6); Hecke reflection maps to Mobius/1/L side with no FE => only (Z/2)^2 acts;
  hull reaches centre only for Re w > 5/2 (convexity). Under GRH these methods give exactly Q^{3/2}. (L^3) sits at the fixed point of every duality.
  Missing lemma (ML_theta'): ratio-type Mobius mean square with error H^{theta'}; theta'=3/4 <=> (L^3); 13/16 <=> application budget.
- G: Theorem Z (PROVED): PO(eta,lambda) => L(s,(k0/.)_4) != 0 for Re s > lambda/2 + max(1-eta/8, 1/2+3eta/8) for all k0 (incl. zeta_Q(i)).
  Corollary: (L^3) uniform in Q with loss Q^{theta'} => zeta and all quartic Hecke L(s,(k/.)_4) zero-free for Re s > 7/8 + theta'/2.
  So (L^3) implies a new quasi-RH half-plane for the quartic Hecke family over Q(i) AND Weyl-strength subconvexity. Method barrier B1: AFE + any
  arbitrary-coefficient large sieve >= Q^{3/2-eps}. alpha-twisted series has no theta pole (proved modulo Langlands bookkeeping; numerics z<=0.44).
- SEXTIC LEAD (my model, sextic/model.py, fastmodel.py, mixed.py, mixed2.py; calibrated: OpenAI inputs -> theta_max 0.87504):
  single-input improvements (plain capacity, inverse amplification, inverse capacity) give at most 0.8702 even combined.
  Mixed moment (MM): sum_{N(u)<=U} |M_u(U^r) S_u(U^m)|^2 << U^{1+eps} for r+m<=1, 0<=m<=1/2 gives ell=1/5, theta=13/15=0.86667 (ceiling).
  Needs full range r+m=1 (r+m<=0.9 gives nothing; m<=0.4 gives 0.8746). Loss tolerance: U^0.005 -> 0.86746, U^0.01 -> 0.86825,
  U^0.02 -> 0.86979, U^0.04 -> 0.87277. Existing marked inverse moment covers only r+2z<=1 (prime slots); plain x plain at total length 1 is proved.
- Exact conditional formula (sextic): if the sextic family satisfies a DH-type row count U^{1+lambda-delta} for all delta, then
  theta = 11/12 - (1-2 lambda)/(4(5+6 lambda)) (lambda=0: 13/15). Matches grid model.
- DH window needed (sextic/dhneed.py, base OpenAI row counts, no Cycle61 refinements):
  0.8745: sigma in [0.647,0.746]; 0.8740: [0.633,0.788]; 0.8730: [0.599,0.820]; 0.8720: [0.581,0.848]; 0.8710: [0.579,0.871];
  0.8700: [0.558,0.870]; 0.8690: [0.539,0.869]; 0.8680: [0.522,0.868]; 0.8670: [0.510,0.867].
  => sextic improvements also need DH-strength zero density for a thin family over a broad middle range (unknown).
- CORRECTION (sextic/mixed3.py): detector witnesses satisfy r+m in [t, t+1/2], not r+m=t. A mixed-moment count 1-delta(r+m) grows with r+m,
  so the adversary escapes it with r+m>1. MM at critical length only (r+m<=1) gives NO gain: theta_max 0.875058 = baseline.
  The earlier "MM => 13/15" was an artifact of fixing r+m=t.
- What gives DH (=> 13/15): MM' = sum_u |M_u(U^r) S_u(U^m)|^2 << U^{1+eps} for ALL r<=1, m<=1/2 (supercritical 1 < r+m <= 3/2 included).
- Why MM' is beyond current methods: Poisson in u turns length s=r+m into a dual with U^{2s-1} rows and U^s columns. For s>1 (rows > columns)
  orthogonality only gives U^s, so the dual sums need more than square-root cancellation, i.e. a reflection. Theta reflection of the Mobius
  (cubic Gauss) part shortens the dual only if 4s-2-r < r, i.e. r + 2m < 1 = exactly the paper's marked constraint r+2z<=1 (intrinsic, not slack).
  Detector pairs always have r+m>=1, so the marked method can never reach them. The plain part's dual is a sextic Kubota coefficient at index h
  (no mu, no theta), reflecting only when r<1/2.
- Required density exponent (sextic/needcurve.py; A defined by #rows with zero Re>=sigma << U^{A(1-sigma)}, DH=2):
  theta 0.8745: A<=2.155..2.28 on sigma in [0.70,0.80]; 0.874: A<=2.14..2.25 on [0.68,0.82]; 0.872: A<=2.09..2.21 on [0.62,0.86];
  0.870: A<=~2.065 on [0.58,0.86].
  Calibration (analogy, different family/aspect): for zeta the best known A in sigma~0.70-0.78 is Guth–Maynard 30/13≈2.308; DH known only
  for sigma>=25/32 (Bourgain). So even theta=0.874 via the sextic architecture needs a density exponent beyond what is known for zeta
  in sigma∈[0.70,0.78], for a thin family whose Halász off-diagonal (chi_{uu'^5}, conductor ~U^2) is much weaker than zeta's.
