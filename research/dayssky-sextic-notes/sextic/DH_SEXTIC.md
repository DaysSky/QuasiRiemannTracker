# The density hypothesis for the sextic twist family: reductions, partial results, barriers

DaysSky with Claude assistance, 16 Oct 2026. Labels: PROVED, MODEL (exponent algebra of the 7/8 architecture), NUMERICAL, OPEN.

Family: rows u (sixth-power-free in Z[ω], N(u) ~ U), characters ψ_u = ν·χ_•(u)^{±1} (ν in the fixed finite group Θ).
DH(window): #{u : some ψ_u has a zero with Re ≥ (1+δ)/2, bounded height} ≪ U^{1−δ+ε} for δ ∈ [1/50, κ].

## 1. Why this is the input for 13/15 (PROVED relative to the 7/8 architecture)

Conditional Theorem D (REPORT §3): a count U^{1+λ−δ} in every bin gives θ = 11/12 − (1−2λ)/(4(5+6λ)); λ → 0 gives
zero-freeness of every Dirichlet L-function for Re s > 13/15. No prime slots are needed once DH holds.

## 2. Four equivalent-strength single inputs, each sufficient (MODEL, sextic/dhforms.py)

Detector at t = 1: every bad row has a tail block of length s ∈ [1, 3/2] with saturated witnesses, in two
representations: (i) M_r·S_m, r ≤ 1, m = s − r ≤ 1/2; (ii) M_{r'}·S_{m'}, r' ∈ [1, s]. Counts available today:
1 − δr (Möbius second moment, r ≤ 1), e(r) − δr with e(r) = (1+5r)/6 (amplified, r > 1), 1 − 2δm (plain pairs).
Adding any ONE of the following gives exactly R(δ) = 1 − δ (DH) for every δ (checked at δ = 0.1 … 0.75):

| input | statement | strength |
|---|---|---|
| MLV | #{u : \|M_u(U^r)\|² ≥ U^{rδ}} ≪ U^{r(1−δ)+ε}, r ≤ 1 | Möbius large values (Montgomery-type) |
| PLV | #{u : \|S_u(U^m)\|² ≥ U^{mδ}} ≪ U^{1−δ+ε}, m ≤ 1/2 | plain large values |
| MM′ | Σ_u \|M_u(U^r)S_u(U^m)\|² ≪ U^{1+ε}, r ≤ 1, m ≤ 1/2 | supercritical mixed moment |
| LM | Σ_u \|M_u(U^s)\|² ≪ U^{1+α′(s−1)+ε} with α′ ≤ δ, s ≤ 3/2 | long Möbius mean square (α′ → 0 for the window) |

Without them: R(0.467) = 0.617 (= Cycle61's ρ), versus DH 0.533. With LM at slope 0.45 only δ ≥ 0.45 is reached.
Remark: MLV holds at r = 1 already (the Möbius second moment gives U^{1−δ}); the whole difficulty is r ∈ (1/2, 1).

## 3. Why each route is blocked with known tools

- LM (PROVED, duality): Poisson in u turns Σ_u|M_u(U^s)|² into a dual of U^{2s−1} cubic-theta rows against U^s
  columns; e(s) = 1 for s > 1 needs the dual average below square-root size, which can only come from the mean-zero
  of the annular weight, and Poisson in h returns the original family (circular; REPORT §7b). Linear copy
  amplification (sixth powers) is optimal and gives the slope 5/6 (REPORT §7a).
- MM′, MLV via theta reflection (PROVED): reflection shortens the dual only if r + 2m < 1; detector pairs have r + m ≥ 1.
- MLV, PLV via Halász–Montgomery/Huxley (PROVED, exponent algebra, sextic/hmfam.py): Gram entries are pair sums of
  χ_{uu′⁵} with conductor ≈ U², the square of the family size (for ζ the conductor equals the family size).
  - With Burgess–Wu for the pair characters (exponent 3/16 of the norm), HM never applies in the window
    (σ = 0.85 needs polynomial length U^{1.88}).
  - Huxley with the pair-family fourth moment gives count exponents 1.01, 1.00, 0.667, 0.429 at
    σ = 0.70, 0.75, 0.80, 0.85, against the current counts 0.683, 0.583, 0.474, 0.352.
  - Even pointwise Lindelöf for the pair family gives DH only for σ ≥ 3/4, as for ζ.
- PLV via moments (MODEL): critical moments of every integer order still leave m ∈ (1/(k+1), 1/k) open (integer
  power obstruction); the supercritical sixth moment needs loss ≤ 0.2δ (≈ 0.09 at the binding δ); the best available
  loss is 2(β*−½)m ≈ 0.34 (pointwise bound inside the contradiction argument × plain fourth moment).
- Bootstrap from the contradiction hypothesis (PROVED, SESSION3 §1): zero-freeness beyond β* never improves a count
  at a zero of real part ≤ β*.
- Averaging the probe over targets (PROVED, SESSION3 §6.8): re-encodes only the multiplicativity of μ.
- Literature (Oct 2026 search): no density-hypothesis-type result for any thin family (quadratic, cubic, sextic);
  for quadratic characters the best is Heath-Brown's Ingham-type A = 3/(2−σ); our family's count already beats it
  (A = 2.275 vs 2.307 at σ = 0.70) thanks to OpenAI's Möbius second moment and prime compensation.

## 4. Numerical evidence (NUMERICAL, sextic/mmnum)

- Möbius fourth moment ratio M4/(2M2²) at r = 1 matches the random multiplicative control (1.09 vs 1.09 at U = 3·10⁵).
- Plain sixth moment at m = 1/2 (U = 10⁶, exceptional λ-power rows removed): 1.126 vs control 1.133.
- Large-value tails of |M_u(U^{0.75})|²: same heavy tail as random multiplicative (t = 16 mean: 64× vs 47× Gaussian),
  max/mean 19 among 3.6·10⁶ rows: no power-size large values. All four inputs MLV, PLV, MM′, LM are numerically
  consistent with the data.

## 5. What a proof needs

A large-values estimate for Möbius or plain polynomials over a family whose pair conductor is the square of its size,
in the range σ ∈ (1/2, 3/4) where even for ζ only Guth–Maynard-type ideas beat the classical bounds. A route would
have to supply cancellation in the pair-family Gram matrix restricted to an arbitrary bad set B (|B| = U^{R}),
i.e. bound ‖G_B − N·I‖ below |B|·N^{1/2}, using the arithmetic of sextic reciprocity rather than size.

## 6. The irreducible core: fractional critical moments (PROVED, exponent algebra)

At t = 1 the adversary's witness pair is (r, m) = (1−m, m). Möbius gives (1−m)δ; plain 2k-th moments at critical
length (km ≤ 1) give k·m·δ with k = ⌊1/m⌋. On m ∈ (1/3, 1/2) only k = 2 exists, and at m = 1/3⁺ both counts equal
(2/3)δ. Hence integer moments of every order, even all proved at critical length, leave R = 1 − (2/3)δ at t = 1.
DH therefore needs, at minimum, the fractional critical moments Σ_u |S_u(U^m)|^{2/m} ≪ U^{1+ε} for m ∈ (1/3, 1/2)
(equivalently PLV there), or a substitute (MLV on r ∈ (1/2, 2/3), MM′, LM). For ζ the analogous fractional large-values
statements are reached only by Halász–Montgomery for σ ≥ 3/4, and for this family §3 shows HM is unavailable.

Even the integer critical sixth moment (three plain factors of total length ≤ 1) is out of reach of OpenAI's plain
recursion: two factors with total > M′ always contain one factor > M′/2, which Hecke reflection shortens; three
factors can exceed M′ with all factors < M′/2, so the per-step width tax σ accumulates with nothing to pay it.
Available substitutes: BGL for S³ loses U^{1/3}; the pointwise bound inside the contradiction argument loses
U^{2m(β*−1/2)} = U^{0.25} at m = 1/3.

## 7. Price list at the binding class (MODEL, sextic/lmsens.py, sextic/bindpt.py)

Binding witness at δ = 0.466, p = ½ (record geometry): t* = 1.162, Möbius piece r* = 0.745, plain piece m* = 0.417,
three-way tie A_I = A_P = long = 0.5935 (DH value 0.534). θ in the Cycle61 grid model (baseline 0.87492) under single changes:

| change | θ |
|---|---|
| LM slope α′ = 0.80 / 0.75 / 0.70 / 0.60 / 0.50 / 0.40 / 0.20 / 0 | 0.87467 / 0.87422 / 0.87372 / 0.87222 / 0.87055 / 0.86940 / 0.86782 / 0.86667 |
| Möbius count −e·δ, e = 0.01 / 0.02 / 0.05 / 0.1 / 0.2 | 0.87468 / 0.87424 / 0.87333 / 0.87172 / 0.86841 |
| plain count −e·δ, e = 0.01 / 0.02 / 0.05 / 0.1 / 0.2 | 0.87481 / 0.87472 / 0.87439 / 0.87388 / 0.87284 |

LM alone at α′ = 0 lands exactly on 13/15 (consistency check of Conditional Theorem D). The Möbius count is the most
θ-efficient single lever; the plain count is about three times less efficient.

## 8. Benchmark against zeta (PROVED exponent comparison, sextic/zetabench.py)

Same σ = (1+δ)/2, counts as exponents of the family size (zeta: A(1−σ) = A(1−δ)/2):

| δ | σ | ρ (p=0) | R (p=½) | Ingham | Huxley | Guth–Maynard | DH |
|---|---|---|---|---|---|---|---|
| 0.40 | 0.700 | 0.6825 | 0.6561 | 0.6923 | 0.7200 | 0.6923 | 0.600 |
| 0.466 | 0.733 | 0.6182 | 0.5935 | 0.6322 | 0.6408 | 0.6162 | 0.534 |
| 0.50 | 0.750 | 0.5833 | 0.5602 | 0.6000 | 0.6000 | 0.5769 | 0.500 |
| 0.60 | 0.800 | 0.4737 | 0.4575 | 0.5000 | 0.4800 | 0.4615 | 0.400 |

The family's current count with prime slots already beats the zeta analogue of Guth–Maynard at every σ of the binding
range. DH for the family at δ ∈ [0.40, 0.56] would therefore be stronger, relative to family size, than every
zero-density exponent known for ζ at the same σ (Bourgain's DH range for ζ starts at σ = 25/32, δ = 0.5625).

## 9. Guth–Maynard and its q-aspect version do not transfer (PROVED, structural)

Cubic trace for M_{u,n} = ψ_u(n): tr((M_B M_B^*)^3) = Σ_{u1,u2,u3 ∈ B} K12 K23 K31, K(u,u′) = Σ_n w(n/N) ψ_u ψ̄_{u′}(n).
Poisson in n mod uu′ and multiplicativity of Gauss sums for coprime moduli give
  τ(ψ_{u1} ψ̄_{u2}) = ψ_{u1}(u2) ψ̄_{u2}(u1) τ(ψ_{u1}) τ(ψ̄_{u2}),
so the cyclic product of the three Gauss sums is N(u1u2u3)·(a residue-class unit from sextic reciprocity and ψ(−1)).
Hence, after splitting B into boundedly many classes,
  S3 = (N/U)^3 · c · Σ_{h1,h2,h3 ≲ U²/N} F_B(h3/h1) F_B(h1/h2) F_B(h2/h3),  F_B(x) = Σ_{u∈B} ψ_u(x),
which is the same cubic trace at the dual length H = U²/N. Equivalently, the Gram matrix obeys the exact reflection
  M_B M_B^* = N·I + (N/U)·D (M′_B M′_B^* − H·I) D^*,  D = diag(normalized Gauss sums),
so the spectrum at length N is an affine image of the spectrum at length H.
- Guth–Maynard (ζ) and Chen–Gupta–Li (arXiv 2507.08296, all characters mod q) gain by NOT evaluating the archimedean
  integral (resp. the common-modulus Gauss sum). The leftover variable (u ∈ [1,2]³, resp. a ∈ (Z/q)³) yields a linear
  constraint on the dual frequencies and continuous L²/L⁴ orthogonality for the R-function. In the sextic family the
  moduli u1u2, u2u3, u3u1 differ. By CRT the residue variables split into the separate Gauss sums, which cancel
  cyclically, so no free variable is left. The dual R-function problem is literally the original problem.
- The discrete L⁴ of F_B on the dual box is |B|²H² + |B|⁴H (no orthogonality beyond conductor H), against the continuous
  E(W) in Guth–Maynard. Multiplicative energy is minimal (≤ |B|²U^ε), but this cannot be used without the integral.

## 10. Heath-Brown identity decomposition gives nothing new (PROVED, exponent algebra)

Decompose M_r into tiny Möbius atoms and smooth atoms (reflected to length ≤ 1/2). Available group moments: k = 2 for a
single smooth atom ≤ 1/2 or a Möbius group ≤ 1/3 (marked family, r1 + 2r2 ≤ 1); k = 1 otherwise. With shares x_G the
count is 1 − δx / Σ_G (1/k_G), and beating 1 − δr needs Σ_G 1/k_G < 1. The adversary picks smooth ≈ 1/2 and Möbius ≈ r − 1/2,
giving Σ 1/k = 1. So any gain needs a k = 3 group (a critical sixth moment), which is again the integer-power obstruction (§6).

## 11. Large-sieve side: Dunn–Radziwiłł and the (KD)^{2/3} term (CITED + PROVED reading)

- OpenAI's sextic large sieve (Lemma sextic-large-sieve, BGL recursion) is K + D + (KD)^{2/3}. The bracket
  [max(v, z_a, 2(v+z_a)/3) − v] in the reflected energy is exactly this (KD)^{2/3} term, and it is where the slot cost 2
  of the inverse lemma comes from (§9.1 of SESSION3).
- Dunn–Radziwiłł (arXiv 2109.07463, Thm 1.4) show that under GRH the cubic term (AB)^{2/3} is attained for A ∈ [√B, B²]
  by sequences β_b g̃(b) (Gauss-sum twisted, β ≥ 0 uncorrelated with cubic symbols). The column sequences produced by
  reflection are Gauss-sum twisted, so this is the extremal case, not an artefact.
- Their corrected sieve (Prop. 9.1/9.2) saves only A^{−ρ}B^{−ρ} with ρ(ε) < 10⁻⁴. It needs β ∈ C_η(B, w), i.e. square-root
  cancellation of β against every cubic Hecke twist (verified only under GRH).
- Consequence: a sharp K + D sieve, which would lower c_I from 2 toward 1, is false for the relevant sequences.
  A pure prime-slot moment needs a sharp sieve range a0 ≥ 0.87 at the binding class (count 1 − 2qa0 < 0.5935). The
  BGL/Heath-Brown range is a0 = 1/2, giving 1 − q = 0.767 (no effect: θ stays at 0.87492 in the grid model, sextic/lsconj.py).

## 12. Using every t at once does not help (PROVED, exponent algebra; paper Prop. detector-witness)

The detector holds for EVERY t ∈ [1, 3/2] with the same zero. Per row, write 𝒟 = {r : |M_r|² ≥ U^{δr−ε}} and
𝒩 = {m ≤ 1/2 : |S_m|² ≥ U^{δm−ε}} (0 ∈ 𝒩 always). The proposition for all t is equivalent to the covering condition
  [1, 3/2] ⊆ ∪_{r∈𝒟} [r, r + m_max],  m_max = max 𝒩.
Costs: c_M(r) = A_I(r) for r ≤ 1 (decreasing, value 1 − δ at r = 1), c_M(r) = (1+5r)/6 − δr for r > 1 (increasing),
c_S(m) = A_P(m). A row costs min(c_S(m_max), min_{r∈𝒟} c_M(r)). Every covering chain either contains r = 1 (DH cost)
or has consecutive points p_a ≤ 1 < p_b with p_b − p_a ≤ m_max. The adversary's value is therefore
  max_{r+m=s} min(A_I(r), A_P(m), long(s)) maximised over s, = max_s min(short(s), long(s)).
Since short decreases and long increases in s, this equals min_t max(short(t), long(t)) = R. The single-t count is
lossless relative to the all-t covering structure. A gain would need a JOINT count of rows saturated at two Möbius
lengths p_a < 1 < p_b. In the explicit-formula picture this is a dip of Σ_j c_j D^{ρ_j−ρ′} at D = U, which two zeros at
height distance ≍ 1/log U already produce, and no tool counts such configurations.

## 13. What a sharp sextic large sieve would buy (MODEL, sextic/lsconj.py, Cycle61 grid, baseline 0.87492)

Hypothesis LS(a0): Σ_{u~U}|Σ_{n≤U^a} c_n ψ_u(n)|² ≪ U^{1+ε}‖c‖² for squarefree c and a ≤ a0, which turns pure prime-slot
moments into the count 1 − 2qa0. Separately, slot costs c_I, c_P are lowered toward the sharp-sieve values.

| a0 | c_P | c_I | θ |
|---|---|---|---|
| 0, 0.5 (BGL), 0.75 | 6 | 2 | 0.87492 |
| 0.9 | 6 | 2 | 0.87450 |
| 1 | 6 | 2 | 0.87398 |
| 1 | 6 | 1 | 0.87237 |
| 1 | 4/3 | 1 | 0.86969 |
| 1 | 1 | 1 | 0.86960 |

Even the full sharp sieve stops at ≈ 0.8696 > 13/15, because neither A_I nor long(t) reaches the DH line at the binding
class. 13/15 needs a DH-type large-values input. By §11 the sharp sieve is false for Gauss-twisted sequences anyway.

## 14. Moment literature for r-th order families (CITED, Oct 2026 search)

- David–de Faveri–Dunn–Stucky (arXiv 2410.03048): mollified second moment, with power saving, for cubic Hecke
  L-functions over Q(ω). Short mollifier; relies on Heath-Brown's cubic sieve.
- Diaconu–Ion–Paşol–Popa (arXiv 2607.27131, Jul 2026): twisted first and second moments for r-th order families via
  multiple Dirichlet series, r ≥ 3, over fields containing μ_{2r}. The r = 6 case needs μ_12, so it does not cover the
  sextic family over Q(ω). The second-moment error is X^{(r+A_r)/(r+1)+ε}, and the twists are polynomially short.
- de Faveri–Dunn–Hoffstein (cited there as [DFDH26]): the cubic (AB)^{2/3} cannot be removed, with a refined prediction
  for r ≥ 4.
- Relevance: MM′ (the mixed Möbius×plain moment, a single sufficient input for DH, §2) is a twisted second moment of
  L(s, ψ_u) with Möbius twists of length U^r, r ≈ 0.745 at the binding class. The known twist ranges are polynomially
  small (for ζ the record twisted length is Conrey's 4/7). So MM′ lies far outside them.

## 15. Near the critical line (PROVED exponent comparison)

Slopes in δ of the count at δ → 0 (DH slope = 1): this family ρ: 2/3; with saturated slots R(δ, ½): ≈ 0.82.
ζ-analogue (Selberg/Conrey mollifier, θ = 4/7): 4/7. The family's current near-line slope already beats the best
ζ-mollifier slope, and DH near the line would need mollifier length equal to the conductor.

## 16. The n-th order sieve is not orthogonal: de Faveri–Dunn–Hoffstein (CITED, arXiv 2607.07911, Jul 2026)

- Unconditional: Ξ_3(A,B) ≫ (AB)^{−ε}(A + B + (AB)^{2/3}) (cubic, BGL is optimal); Ξ_4 ≫ A + B + A^{3/4}B^{1/2} + A^{1/2}B^{3/4}.
- Conjecture (1.8): Ξ_n(A,B) = (AB)^{o(1)}(A + B + A^{1−1/n}B^{2/n} + A^{2/n}B^{1−1/n}). They call any power saving on
  (AB)^{2/3} for n ≥ 4 a substantial advance.
- Sextic consequence (rows A = U, columns B = U^a): worst-case sieve U^{max(1, a, 5/6 + a/3)}, so it is sharp only for
  a ≤ 1/2, even conjecturally. For untwisted non-negative sequences (pure prime-slot products), duality puts the extremal
  twisted case at (B²/A, B) and gives A^{1/n}·B = U^{1/6 + a}: sharp only for a ≤ 5/6.
- With §13 (a pure-slot count needs a0 ≥ 0.87 at the binding class), the pure-slot lever is closed even under the optimal
  conjectured sieve. Möbius-signed sequences (marked M·Q) have no bias main term, so the cost-1 direct sieve for M·Q is
  not excluded by this obstruction. Its value is modelled in sextic/lsconj2.py.

## 17. Session 14:37Z: restricted Gram numerics (NUMERICAL, mmnum/gramB.py)
U = 10⁵, D = U^{0.75} (793 squarefree columns, Σw² = 388). Large-value row sets B (top |M_u|) vs random sets:
s1(M_B)² = 1307 / 2635 / 6337 (top 100/300/1000) vs 750 / 1145 / 2242 (random); S2_B/(|B|²Σw²) = 1.00 / 0.92 / 0.85
(top) vs 0.70 (random). The top sets sit at |M_u|² ≈ 10× the mean, far below the DH-relevant threshold U^{δr} ≈ 10^{1.7}.
At reachable U no row comes near the bins that matter, which is the expected (GRH) behaviour, so the restricted-Gram
question at the binding sizes (|B| = U^{0.59}, N = U^{0.745}) cannot be probed numerically. Needed there: s1(M_B)² ≤ U^{0.94}.
Random-matrix size is U^{0.75}, the full-family sieve gives U^{1+}, and pointwise Lindelöf for pairs gives |B|N^{1/2} = U^{0.97}.

## 18. Jensen/Cauchy dichotomy at the zero (PROVED reduction + MODEL, sextic/jensen.py)
Fix the bin zero ρ. Case A: |L| ≤ U^{2η} on [1/2, σ] × (height window). Shifting the plain Mellin contour to Re = 1/2
gives |S_m|² ≤ U^{4η}, so saturation forces m ≤ 4η/δ. The detector at t = 1 then has r ≥ 1 − 4η/δ, giving count
1 − δ + (δ−q)(4η/δ). Case B: |L(s1)| ≥ U^{2η} at some Re s1 ∈ [1/2, σ]. The AFE gives a plain piece large at s1, and the
plain 4th moment gives count 1 − 8η (worst point Re s1 = 1/2). Optimising: δ = .466 → 0.6278; δ = .3 → 0.7606;
δ = .6 → 0.5207. These are close to, but above, OpenAI's ρ(δ) (0.618 at .466), and they never fall below the record's
R(δ,½) = 0.5935: case A needs η < 0.030 while case B needs η > 0.051. It is an independent route to ρ-type counts from
the same two inputs (Möbius second moment, plain fourth moment). It would beat R only with central large-value counts
stronger than the fourth moment (e.g. a 6th moment with loss < 0.17).

## 19. Explicit-formula / Soundararajan / Turán routes (exponent analysis)
- Soundararajan-type upper bound for log|L(1/2+it)| at a bad row: the bad zero at β = σ must lie to the LEFT of the
  auxiliary abscissa σ0. Otherwise its x^{β−σ0} term destroys the bound, since prime-sum length x ≥ (log U)^C makes it
  a power of log U. With σ0 ≥ σ the convexity part (σ0 − 1/2)/2·log U equals δ/4·log U. That is exactly the per-row
  bound |L(1/2)|² ≤ U^{δ/2} already used by Cycle61/OpenAI (§9.2 of SESSION3). It gives nothing new.
- Rep (ii) via the explicit formula: a long saturated Möbius witness with no saturation near r = 1 needs interference
  among zero terms c_j D^{ρ_j−ρ′}, c_j = W̃(ρ_j−ρ′)/L′(ρ_j). Nazarov–Turán needs ≳ η log U active terms for a U^{−η} dip.
  Jensen allows O(log U) zeros in a unit box, and the residues 1/L′(ρ_j) are uncontrolled, so no power saving follows.
- Jensen/Cauchy dichotomy (§18) and Heath-Brown-identity decompositions of witnesses (M_r ⇒ M_{r−a}P_a, marked
  M·P needs r + a ≤ 1 − a) only add witnesses with worse counts. The closure of the available toolkit stays at R(δ,p).
- Slot decoupling: the marked moment equals the mean square of b(m) = −μ(m)ω_P(m) at length U^{r+z}. A sharp (cost 1)
  bound is not excluded by Gauss-sum bias, because b is Möbius-signed. In theta reflection the slot primes sit in the
  Poisson modulus, giving a dual length ∝ P². A cost-1 proof would need square-root cancellation in the slot variable of
  Σ_p g̃(p)χ_p(h)B_{h,p}, a Type II sum for sextic Gauss sums with one very short variable. BGL's sieve in h is sharp
  for the short variable alone, but the coupling through B_{h,p} blocks Hölder/Cauchy (needs the dual 4th moment of
  Möbius–Gauss sums at length D² ≫ H_d). This is the most concrete open sub-problem found. The θ value of
  c_I = 1.9 / 1.8 / 1 is 0.874714 / 0.874445 / 0.87289 in the grid model (baseline 0.87492).

## 20. Literature check, Oct 2026 (CITED)
- Baiying Liu, arXiv 2610.12234 (12 Oct 2026, Lean-formalised): 0.874957019…, varying (b, ℓ, κ) in OpenAI's 7/8 proof.
  Its new analytic input (Prop. 28.1) extends the plain moment to κ ∈ [13/18, 1] with slot condition n1 + n2 + 6κz ≤ M.
  This is the κ-dependent plain slot cost already used by Cycle61 and the record (z_P = (1−2m)/(6κ), κ = 2θ−1).
  Count (29.2): R_κ = 1 − δ + (α−δ)δP/(2J) with α = 5/6, the same short/long crossing. Optimality is proved inside
  that comparison system. Nothing new for the record (0.8748389 < 0.874957).
- Structure confirmed: short(3/2) = 1 − δ and long(1) = 1 − δ. The excess over DH comes only from the long slope α = 5/6
  meeting the short slope. LM slope α′ < 5/6 is the "variance of μ against the family beyond its conductor" problem. For
  the full Dirichlet family mod a single q this is the individual-modulus variance Σ_{a mod q}|Σ_{n≤x, n≡a}μ(n)|² ≪ x^{1+ε}
  for q^{1+ε} ≤ x ≤ q^{3/2}, which is also not known by current methods.

## 21. Toolkit closure (PROVED by explicit adversary, sextic/closure.py)
Proposition. Let the available inputs be: the marked inverse moment (A_I, r ≤ 1, slot cost 2); the plain lemma (A_P,
m ≤ 1/2, slot cost 6κ / Thm 2.2); sixth-power copies for ANY polynomial of length s > 1 (e(s) = (1+5s)/6, including
products of witnesses); Möbius 4th moment via copies; marked Möbius×Möbius (r1 + 2r2 ≤ 1); theta-reflection MM′
(r + 2m < 1); pure slot moments in the BGL range (a0 = 1/2); the pointwise plain sixth moment (loss 0.375); and all
Hölder interpolations of these. Feed in the detector for every t ∈ [1, 3/2] (covering form, §12). Then the bad-row count
exponent equals R(δ, p) exactly. Witness: the adversary row with Möbius saturated at {r*, t*} and plain at m* = t* − r*
(chain p_a = r* ≤ 1 < p_b = t*), checked at δ = .3, .4, .466, .55, .65.
- At δ = 0.466: A_I = A_P = long = 0.5935. M4-copies 0.714, (r*, t*) product copies 0.867, pure slots 0.767,
  plain-6th pointwise 0.792, MM′ and M×M-marked inapplicable.
- Hölder cannot go below the minimum of its vertices.
So the record count is the closure of the toolkit. Any gain requires a genuinely new mean-value input (list in §7/§13).

## 22. The long slope 5/6 is an exact fixed point of Poisson + copies (PROVED, exponent algebra)
MS(U;D) := Σ_{u~U} w(u/U)|M_u(D)|². Poisson in u (Gauss sums factor) gives the exact involution
MS(U;D) = UDŵ(0) + (U/D)·Σ_{h≠0} ŵ(hU/D²)(|B_h(D)|² − D),  B_h(D) = Σ_n μ(n)g̃(n)χ̄_n(h)W(n/D).
The dual of the dual is the original with coefficient exactly 1. So for U < D the statement "MS(U;D) ≪ UD" is
equivalent to an asymptotic dual large sieve with power-saving error O(D²) at K = D²/U > D, and that error's dual is
again MS(U;D). Copies give the only contraction (factor P^{−1} = (H/U)^{−1/6}, needing H ≥ D), hence the slope
(1+5s)/6. Breaking it needs a third ingredient:
(i) sub-square-root Möbius–Gauss sums on average, impossible since the dual diagonal forces Σ_h|B_h|² ≍ KD;
(ii) spectral/Kloosterman input for Möbius-twisted coefficients of the sextic metaplectic Eisenstein series;
(iii) GRH-type input.

## 23. Cycle61 Theorem 2.2 window extension (PROVED reserves + MODEL θ): no gain
- The validity ranges z ∈ [0.034, 0.04], n ∈ [0.41, 0.44] are not used in the body of the proof beyond the cell bounds.
  Recomputing the affine reserves exactly gives:
  - error fixed-inducing 1 − A_j − 9ā_j = 0.174703 − 3.752246·b > 0 iff b < 0.046560;
  - second zero (k = 5) holds for b < 0.1386; main fixed-inducing (k = 4) is increasing in b;
  - comparison budget (2.12) holds for z < 0.0556; support gate 1/4 − 7ā_max > 0 holds;
  - low and high reserves do not depend on b; g_j > 0 iff b > 0.024093.
  So Theorem 2.2 extends (pending the "original support gates", including A_max ≤ 0.9157) to z ∈ (0.02409, 0.04656),
  i.e. a detector plain window m ∈ [0.4005, 0.4458] instead of [0.4137, 0.4235].
- θ effect (record optimizer c61/thoptw.py with configurable window): 0.8748446 for the current window, z-extended, and
  z+n-extended alike (identical to 7 digits). The binding plain witness stays at m* ≈ 0.417 inside the old window, and only
  the capacity there, (B − 2m*)/Λ, matters. Window width is not a lever. The child-density bootstrap is also dead:
  R(δ, 0) = ρ(δ) because slot gains vanish at p = 0.

## 24. Long witness via Vaughan/Heath-Brown identities + in-bin pointwise caps (PROVED, exponent algebra)
Setting: row u in bin (i, a), δ = 2a − 1. Lemma buffered-bins gives |L|, |1/L| ≪ U^ε on Re s ≥ a + 6e (height ≤ (3i+2)T1)
and the reflected bound on Re s = 1 − a − 6e. Consequences for normalized excesses (|X|² = U^{len + exc}):
- (C1) every Möbius piece of length b (any b): exc ≤ δb (Perron, shift to a + 6e);
- (C2) every plain piece of length c ≤ 1: exc ≤ δ·min(c, 1 − c) (direct for c ≤ 1/2; Hecke reflection for c > 1/2);
  plain pieces of length c > 1 + ε are O(U^{−A}).
So a piece of a decomposition can carry the full excess tδ of a saturated M_t only if every factor is saturated and
every plain factor has c ≤ 1/2.
- Vaughan (parameters y', z' with y' + z' < t − 1): M_t = (Type I, negligible since its plain factor has length > 1)
  − Σ μ_{>y}(b)·λ_z(c). Saturation forces the Type II piece M(β)·F(d)·S(m) with m ≤ 1/2, d ≤ z', β ∈ [t − 1/2 − z', t − z'].
  Iterating gives a chain of saturated Möbius witnesses descending by ≤ 1/2 per step, each step bridged by a saturated
  plain witness, down to length ≤ 1 + y' + z'. This is exactly the detector's covering chain (§12); the adversary's
  configuration {r*, t*}, m* = t* − r* is already closed under it. No new constraint.
- Heath-Brown identity (k-fold): pieces with two plain factors (c1 + c2 ≈ t > 1) are never saturated by (C2);
  saturated pieces need ≥ 3 plain factors ≤ 1/2 or Möbius-short products. Three-plain pieces give count
  min_i A_P(c_i) ≤ 1 − 2tδ/3 − (slot gain), which beats copies only for t > 5/(5 − 2δ) (= 1.229 at δ = .466).
  All-Möbius-short pieces (μ_{≤y})^{*k} reproduce the original problem and fall back to copies. Net: long(t) unchanged
  on (1, 1.229], the adversary keeps p_b = t* = 1.162, and R is unchanged.
- Structural by-product: modulo the Möbius-short products, the long witness is controlled by plain moments of order
  ≥ 6 at supercritical total length (three plain factors ≤ 1/2, total t ∈ (1, 3/2]). The LM lever and the plain-sixth
  lever are the same obstruction seen twice.

## 25. Slot cost c_I = 2 as a twisted-moment equidistribution problem (PROVED reformulation)
|M_u(N)P_u(Z)|² = Σ_{p1,p2} ψ_u(p1)ψ̄_u(p2)|M_u|². The diagonal is Z·Σ_u|M_u|² ≪ UNZ. The cost-1 statement (NZ ≤ U)
is equivalent to the off-diagonal bound
  Σ_{p1 ≠ p2 ~ Z} Σ_u χ_{p1}χ̄_{p2}(u)|M_u(N)|² ≪ U^{1+ε}NZ,
i.e. equidistribution of the Möbius mean-square weights |M_u(N)|² over residue classes u mod q = p1p2 (q ≍ Z², tiny)
with relative error q^{−1/2} per class. Poisson in u with the extra modulus q raises the dual length from N²/U to
qN²/U and multiplies the prefactor by √q. Requiring dual length ≤ N gives NZ² ≤ U, which is exactly c_I = 2.
Applying Poisson to the dual (mod q·n1n2) returns the twisted original with coefficient 1 (same involution as §22,
now with a tiny-modulus twist). Per-class upper bounds (Möbius second moment on the sublattice a + qO) give no saving.
A saving needs an asymptotic Σ_{u≡a(q)}|M_u|² = main + O((U/q)N q^{−1/2}). The dual off-diagonal is of the same
order UN as the diagonal (large-sieve N² term at K = N²/U < N). Small h carry a possible Patterson-type bias
Σ_{n≤N} μ(n)g̃(n)χ_n(h) ≈ c(h)N^{2/3}, so no asymptotic is available without control of Möbius-twisted
sextic theta coefficients.

## 26. Where r < 1 is used in OpenAI's Möbius second moment (PROVED reading of the initialization, paper 11601–12340)
Lemma marked → canonical family by one masked Poisson in u. With columns c_i = C z_i (C = gcd), d′ | C from the mask,
P1 = B − θ, residual column length D′ = r + z − 2G, the child has
  F_init = D′ − P1,  M_init = 2D′ − m − 2P1,  puncture Q_init ≤ P1 + G,
and the two canonical invariants become
  F − M − Q − z0 = m − r − 2z + 2G  (≥ c1),   4F − 3M − 6z0 = 3m − 2r − 8z + 10G + 2P1  (≥ c2).
P1 cancels in the first, which is the binding one. Without slots it reads m − r > 0: after Poisson, the dual row range
M_init = 2r − m must stay below the effective column length F_init = r. At r = m the child sits exactly at the self-dual
point M = F. For r > m it lies in the rows-dominant regime, where only the diagonal bound Z^{M} is available, and the
parent inherits the (U + D)D large-sieve shape (slope 1). Copies improve this to 5/6.
Consequence: extending the raw moment to H ≤ D^{1+c0}, which would give slope 5/(6(1+c0)), so c0 = 0.042 already
beats the record, is equivalent to an off-diagonal (asymptotic) estimate for the canonical family at M > F0, with relative
saving Z^{−(M−F0)}. This is the same third ingredient as §22. The puncture invariant is a symptom, not the cause:
P1 cancels from it identically.

## 27. Convolution identities give only descending chains (PROVED, exponent algebra)
Exact identities on every row: for s > 0, Σ_j Σ_{a,b} μ(a)φ_j(a)ψ(ab)w(ab/U^s) = 0 (μ∗1 = δ), and
μ·log = −μ∗Λ. After Mellin separation each dyadic piece is a twisted product M(r)·S(s − r), resp. M(r1)·P_Λ(r − r1).
Constraint: at each length s, the number of large pieces is never exactly one, since a single large piece cannot be
cancelled by O(log U) pieces each smaller by U^{−η}. With the in-bin caps (§24) this gives:
- every saturated M(r) has a partner M(r′)S(r − r′) with r′ < r, r − r′ ∈ 𝒩 (≤ 1/2): a descending chain to 0;
- the detector's two representations (r ≤ t and r′ > t at the same block length s) are the case "split point = t";
- μ·log = −μ∗Λ forces a saturated M(r1)·P_Λ(r − r1). Prime sums have no reflection cap, so P_Λ may be saturated at any
  length, but long prime pieces only cost more (marked lemma: r1 + 2z ≤ 1).
The adversary extends consistently: 𝒟 = {0, r* − m* = 0.328, r* = 0.745, t* = 1.162, r″ ∈ (1.5, 1.579]},
𝒩 = [0, m*], with large pieces paired at s = 0.328, 0.745, 1.162, r″. Every added witness costs more than R
(A_I(0.328) = 0.69, long(1.55) = 0.74). So the closure value R(δ, p) survives every linear convolution identity among
μ, 1, Λ. A pair saturated as witnesses does not force its piece to be large (the Mellin τ-integral may cancel), so no
"every saturated pair needs a partner" constraint is available.

## 28. Numerics: the long Möbius mean square is flat to s = 3/2 (NUMERICAL, mmnum/lm.py)
MS(s) = mean_u |D^{−1/2} Σ_n μ(n)χ_n(u)W(N(n)/D)|², D = U^s, all sixth-power-free u with N(u) ≤ U, ratio to the diagonal:
| U | s = 1.0 | 1.1 | 1.2 | 1.3 | 1.4 | 1.5 | cols_eff/rows at s_max |
|---|---|---|---|---|---|---|---|
| 3·10³ | 0.946 | 0.938 | 0.926 | 0.931 | 0.953 | 0.928 | 1.05 |
| 10⁴ | 0.934 | 0.953 | 0.944 | 0.933 | – | – | 0.30 |
The random multiplicative ±1 control tracks μ exactly (0.926–0.946). The plain squarefree control grows (×2.3 at s = 1.5),
driven by a few rows (max |·|² ≈ 10³), as expected for plain sums beyond the conductor. Consistent with α′ = 0 (LM):
no off-diagonal mass is visible. See §44 for how far beyond the conductor these sizes actually reach (modest).

## 29. LM in the literature, quadratic analogue (CITED) and the size of the needed dual saving (PROVED arithmetic)
- Gao–Zhao, arXiv 2212.09241: S(X,Y) = Σ*_{d≤X} M(Y, χ_{8d})² = (4/π²)XY·h̃1·Z2(1) + O(X^{1/2}Y^{3/2+ε} + XY^{1/2+ε})
  UNDER GRH, valid only for Y ≪ X^{1−ε}. Beyond the conductor (Y > X), which is the LM range, there is no asymptotic
  even under GRH, and unconditionally there is no upper bound below the large-sieve/copies shape. For quadratic
  characters, copies by squares give slope 1/2 (Heath-Brown's sharp quadratic sieve as base). So LM is open for every
  family of fixed-order characters, quadratic included.
- Copies, read on the dual side: MS = UD + (U/D)E. Copies give E ≪ KD·(K/D)^{−1/6} (relative saving U^{−(s−1)/6}).
  LM (α′ = 0) needs E ≪ D² (relative saving D/K = U^{−(s−1)}). Any α′ < 5/6 needs relative saving (K/D)^{−1/6−η}.
  BGL alone gives Σ_h|B_h|² ≤ (KD)^{2/3}D, which is above KD (since K < D²). In the dual even "E ≪ KD" is not
  available from the large sieve; it comes only from copies.
- [Corrected in §56: via Hasse–Davenport the coefficients are twisted cubic Gauss sums.] A spectral (metaplectic Kuznetsov) approach would run through Type I sums Σ_m g(h, dm)W(dm/D), i.e. Fourier
  coefficients of the sextic theta function. Unlike n = 2, 3, these are not explicitly known for n = 6 over Q(ω), so even
  the main terms of the Type I contribution are not available.

## 30. Two further negative checks (PROVED structural / NUMERICAL)
(a) Möbius-heavy pieces are unavoidable in the standard combinatorial decompositions of μ (explicit for Vaughan,
Heath-Brown and Linnik; heuristic in general). Every such identity expressing
μ·1_{[D,2D]} through convolutions of 1, μ_{≤y}, μ_{>y}, Λ (Vaughan, Heath-Brown k-fold, Linnik) contains a piece
with no smooth factor of length ≥ U^{c}, e.g. (μ_{≤y})^{∗k} with y^k ≍ D. Otherwise 1/ζ would be a finite
combination of series each divisible by a long ζ-piece. For such pieces, LM is the 2k-th moment of a short Möbius
polynomial at total length y^k = D > U, which is again beyond the conductor (for ζ it is the analogous open
high-moment problem for partial sums of 1/ζ). So no identity-based proof of α′ < 5/6 exists without a new
beyond-conductor input. This complements §24 (large-values version) at the mean-square level.
(b) Prime-side explicit formula (task #44). At the binding p = ½ the prime slots are at their in-bin cap (2qz = δz,
which is why p ≤ ½). The prime explicit formula has residues 1, not 1/L′(ρ), so prime saturation at slot scales
(z ≤ η) only says that the zero cluster within height 1/(z log U) of γ0 adds coherently. Möbius saturation at r*, t*
with a dip at r = 1 is the damped-exponential sum Σ_ρ (W̃(ρ)/L′(ρ)) e^{(ρ−c0)ξ} vanishing at ξ = log U. Two zeros at
height spacing ≍ 1/log U with residues of opposite sign (forced for close pairs, c_{ρ1} ≈ −c_{ρ2}) already do this.
Prime saturation at small scales does not see the spacing ≍ 1/log U, so it gives no constraint. Plancherel against the
in-bin cap |1/L| ≤ U^ε on Re s = a + 6e is also slack: the spectral mass of a saturated M(U^r) there is
U^{−12er} ≪ the L² budget.
(c) Sextic Gauss sums, bias test (NUMERICAL, mmnum/gbias.py; twisted multiplicativity
g(np) = χ_n(p)χ_p(n)g(n)g(p) verified against brute force for all squarefree composites of norm ≤ 217).
Σ_{N(n)≤X} μ(n)g̃(n) = −43.9, +20.9, +72.5 at X = 1875, 7500, 30000 (|S|/X^{1/2} = 1.01, 0.24, 0.42): the sign
changes and the size stays square-root. Σ_p g̃(p)log p: |S|/(X^{1/2}log X) = 0.50, 0.03, 0.15. No Patterson-type
X^{2/3} bias is visible for sextic Möbius–Gauss sums. All partial sums are real (conjugation symmetry), and
g̃(p) = 1 for rational p = ππ̄. So the small-h dual contributions are square-root sized, and c_I = 2 is not explained by a
bias. It is the large-sieve (KD)^{2/3} term of the method (§25, SESSION3 §6.3).

## 31. Slot costs are artefacts numerically: c_I and c_P look natural (NUMERICAL, mmnum/slotnum.py, slotplain.py)
Ratio ρ = mean_u|X_u P_u|² / (mean_u|X_u|² · mean_u|P_u|²) over all sixth-power-free rows N(u) ≤ U, with P_u a prime slot
polynomial of length U^z and coefficients a(p) ∈ {1, conj(g̃(p)) (Gauss-correlated), random phases}.
Under natural (cost-1) capacity ρ stays O(1); a genuine cost-2 (resp. 6κ) obstruction would make ρ grow like a power
of U beyond the paper's ranges.
| X | U | (r or 2m, z) | paper range? | ρ (one / gconj / rand) |
|---|---|---|---|---|
| M_u(U^0.5) | 10⁵ | (0.5, 0.40): r+2z = 1.30 | no (needs ≤ 1) | 1.0013 / 0.9980 / 0.9997 |
| M_u(U^0.5) | 10⁶ | (0.5, 0.40) | no | 0.9993 / 1.0004 / 0.9996 |
| M_u(U^0.5) | 10⁶ | (0.5, 0.45): r+2z = 1.40 | no | 1.0000 / 0.9997 / 1.0006 |
| S_u(U^0.35)² | 10⁵ | (0.70, 0.30): 2m+4.5z = 2.05 | no | 1.0003 / 0.9979 / 0.9981 |
| S_u(U^0.35)² | 10⁶ | (0.70, 0.25): 2m+4.5z = 1.83 | no | 1.0001 / 0.9998 / 0.9996 |
| S_u(U^0.35)² | 10⁶ | (0.70, 0.30): 2m+z = 1.00 | no | 0.9994 / 1.0019 / 1.0000 |
| M_u(U^0.5) | 3·10⁶ | (0.5, 0.40) / (0.5, 0.45) | no | 1.0005 / 0.9993 / 0.9997 ; 1.0000 / 1.0000 / 1.0001 |
| S_u(U^0.35)² | 3·10⁶ (1.09·10⁷ rows) | (0.70, 0.25) / (0.70, 0.30) | no | 0.9996 / 0.9988 / 1.0004 ; 0.9992 / 0.9994 / 0.9995 |
(z = 0.2 rows have 1–2 primes in range and are degenerate; omitted.) From U = 10⁵ to 3·10⁶ there is no trend: a genuine
power loss U^{c} would show up as a drift ≥ U^{c}·log(30); the data bounds any such c by about 10⁻³. Even Gauss-correlated slot coefficients, which
make the dual coherent for positive sequences (Dunn–Radziwiłł's extremal mechanism, §30c), show no excess. For the
Möbius case this matches the structural reason: μ(n) multiplies the marked sequence b(m) = −μ(m)a(p_m), so the
dual-of-dual sixth-power rows (where all characters are trivial) see Σ μ(m)a(p_m), which cancels.
θ values of these critical-length conjectures (Cycle61 grid model, c61/capc3.py; record-equivalent baseline 0.87492):
  c_I = 1 → 0.87289;  c_P κ = 2, c_I = 2 → 0.87339;  c_P κ = 2, c_I = 1.5 → 0.87215;  c_P κ = 1, c_I = 1 → 0.86979.
These are the cheapest conjectures with a large θ-value. They are critical-length statements (total length ≤ U),
unlike LM/MLV, and the data says they are true.

## 32. Source of c_I = 2 inside the reflected energy (PROVED reading, paper 7657–7900; corrects a first draft)
Lemma reflected-energy:
  E_ref = O/2 + max(H, v + ℓ_b) − S0 − B0 + z_a − u − ℓ_b − 2e_λ/3 − ½(T_d − v − 3ℓ_b − e_λ)_+,
  u = min{v, z_a, (v + z_a)/3},  T_d = 2H + 2A0 + 2z_a − N* − N0 − 3B0  (retained dyads: v + 3ℓ_b + e_λ ≤ T_d).
Two separate slot effects:
(a) Structural, and the binding one. Active slot primes P sit in the reflection conductor, so the maximal dual length
    T_d grows by 2z_a. The saving term is capped, u ≤ z_a. At the maximal dual (v ≈ T_d, so u = z_a) this gives
    E ≈ 2H + 2z_a − N*, an extra +z_a net. This is the second copy of z in r + 2z ≤ 1 (SESSION3 §6.3). It is forced by
    the geometry of Poisson in u: the conductor contains both copies p1, p2, and no saving can exceed the slot length.
(b) Small-dual dyads (2z_a > v): the BGL cubic term gives u = (v + z_a)/3 < z_a, an additional penalty
    (2z_a − v)/3. This is not the binding case.
So within the theta-reflection method, c_I = 2 is intrinsic: a sharper bilinear (n, P) sieve would only remove (b).
Removing (a) needs a non-reflection treatment of the slot primes, for example detecting them after Poisson (as marks
outside the conductor) or by an asymptotic (equidistribution) input as in §25/§33. The numerics (§31) say the true
marked energy has cost 1, so the loss is the method's. Relevant literature: Dunn–Radziwiłł's extremal configurations
concern (b) only. No large-sieve improvement for prime-supported cubic/sextic sequences is known (search, Oct 2026).

## 33. Cost-1 slots from the Möbius fourth moment (PROVED reduction)
For one prime slot P_u = Σ_{p~Z} a(p)ψ_u(p) with Z⁴ ≤ U (z ≤ 1/4), and any N:
  Σ_u |M_u(N)P_u|² ≤ Σ_p|a(p)|²·Σ_u|M_u|² + Z·(U·Σ_u|M_u(N)|⁴)^{1/2}·U^ε.
Proof: write the off-diagonal as T = Σ_{p1≠p2} a(p1)ā(p2) Σ_u χ_{p1}χ̄_{p2}(u)|M_u|². Apply Cauchy over the ≤ Z² pairs and
the large-sieve duality for the set 𝒮 = {χ_{p1}χ̄_{p2}}. Take the u-sum over ALL elements with a smooth weight; this is
allowed by positivity. The hypothesis below is then stated for all elements u ≤ U (rows d⁶v′ are copies, so it follows
from the sixth-power-free version at every smaller family size V = U/q_d⁶, including V < N). Two distinct members of 𝒮 differ by a nonprincipal character of modulus
≤ N(p1p2p3p4) ≤ Z⁴ ≤ U^{1−ε}, whose smooth lattice sum over u ~ U is O(U^{−A}) (Poisson, §25). So the dual Gram is
U·I + O(U^{−A}). If only sixth-power-free u were allowed, the off-diagonal Gram entries would be O(U^{1/6}Z^{10/3}),
which forces z ≤ 5/32.
Consequence: if the Möbius fourth moment is sharp (Σ|M_u(U^r)|⁴ ≪ U^{1+2r+ε}), the marked moment holds with slot cost 1
for slots of total length ≤ 1/4: count 1 − δr − δz.
With loss λ (Σ|M|⁴ ≪ U^{1+λ+2r}) the count is 1 + λ/2 − δ(r + z), z ≤ 1/4. At r* = 0.745 this beats A_I(r*) = 0.5935 iff
λ < 0.114. Copies give λ = (1 + 10r)/6 − 2r = 0.408, and the marked M×M route needs r ≤ 1/3.
So partial cost-1 slots are equivalent in strength to a Möbius fourth moment at r ≈ 0.75 with loss < 0.11:
an MLV-type input again.

## 34. LM is the family's weak Mertens conjecture (PROVED reformulation, explicit formula)
For any D, shifting Perron's integral for M_u(D;W) = (1/2πi)∫W̃(s)D^s L(s,ψ_u)^{−1}ds to Re s = −ε, where the
functional equation gives |1/L(−ε+it)| ≍ U^{−1/2−ε}|1/L(1+ε−it)| ≪ U^{−1/2}log U, yields (simple zeros)
  M_u(D;W) = Σ_ρ W̃(ρ) D^ρ / L′(ρ, ψ_u) + O(D^{−ε}U^{−1/2}log U),
with ρ over all nontrivial zeros of L(s, ψ_u). Hence
  Σ_u |M_u(D)|² = D·Σ_u |Σ_ρ W̃(ρ) D^{ρ−1/2}/L′(ρ)|².
LM (α′ = 0) therefore needs (i) no zeros with β − 1/2 ≫ 0 in bulk (DH-type), and (ii) the family-averaged negative
moment Σ_u Σ_{ρ_u} |W̃(ρ)|²/|L′(ρ, ψ_u)|² ≪ U^{1+ε} (Gonek–Hejhal J_{−1} on average), plus off-diagonal decay of
D^{i(γ−γ′)}. For ζ the analogue is the weak Mertens conjecture ∫^X |M(x)|²x^{−2}dx ≪ log X, known only under RH together
with J_{−1}(T) ≪ T (Ng). Under GRH for the family, LM follows directly from Soundararajan-type bounds
M_u(D) ≪ D^{1/2+ε}. This explains structurally why Poisson/large-sieve methods stall at the copies slope: they see
neither small values of L′ at zeros nor zero locations.

## 35. θ under combined conjectures (MODEL, c61/capc4.py; Cycle61 grid, baseline 0.87492 ≈ record)
| slot costs \ LM slope α′ | 5/6 (now) | 0.70 | 0.50 | 0.30 | 0 |
|---|---|---|---|---|---|
| c_Pκ = 4.5, c_I = 2 (now) | 0.87492 | 0.87372 | 0.87055 | 0.86844 | 0.86667 |
| c_Pκ = 2, c_I = 2 | 0.87339 | 0.87267 | 0.87040 | 0.86847 | 0.86667 |
| c_Pκ = 4.5, c_I = 1 | 0.87289 | 0.87150 | 0.86922 | 0.86793 | 0.86667 |
| c_Pκ = 2, c_I = 1 | 0.87164 | 0.87082 | 0.86904 | 0.86780 | 0.86667 |
LM at α′ = 0 alone gives exactly 13/15 whatever the slot costs. The slot conjectures (numerically supported in the tested ranges, §31, §44) are worth up
to −0.0033 on their own, but they cannot reach 13/15: at α′ = 5/6 the floor with natural slots is ≈ 0.8698 (§13).

## 36. LM is harder than its ζ-analogue, which is open (CITED + PROVED comparison)
For ζ in the t-aspect, LM reads ∫_T^{2T}|Σ_{n~D}μ(n)n^{−it}|²dt ≪ T^{1+ε}D for D ∈ (T, T^{3/2}]. By Gallagher's lemma this is
the variance of μ in short intervals of relative length 1/T (H = D/T) with square-root cancellation on average.
Unconditionally only logarithmic savings are known there:
- Ramachandra (1976): V_μ(X,H) ≪ H²(log X)^{−A} + X^{1/6+ε};
- Matomäki–Radziwiłł (2016): o(H²).
See Gorodetsky's thesis §1, people.maths.ox.ac.uk/gorodetsky/PhD.pdf. The ζ family has no sixth-power copies, so its
current LM slope is 1. The sextic family's slope 5/6 (copies) is already ahead of ζ, and any α′ < 5/6 would be a
power-saving Möbius-variance result of a kind not known for ζ either. With §34 (weak Mertens) and §29 (quadratic
analogue open even under GRH for the asymptotic), LM sits at the level of the classical open problems on Möbius
randomness. It is the single input that gives exactly 13/15 (§35), and the numerical tests (§28, §41) are consistent with it, within the modest range of §44.

## 37. One constant behind two levers: copy efficiency 1/6 (PROVED structural observation)
- LM slope: averaging a row u over its P sixth-power copies u a⁶ (q_a ≤ P) costs a row extension P⁶. Net gain per unit of
  extension is 1/6, hence α′ = 1 − 1/6 = 5/6 (§22, paper Lemma inverse-amplification).
- Plain slot cost: the plain lemma pays the child excess by pointwise slot removal, financed by averaging the Gauss row
  over h ↦ hp⁶ (pool Z^ℓ). That costs a row extension 6ℓ per allowance ℓ, the same 1/6 efficiency, hence
  c_P = 6 (×κ from the pointwise cap) (SESSION3 §6.4; paper 14180–14300).
- Sixth powers are the only multipliers fixing every column Gauss sum. For k < 6, u ↦ up^k twists the row by the small-
  conductor character χ_•(p)^k, which gives new rows, not copies. So within linear row amplification, 1/6 is forced.
- c_I = 2 has a different source: the BGL cubic term in the (dual theta ideal) × (active slot) bilinear form (§32).
So a better copy mechanism, i.e. any way to replicate a row with efficiency > 1/6 (for instance through metaplectic Hecke
operators acting with eigenvalues rather than by exact invariance), would improve both α′ and c_P at once. With
efficiency e: α′ = 1 − e and c_P = 1/e.
Remark (MODEL, speculative). In an n-th order theta version of the architecture, the Poisson frequencies are ua^{2n} and the
copy efficiency would be 1/(2n). n = 2 has trivial Gauss sums (no metaplectic phases for the probe), so n = 3 (rows
ua⁶, efficiency 1/6) is the smallest workable order. The quartic analogue (ua⁸, efficiency 1/8) would be worse on both
α′ and c_P. The sextic family is plausibly the optimal family for this architecture.

## 38. Cross-scale correlation of Möbius large values (NUMERICAL, mmnum/xscale.py; U = 3·10⁴, 108 606 rows)
Correlation of log|M_u(U^r)|² across r ∈ {0.5, 0.745, 1.0, 1.162}: 0.07–0.10 for every pair, the same for μ and for the
random multiplicative control. The tails are much more strongly linked:
- rows in the top 1% at both r* = 0.745 and t* = 1.162: 84, against 10.9 if independent (×7.7; control ×9.0);
- their |M(U¹)|² averages 3.2× the mean (control 3.6×), and only 7% (control 9%) are below the median at r = 1.
So the adversary's pattern (saturated at r*, t*, no saturation at 1) is atypical for typical large values: tails at
different scales share the small-prime (multiplicative-chaos) factor.
- In the explicit-formula picture (§34), a dip at D = U between saturated scales needs ≥ 2 zeros near the bin's maximal real
  part, with comparable residues W̃(ρ)/L′(ρ) and height spacing tuned to ≍ 1/log U.
- A single dominant zero would make M saturated at every scale, including r = 1, which is the DH cost A_I(1) = 1 − δ.
- So the non-DH bad rows are exactly the multi-zero rows. An input bounding rows with ≥ 2 comparable off-line zeros by
  U^{1−δ} (a pair-density statement) plus control of the far-zero tail Σ_j c_jD^{ρ_j} would give DH.
- Neither is available: the tail involves 1/L′ at zeros left of the bin (the same negative-moment issue as §34).

## 39. The core problem, unified, and its numerics (PROVED reductions §33/§36 + NUMERICAL, mmnum/lm2.py)
Every lever in the architecture reduces to beyond-conductor mean squares of inverse-power coefficients:
  (LM_k)  Σ_u |Σ_{n~D} a_k(n)ψ_u(n)|² ≪ U^{1+ε}D,  D ∈ (U, U^{3/2}],  a_k = coefficients of L^{−k}.
- LM_1 at slope 0 gives 13/15 exactly (§35).
- LM_2 at D = U^{2r}, r ≈ 0.75 (equivalently the Möbius fourth moment at r) with loss < 0.11 gives slot cost < 2 (§33).
- MLV at r* is a large-values form of LM_1 below the conductor, enhanced by slots.
Known: slope 5/6 for every k (copies). Numerics, ratio of mean to diagonal for a_2 = μ∗μ:
| U | s = 0.5 | 0.75 | 1.0 | 1.1 | 1.2 | 1.3 | 1.4 | 1.5 |
|---|---|---|---|---|---|---|---|---|
| 3·10³ | 0.903 | 0.878 | 0.877 | 0.879 | 0.860 | 0.845 | 0.869 | 0.844 |
| 10⁴ | 0.910 | 0.891 | 0.866 | 0.870 | 0.860 | 0.849 | – | – |
The random multiplicative control (ε∗ε-type) gives 0.87–0.91. The data are flat, consistent with LM_2 (same §44 caveat on range).
At U = 3·10⁴ (4000 sampled rows, mmk2), LM_2: 0.855, 0.856, 0.859, 0.854, 0.829 (±0.02) for s = 1.0 … 1.4. This is
flat. Effective columns/rows reaches only 1.18 at s = 1.4.

## 40. Nazarov–Turán + local Jensen: the dip at r = 1 and what it costs (PROVED reduction with explicit constants)
The adversary needs M_u unsaturated on the whole interval r ∈ (r*, t*) = (0.745, 1.162) (§12: a saturated witness at
r ∈ (1 − 0.255, 1] already costs < R), yet saturated at both ends. Write f(ξ) = e^{−aξ}M_u(e^ξ; W). Shifting the Mellin
contour to a zero-avoiding line Re s = c ∈ [a − x − w, a − x] (c chosen by averaging Σ_j log(1/|c − β_j|)) gives
  f(ξ) = Σ_{j≤n} c_j e^{(ρ_j − a)ξ} + E(ξ),  c_j = W̃(ρ_j)/L′(ρ_j),
where n is the number of zeros in a disk of radius ≍ 4x around s0 = a + 6e + iγ.
- Jensen with |L(s0)| ≥ U^{−ε} (bin) and max_disk |L| ≤ U^{μ0} gives n ≤ μ0 log U/log 2.
- Titchmarsh's local factorisation plus the averaged choice of c gives |E| ≤ U^{−xr}·U^{μ0(1.44(log(1/w)+1) + C_T)}.
- Nazarov's Turán lemma on I ⊃ E_dip (|I|/|E_dip| ≈ 1): sup_{E_dip}|f| ≥ U^{−0.42(x+w)}·C^{−n}·sup_I|f|.
For the dip to be impossible at the depth that matters (|M(U)| within U^{−0.03} of saturation, which would give count
A_I(1) + 0.06 < R = 0.5935) we need n·log C ≤ 0.03 log U, i.e. n ≲ 0.02 log U, and |E| ≪ U^{−0.03}. Hence:
  DH at the binding class ⟸ for bad rows, max|L| ≤ U^{μ0} on a disk of radius ≈ 0.25(a − 1/2) around the bin, with μ0 ≲ 0.014.
- Available: convexity/Phragmén–Lindelöf from the bin gives μ0 ≈ (radius)/2 ≈ 0.03–0.06. Excluding rows by the plain fourth
  moment at Re s > 1/2 gives μ0 = δ/4 = 0.117 at DH quality.
- Reaching μ0 = 0.014 at DH quality by moments needs Lindelöf-on-average of order 2k with k ≈ δ/(2μ0) ≈ 17.
So the explicit-formula/dip route shows, modulo the remainder E (whose control is itself the obstruction), why
single-zero rows would be DH-cheap (n = 1: no dip possible) and why the adversary needs ≍ log U zeros near the bin.
Closing it needs high moments, as in the classical LH ⟹ DH.
Addendum to §40 (MODEL exponent algebra). Off-line moments help, but not enough. At σ1 ≈ 0.70, the approximate functional
equation for L(σ1)^k (main part with weights n^{−σ1}, handled by copies up to length U^{ν}, ν ≤ (5/6)/(11/6 − 2σ1) ≈ 1.92;
dual part damped by U^{k(1/2−σ1)}) gives Σ_u|L(σ1,ψ_u)|^{2k} ≪ U^{1+ε} for k ≤ 1 + 1.48ν ≈ 3.8, i.e. a 7.7-th moment.
- Many-zero rows (n ≥ c0 log U in a disk of radius ≈ 0.03) have max_disk|L| ≥ U^{0.69c0}, so they number ≤ U^{1−5.2c0}.
  This is DH quality iff c0 ≥ 0.09.
- Few-zero rows are dip-free iff c0·log C ≤ 0.03, with Nazarov's constant log C ≈ 3, i.e. c0 ≤ 0.01.
The gap is a factor ≈ 9, so the route needs a ≈ 20th-order off-line moment (or a 9× better zero-cluster bound).

## 41. LM at larger scales and the dual side (NUMERICAL; mmnum/mmk2.c table-free kernel, lm_big.py, dualB.py)
Kernel mmk2 computes χ_π(u) by direct powering (validated: identical output to mmk at U = 10³, s ≤ 1.5).
Rows are a uniform random sample of the sixth-power-free rows; values are mean|M|²/diag ± standard error.
| U (sample) | s = 1.0 | 1.1 | 1.2 | 1.3 | 1.4 | 1.5 | cols_eff/rows at s_max |
|---|---|---|---|---|---|---|---|
| 10⁴ (6000) | 0.919±.013 | 0.952±.013 | 0.950±.013 | 0.930±.013 | 0.924±.013 | 0.930±.013 | 1.92 |
| 3·10⁴ (4000) | 0.938±.015 | 0.923±.016 | 0.934±.016 | 0.953±.017 | 0.913±.016 | – | 1.18 |
| 10⁵ (4000) | 0.934±.016 | 0.937±.016 | 0.914±.015 | – | – | – | 0.19 |
| 3·10⁵ (4000) | 0.927±.016 | 0.941±.016 | – | – | – | – | 0.07 |
| 10⁶ (4000) | 0.921±.016 | 0.944±.017 | – | – | – | – | 0.08 |
Random multiplicative control: 0.92–0.97. Max/mean ≤ 19 throughout, so there are no large-value outliers.
Robustness to the weight (lm_big_w2.py): a flat-top bump exp(1 − 1/(1 − s⁸)), truncated sharply at N(n) = D/2, gives
0.933, 0.947, 0.948, 0.925, 0.928, 0.928 at U = 10⁴, s = 1.0, …, 1.5. The flatness is not an artefact of the
smooth weight.
Dual side (dualB.py, U = 10³, s = 1.3, D = 7943, all 228 900 elements h with N(h) ≤ K = D²/U):
B_h(D) = Σ μ(n)g̃(n)χ_n(h)W(N(n)/D) has mean|B_h|²/diag = 0.941, max 10.7, and the same mean on N(h) ≤ K^{1/4} (0.950),
(K^{1/4}, K^{1/2}] (0.917) and (K^{1/2}, K] (0.941). There is no bias at small h and no structured large values. The
involution's dual is as random as the primal, so the off-diagonal E is invisible in both representations, as LM predicts.

## 42. Slot test at the binding Möbius length r* = 0.745 (NUMERICAL, mmnum/slotbig.py; 20 000 sampled rows each)
ρ = mean|M_u(U^{0.745})P_u(U^z)|²/(mean|M|²·mean|P|²) with slots beyond the paper's range (r + 2z > 1, r + z ≤ 1):
| U | z = 0.20 (r+2z = 1.145): a ≡ 1 / random | z = 0.25 (r+2z = 1.245): a ≡ 1 / random |
|---|---|---|
| 10⁷ | 1.043±.011 / 0.992±.011 (5 primes) | 1.067±.013 / 1.009±.012 (6 primes) |
| 10⁸ | 1.051±.013 / 0.974±.011 (5) | 1.074±.015 / 0.996±.012 (10) |
| 10⁹ | 1.051±.013 / 1.005±.012 (6) | 1.064±.015 / 0.990±.013 (17) |
- Random-phase slots: ρ = 1 within error. No penalty at all.
- For a ≡ 1 there is a stable +4–7% excess that does not grow from 10⁷ to 10⁹. It is explained exactly by the Möbius
  pair terms n1 = kp2, n2 = kp1, where μ(kp1)μ(kp2) = μ(k)² > 0. These give relative size ≈ (n_P − 1)/Z
  ≈ 1/(2 log Z) (n_P = number of slot prime ideals): predicted 0.05–0.09 for the rows above.
This is a logarithmic correlation, not a power of U. A genuine cost-2 obstruction at r + 2z = 1.245 would show as growth
∝ U^{c·0.245}, and none is seen across a factor 100 in U.

## 43. What natural slot costs would and would not do (MODEL, c61/capc3.py)
With c_I = 1 and full saturation (p = ½, 2q = δ), the marked count is A_I(r) = 1 − δr − δ·min(1 − r, ℓ/d), where ℓ/d ≈ 0.26 is
the architecture's slot supply (ℓ = 11/3 − 4θ < 1/5). With c_Pκ = 1, A_P(m) = 1 − 2δm − δ·min(1 − 2m, ℓ/d).
- Both equal the DH value 1 − δ once r ≥ 1 − ℓ/d, resp. m ≥ (1 − ℓ/d)/2. The adversary can then avoid DH only at block
  lengths t < (1 − ℓ/d)·3/2 ≈ 1.11, so R − (1 − δ) ≈ (5/6 − δ)(0.11) ≈ 0.02-0.03. Measured at the binding:
  δ = 0.536, R = 0.4869 against DH 0.4638, giving θ = 0.86979.
- Removing the supply cap altogether (ℓ/d = ∞) gives θ = 0.86921, not 13/15. At p = ½ the count is then exactly DH, and
  the binding moves to partially saturated bins p < ½. There the slot gain 2δp(1 − r) < δ(1 − r), and the long witness
  (slope 5/6) is the bottleneck again.
So even an optimal slot theory leaves the p < ½ bins needing LM. That confirms LM as the single essential input for 13/15.


## 44. Correction: how far beyond the conductor the LM numerics really reach (NUMERICAL bookkeeping)
The tables in §28, §39 and §41 originally compared with "copies would allow ×U^{5(s−1)/6}". That is an asymptotic exponent
with unknown constants, so the comparison overstated the evidence. The honest finite-U benchmark: the family has
rows ≈ 3.62·U (sixth-power-free elements with units), while the Möbius polynomial has effective column count
cols_eff = Σ_n μ²(n)W(N(n)/D)² (a small multiple of U^s; the weight W is narrow), measured directly as cols_eff/rows:
  U = 3·10³: 0.02 (s = 1), 0.09 (1.2), 0.47 (1.4), 1.05 (1.5);  U = 10⁴: 0.30 (1.3), 0.76 (1.4), 1.92 (1.5);
  U = 3·10⁴: 1.18 (1.4);  U = 10⁵: 0.19 (1.2);  U ≥ 3·10⁵: ≤ 0.08.
Even a sharp large sieve allows a mean square up to (1 + cols_eff/rows)·diag. So only points with cols_eff/rows ≳ 1
test LM genuinely: U = 10⁴, s = 1.5 (sharp-sieve room ×2.9, observed 0.930); U = 3·10⁴, s = 1.4 (room ×2.2, observed 0.913);
U = 3·10³, s = 1.5 (room ×2.05, observed 0.928). The other points are rows-dominant (the OpenAI raw-moment regime in
counts), where ratio ≈ 1 is expected anyway.
So the numerical support for LM is real but modest: no off-diagonal mass at up to twice the effective conductor.
The slot tests (§31, §42) are likewise in the rows-dominant regime for the product M·P. There they show that no
Gauss-sum coherence (Dunn–Radziwiłł mechanism) creates excess, which is the relevant point for c_I, but they cannot
by themselves distinguish cost 1 from cost 2 asymptotically.

## 45. LM far beyond the conductor (NUMERICAL, mmnum/lm_beyond.py; mmk2, all 3636 rows at U = 10³)
This addresses §44: points with cols_eff/rows ≫ 1, where the sharp large sieve allows (1 + cols_eff/rows)·diag.
| s | 1.0 | 1.3 | 1.6 | 1.8 | 2.0 | 2.2 |
|---|---|---|---|---|---|---|
| cols_eff/rows | 0.02 | 0.15 | 1.20 | 4.79 | 19.1 | 76.0 |
| μ: mean\|M\|²/diag | 0.930±.016 | 0.934±.016 | 0.954±.017 | 0.945±.017 | 0.893±.016 | 0.960±.018 |
| random multiplicative | 0.941 | 0.927 | 0.949 | 0.979 | 0.969 | 0.951 |
| squarefree indicator | 0.95 | 1.09 | 2.5 | 7.7 | 29 | 115 |
For μ, max/mean ≤ 12.3 at every s. The squarefree control grows like cols_eff/rows, but this comes from a single row: max/mean
≈ 3600 ≈ number of rows, so it is the unit-type row with trivial character. The plain sum there is ≈ D, as it should be.
**Reading.** At fixed U and s → ∞, Proposition C′ predicts mean|M|²/D → avg_u Σ_ρ|W̃(ρ)|²/|L′(ρ)|². This is a constant if the
family has no off-line zeros and has Gonek–Hejhal-size J_{−1}. The data show exactly that plateau, at 0.89–0.96 up to 76×
past the effective conductor. So the beyond-conductor flatness is the GRH/J_{−1} behaviour of the 3636 individual
L-functions. It supports LM, but it cannot separate LM from GRH-on-average, as no finite computation can.
Addendum, U = 3·10³ (3000 of 10 872 rows sampled; the unit-type row is not in the sample):
| s | 1.5 | 1.7 | 1.9 | 2.0 |
|---|---|---|---|---|
| cols_eff/rows | 1.05 | 5.20 | 25.8 | 57.4 |
| μ | 0.930±.019 | 0.939±.018 | 0.878±.017 | 0.963±.020 |
| random multiplicative | 0.963 | 0.980 | 0.953 | 0.979 |
| squarefree indicator | 0.721 | 0.499 | 0.253 | 0.145 |
With the trivial row absent, the squarefree control falls steadily below the diagonal. That is the expected behaviour:
μ² = Σ_{d²|n}μ(d) makes the sum mostly smooth, and smooth sums vanish beyond the conductor. So the test separates
coefficient types clearly, and μ stays on the random-multiplicative plateau at both U.

## 46. Formulation check: LM has to be stated for sixth-power-free rows (PROVED, short argument)
Theorem C's sum runs over all elements. Over all elements u ≍ V, the rows εa⁶ (ε a unit, N(a) ≍ V^{1/6}) are copies of the
unit-type rows: ψ_{εa⁶}(n) = ψ_ε(n)·1_{(n,a)=1}, where ψ_ε has conductor O(1). Averaging over a removes the coprimality, and
positivity gives
  V^{1/6}·|M_ε(D)|² ≲ MS_all(V; D).
- The critical-length all-elements second moment (V = D^{1+ε}, which the copies step of the record already uses) gives
  |M_ε(D)| ≪ D^{11/12+ε}, i.e. the zero-free exponent 1 − 1/(2·6). This is the architecture's base exponent 11/12.
- An all-elements bound MS_all(V; V^s) ≪ V^{1+α′(s−1)+ε}V^s gives |M_ε(D)|² ≪ D^{1+(5+6α′(s−1))/(6s)}, i.e.
  zero-free for Re s > 1/2 + (5 + 6α′(s−1))/(12s). Here α′ = 5/6 returns 11/12, which is the copies fixed point again.
  All-elements LM (α′ = 0, s = 3/2) would give 7/9 for these L-functions directly, and α′ < 0.58 at s = 3/2 would already
  beat 0.8748.
- So all-elements LM contains a quasi-GRH statement for the small-conductor members. It is circular as a route, since its
  principal-type part is as large as D^{2θ_true}·V^{1/6}.
The LM of DH_NOTE §1/§11 is over sixth-power-free rows with N(u) ≍ U and s ≤ 3/2. For a fixed row it does not reach D → ∞,
so it implies no zero-free region for any single L-function and is strictly weaker. That is the version the bins need,
and the version tested numerically.

## 47. The dip as a divided difference (PROVED identity; HEURISTIC consequences, dip/prog.py)
Take a cluster ρ_1, …, ρ_n of zeros of L(s, ψ_u) and write L(s) = Π_j(s − ρ_j)·g(s). Then L′(ρ_j) = g(ρ_j)Π_{k≠j}(ρ_j − ρ_k), so
the cluster's part of f(ξ) = e^{−aξ}M_u(e^ξ; W) is exactly a divided difference:
  f_cl(ξ) = Σ_j F(ρ_j)/Π_{k≠j}(ρ_j − ρ_k) = F[ρ_1, …, ρ_n],  F(s) = W̃(s)e^{(s−a)ξ}/g(s).
Hermite–Genocchi (valid for complex nodes) gives f_cl = ∫_{simplex} F^{(n−1)}(Σ t_jρ_j)dt. If h = W̃/g has nearly constant phase on
the hull of the cluster, this becomes
  f_cl(ξ) ≈ h·ξ^{n−1}/(n−1)!·ν̂(ξ),  ν̂(ξ) = ∫ e^{(w−a)ξ}dν(w),
where ν ≥ 0 is the push-forward of Lebesgue measure on the simplex.
Consequences:
- With ν ≥ 0, a pure cluster is a Laplace transform of a positive measure. Its zeros near the real ξ-axis have density at most
  about X/(2π), where X is the imaginary width of the hull (Pólya; the bound is attained by a vertical progression, for
  which |f_cl| ∝ |2 sin(Δξ/2)/Δ|^{n−1}). A dip of depth U^{−0.03} over the window r ∈ [0.88, 1.08] needs ≳ 0.015 log U zeros of f
  near the window (slope argument: log|f| changes slope by +0.3 across the bottom). So the relevant zeros
  (Re ρ ≥ a − 0.03) must spread over heights X ≳ 0.2–0.5. §40's local disk of radius ≈ 0.06 cannot dip on its own. The
  dip needs interference between zeros at heights ≥ 0.2 apart.
- The vertical progression (spacing 15/log U, n = c0 log U, X = 15c0) gives the profile η(r) = −c0 log|cos(π(r − r*)/0.417)|.
  It is quadratic at the saturated ends, while the cost profile needs linear decay (0.1165 per unit r after r*). So with
  exact saturation at r*, t* the adversary loses for every c0 (prog.py: it would need its maxima moved 0.2–0.3 outside
  [r*, t*], which breaks the covering condition). The adversary therefore needs a non-uniform zero pattern with sharp edges.
- Localising in height is convolution in ξ. A weight whose Mellin transform separates clusters at height distance X_sep
  to precision U^{−κ} has log-support ≈ κ log U/X_sep. It costs ≈ (1 + α′)κ/X_sep in the count (longer polynomial). With
  κ = 0.03 and α′ = 5/6 this beats R at the binding only if X_sep ≳ 0.9.
Why this does not close the route:
(i) positivity needs arg(W̃/g) nearly constant on the hull, i.e. |L′/L| ≲ 1/X off the cluster. That is a large-values
condition on L′/L, the same kind of input as §40.
(ii) the adversary can spread its relevant zeros continuously (no gaps of size X_sep), so localisation cannot isolate a
sub-cluster without tapered edges, and tapering destroys positivity.
Net: the dip needs many off-line zeros spread over heights ≳ 0.2 with no large gaps, plus control of L′/L. Bounding
such rows is a multi-zero (pair-correlation-type) density statement, consistent with §38 and §40.

## 48. Literature update, 16 Oct 2026 (CITED; downloaded HTML read as data, scratchpad/lit2610/)
- **Xu–Sang, arXiv 2610.15291 (v2, 15 Oct), B1 = 0.8748867.** Its levers are: the low estimate keeps the rescaled-tuple
  decay (ℓ up to 1/5); an amplified first moment of the reflected numerator (Prop. 5.1: plain fourth moment with prime
  slots plus Hölder, giving K = (3R+1)/4 and contour shift ϑ = q/(6κ)); and κ-feedback (Lemma 2.1, κ ≥ 13/18).
  Its endpoint (28) is E = C0 + δ/2 + q{ℓ + c(h/4 − l_y/2)} − (3h/4)(1 − R), with C0 = −1/4 + b/6 + 5ℓ/4 = −x/3 − 1/12 + 13ℓ/12.
  This is exactly the high side already in our Cycle61-type model (c61/c61model.py: d(1/4 + 3R/4 − 17/50) − (q/2)C_d,
  with C_d = min(d, 2y−d)/(6κ)). So it is not a new lever relative to our pipeline, and 0.8748867 > 0.8748389.
  Its row count R_κ(δ,x) = 1 − δ + (α−δ)δP/(2J) is the closed form of the OpenAI count; there is no DH or LM input.
- Its ceiling is the same: B(ℓ) = 11/12 − ℓ/4 with ℓ ≤ 1/5, and C0 = b/6 > 0 at ℓ = 1/5. With R = 1 − δ + λ, the λ-tolerance at
  a given ℓ is λ ≲ 0.3(1 − 5ℓ), which is essentially Conditional Theorem D. DH is still what 13/15 needs.
- **Liu, arXiv 2610.12234: B_L = 0.874957.** Parameter optimisation of (b, ℓ, κ) with an optimality proof inside that system.
  No new input.
- **Xu, companion manuscript [Xu26], as described in 2610.15291 §1.1.** A lower bound ≫ A^{5/6}B^{1/3} for the squared norm of
  the sextic residue-character matrix with squarefree indices, for B^{16/11+δ} ≤ A ≤ B^{2−δ}. In LM terms (A = D = U^s,
  B = U, s ∈ [16/11, 3/2]) the general-coefficient large-sieve constant is Δ ≥ D^{5/6}U^{1/3}. That is not sharp: it
  exceeds max(U, D) = D by (U²/D)^{1/6}. It also exceeds what copies give for μ: Σ_u|M_u|² ≪ U^{1/6+5s/6}·U^ε·D, compared
  with Δ·D ≥ U^{1/3+5s/6}·D, a gap of U^{1/6}. So even the existing slope 5/6 is Möbius-specific (it rests on the
  critical-length Möbius moment) and cannot come from operator-norm bounds. Any improvement must also use the arithmetic
  of μ. The exceptional direction is the Gauss-sum vector itself (h = k⁶ in the dual, coefficient
  B_1(D) = Σ μ(n)g̃(n)W). For μ it is harmless if |B_1(D)| ≪ D^{5/6}U^{1/12}, and §30(c) observes |B_1| ≍ D^{1/2}.
  Removing it leaves the dual asymptotic with relative error U/D, which is LM itself (Theorem C).
- No 2026 result on density hypotheses for thin families (quadratic, cubic, sextic) or on Möbius mean squares beyond the
  conductor turned up. Gao–Zhao (2212.09241) is still the closest, and it is below the conductor and under GRH.

## 49. The DH gap in closed form (CITED Xu–Sang (21); identity with our model PROVED numerically to 1e−16)
Xu–Sang's closed form of the OpenAI row count equals our c61/fastR.py R(δ, p, κ, 'old') with x = p, to machine precision,
at 40 test points across δ ∈ [0.05, 0.74] and p ∈ [0, 1/2]:
  R = 1 − δ + λ,  λ(δ,x) = (α − δ)·δ·P(x)/(2J(δ,x)),  P = (2 − 2cx)(1 − x),  J = (α − δ)(3 − (1+2c)x) + δP,
with α the long-witness (LM) slope and c = 1/(3κ) = 2/c_P (c_P the plain slot cost). Consequences:
- **The DH gap is proportional to (α − δ).** For c_P > 1, P > 0 on x ∈ [0, 1/2], so the architecture reaches DH at bin δ iff
  δ ≥ α. With α = 5/6 this holds only for σ ≥ 11/12, outside the needed range. With an LM slope α′ it holds for
  σ ≥ (1 + α′)/2. So a slope α′ gives a **density hypothesis on [(1+α′)/2, 1]** for the sextic family, and the target range
  (1/2, 0.875] needs α′ → 0 (α′ ≤ 0.02 already for σ ≥ 0.51). That is LM in substance.
- Slot costs enter only through P and J. Natural plain cost c_P = 1 (c = 2) gives P(1/2) = 0, so DH holds at fully
  saturated slots (x = 1/2) for every δ, as in §43. The p < 1/2 classes still carry λ ∝ (α − δ).
- Values at the record (α = 5/6, κ = 0.74968): λ = 0.0596 at the binding (δ, x) = (0.466, 1/2); max over x = 0 is 0.084
  near δ ≈ 0.45; λ ≈ δ/3 for small δ.
- Since λ ≤ δ·P/(2D) ≤ δ/3, bins with δ ≤ 3λ_target are automatic. Every other bin needs α′ ≲ δ. That is why θ moves
  continuously with α′ (§35) and reaches 13/15 only as α′ → 0.
So the whole DH problem for this architecture is the single parameter α (long Möbius slope), apart from the P-weight of
the slot terms. This is the sharpest form of the reduction: LM slope α′ ⟹ DH on σ ≥ (1+α′)/2 within the architecture (given the other
inputs). The converse is not claimed.

## 50. Correction: the bin floor a0 = 51/100 caps DH-based gains at 0.86979 (MODEL, sextic/floor/floor.py; paper text)
The buffered-bin lemma starts at a0 = 51/100 (dl_qrh2/p.txt l.1082; region D1: sr ≥ 51/100 in (4.11), where the Euler
correction H_{η,u} is holomorphic). The floor class a = a0 (δ0 = 1/50) is counted trivially (R = 1) and contributes
E(h) ≤ C0 + (1/2 + 3ℓ/2)δ0 (paper: −m_F = C0 + (1/2+3ℓ/2)/50). In Conditional Theorem D's language this is λ = δ0 = 0.02
at the floor, whatever the counts above it. OpenAI high side, closed-form counts below the DH range, trivial floor
(floor.py; the no-DH run reproduces Liu's 0.874957 as a check):
| floor σ | DH assumed from σ | θ |
|---|---|---|
| 0.51 (paper) | none | 0.874957 |
| 0.51 | 0.51 | 0.869792 |
| 0.505 | 0.51 or 0.505 | 0.868248 |
| 0.5025 | 0.51 | 0.867462 |
| 0.501 | 0.51 | 0.867269 |
| 0.501 | 0.501 | 0.866986 |
So **DH on [0.51, 0.867] alone gives 0.86979, not 13/15.** 13/15 needs the floor a0 → 1/2 together with DH down to
1/2 + O(λ). Bins with δ ≲ 3λ are already covered by the existing count (λ ≈ δ/3 there, §49). Earlier statements of the
target window "σ ∈ [0.51, 0.867]" (DH_NOTE §0, DH_SUMMARY §1, task #37 title) are corrected accordingly.
Lowering the floor needs the correction region D1 (sr ≥ 51/100, zr ≥ 17/50, wr ≥ −1/100, sr + wr ≥ 1 + ε0) to extend
to sr ≥ 1/2 + ε. The Euler correction is a ratio of local factors whose product converges absolutely for sr > 1/2, so
this is plausible with ε-dependent constants. Partial check (Liu 2610.12234, §34.1 local table, dl_qrh2/p.txt
l.3750): for p ∤ u the local error exponents are −sr−ω, −sr−6ζ, −ω−6ζ, 4−6sr−6ζ, 1−sr−ω−6ζ, and all must be < −1. At
ζ = 17/50 the only one involving sr alone is 4 − 6sr − 2.04 < −1, i.e. sr > 0.4933. The central contours are Re s = a+16e and
Re w = 1−a−6e, so sr + wr = 1 + 10e for every a. Nothing in this table needs sr ≥ 51/100: the floor is a margin choice,
and a0 = 1/2 + ε is compatible with it. The other occurrences of 51/100 are checked in §53; all are margins. At the current record the floor is not binding
(slack ≈ 0.01 in E), so this only matters near 13/15. The same caveat applies to "LM ⟹ 13/15" (§35, LM_PROBLEM.md):
the capc4 grid starts at δ = 0, so it implicitly assumes the floor at 1/2.

## 51. Row scales below 1/2: needs the record's counts, not Liu's cruder one (MODEL, sextic/floor/scales.py)
OpenAI/Liu's exponent at row scale U = Z^d is E(d) = a − B + h(17/50 − 1/6) − a·l_y − ℓ/2 + qℓ + d(R_d + δ/2 − 17/50).
The intermediate norms d ∈ [d_min, 1/2] use the unselected count R = 76/75 − 2δ/3 (Liu l.3421ff; affine bound
[(225ℓ + 50 − 75b)δ + 78ℓ − 49b − 73]/300). At the 13/15 parameters (ℓ = 1/5, b → 0, δ up to κ = 0.733) this equals
+0.041 > 0, so these rows fail. REPORT §3's Conditional Theorem D assumed DH only for d ∈ [1/2, h+ζ]:
| floor σ | DH at d = h only (Liu's crude count for d ≤ 1/2) | DH at all d |
|---|---|---|
| 0.51 | 0.872248 | 0.869792 |
| 0.501 | 0.870121 | 0.866986 |
With DH at d = 1/2 as well, E(1/2) = −0.198 + 0.15δ ≤ −0.088. So the correct hypothesis is the DH count for the sextic family
at every scale U (which is how DH is naturally stated), with the floor a0 → 1/2 (§50). REPORT.md and DH_NOTE §1 are
amended. LM at every scale gives DH at every scale, so "LM ⟹ 13/15 (with a0 → 1/2)" stands.
**Softening (re-check, same script).** The failure above needs Liu's crude count R = 76/75 − 2δ/3 at d ≤ 1/2. The record's own
detector count (the closed form of §49, a theorem at every scale) is much better there. With DH only at d = h and the
closed-form count at d = 1/2 (floor 0.501), θ = 0.866986, the same as with DH at all d. The tolerance at d = 1/2 is
(h/d − 1)(0.66 − δ/2) ≥ 0.17, while the existing λ ≤ 0.084. So REPORT's formulation (DH for d ≥ 1/2) is fine in substance,
provided the existing counts are used below. Only the floor correction (§50, §53) is essential.
**Third technical condition (principal Euler region).** The contradiction step needs H_η holomorphic and ≈ 1 on
Re s > B. OpenAI prove this on Re s > 7/8; Xu–Sang extend it (App. C.4, region D2′) to Re s ≥ 87/100, with
|D| ≤ Q^{−87/100}, |V| ≤ Q^{−99/100}, |R| ≤ Q^{−221/100} and H_p − 1 = O(Q^{−181/100}) for p ∤ u. For B = 13/15 = 0.8667 this has
to go slightly below 0.87. The edge exponent −1.81 leaves room of 0.81 against the convergence threshold −1, so a small
extension is plausible. Checked in §52: holomorphic for Re s > 2/3.
**Summary of what 13/15 needs, precisely.**
(i) DH counts for the sextic family at every row scale and every bin σ ∈ (1/2 + O(λ), 0.867];
(ii) the bin floor a0 lowered from 51/100 to 1/2 + O(λ) (local Euler exponents allow sr > 0.4933);
(iii) the principal Euler region extended to Re s > 13/15 (currently 87/100).
(ii) and (iii) look like routine margin changes; (i) is the open problem, and LM at every scale implies it.

## 52. Conditions (ii) and (iii) checked against the exact local Euler formulas (PROVED at table level; Liu (17.4) = OpenAI Lemma 7.1)
Local data (Cycle61 paper "Literal local Euler data"; Liu 2610.12234 proof of Lemma 17.1):
- V = Q^{−6z}, R = a_p²Q^{4−6s−6z}, W = χ_p(u)Q^{−w}, D = η(p)χ̄_p(u)Q^{−s};
- H_p − 1 = [D(V + W − VW) − VW + (1−V)(1−W)F]/(1 − D), with F = P* + D.
For p ∤ u, the exponents occurring are −sr−wr, −sr−6zr, −wr−6zr, 4−6sr−6zr, 1−sr−wr−6zr. Normal convergence of Π_{p∤u}H_p
needs all of them < −1.
**(ii) Floor.** On the central contours Re s = a+16e, Re w = 1−a−6e, Re z = 17/50 at a = 1/2 the exponents are
(−1−10e, −2.54, −2.54, −1.04−96e, −2.04). All are < −1, and the binding one 4−6sr−6zr < −1 holds for a > 0.4933. At a row
prime p | u (W = D = 0) the J_j terms have exponents −1, 1.5−3a, 1−2a, 2−4a, 1.5−3a, 3−6a. These are 0 at a = 1/2 and so
O(1) per row prime, giving a finite product ≪ U^ε. So the Euler correction allows the floor a0 = 1/2 + ε. Not re-derived
here: the paper's weighted-mass accounting at row primes (inequalities B0.42–B0.44). Its displayed j = 2 example,
Q^{1−a}Q^{−3/2}Q^{(2a−1)/2} = Q^{−1}, is independent of a.
**(iii) Principal region.** At the principal residue point (w, z) = (1, 1/6) the exponents are (−1−s, −1−s, −2, 3−6s, −1−s).
They are < −1 for Re s > 2/3, so H_η is holomorphic on Re s > 2/3. At s = 13/15 the maximum is −1.867, so the tail after a
prime cutoff is < 1/2. At the D2 corner (Re w = 19/20, Re z = 33/200) with s = 13/15 the maximum is −1.807. The selected
ratio G_p/H_p = −1 + O(Q^{−min(sr, 1, 5sr−3)}) = −1 + O(Q^{−0.867}). So the principal region extends to Re s ≥ 13/15 with
room to spare. Xu–Sang's 87/100 is a choice, not a limit.
So, at the level of the exact local tables, (ii) and (iii) are margin changes. The conditional statement is:
**DH (count U^{1+o(1)−δ}) for the sextic family at every scale and every σ ∈ (1/2, 0.875] ⟹ θ = 13/15**
(the upper end is β* ≤ 0.8748389, which is known, because in the contradiction step the bins run up to δ ≤ 2β* − 1), up to re-checking the
row-prime mass inequalities at a0 → 1/2. [Resolved by §61: B0.42–B0.44 belong to the Cycle61 plain witness, i.e. to the record's
detected count, which the final form of D′ no longer uses. D′ uses only the DH hypothesis at d ≥ 2/5 and the OpenAI paper's own
counts. The paper's row-prime bounds are covered in §53 (Lemma 7.1, Prop. 16.1).] Item (i), DH itself, is the open input.
Cross-check (scales.py with a uniform λ at all scales, floor 0.50005): θ = 0.866683, 0.869792, 0.872774 for λ = 0, 0.02, 0.04,
against the Conditional Theorem D formula 0.866667, 0.869792, 0.872774. With the corrected hypotheses the formula is exact.

## 53. Every use of the floor 51/100 in the OpenAI paper is a margin (PROVED by inspection; oai_pdf/oai5.txt, public PDF)
Occurrences, with what they need at a0 = 1/2 + ε_f (ε_f > 0 fixed before Z, constants allowed to depend on it):
- **Lemma 7.1, region (7.14).** Good primes: the paper uses 4 − 6x_r − 6z_r ≤ −11/10. At x_r ≥ 1/2 + ε_f it is ≤ −1.04 − 6ε_f.
  For x_r < 51/100, the condition x_r + w_r ≥ 1 + ε0 forces w_r > 0.49, so ϑ = (−w_r)₊ = 0 and E_p ≪ Q^{−1.04}; DW ≤ Q^{−1−ε0}.
  Row primes: the closest exponent 3/2 − 3x_r ≤ −3ε_f gives H_p − 1 = O(Q^{−3ε_f}), so H_p is bounded and the divisor-product
  bound still gives O(q_u^ε). Hence ε_H = min(ε0, 3ε_f) > 0.
- **§8 buffered bins (pigeonhole base) and Prop. 8.3 (detector).**
  - Re(ρ − 1/4) ≥ 1/4 > 0 (paper: 26/100).
  - Gamma-integral exponent −3 + (3/2)(5/4 − σ) ≤ −15/8 (paper: −189/100).
  - "δ ≥ 1/50 forces m ≤ 1/2 + O(ε)" becomes m ≤ 1/2 + O(ε/δ_f), fine since ε is chosen after δ_f = 2ε_f.
- **Prop. 9.2.** The floor bin uses the trivial count R = 1. Nothing else changes.
- **Prop. 16.1 (dynamic local errors), valid for a ≥ 1/2 + ε_f.**
  - Good-prime exponents −a−4e, −51/25+12e, 49/25−5a−68e, a−51/25+12e: the maximum is −a−4e < −1/2. The paper's Q^{−51/100}
    becomes Q^{−1/2−4e}, which still gives the stated P_i^{z_r−1/2+O(e)}.
  - Row-prime exponents (a−2, 1/2−2a, −a, 1−3a, 1/2−2a, 2−5a) are ≤ −1/2 at a = 1/2, and strictly below with the
    e-corrections (−32e, −26e/−48e, −42e, −80e) and ε_f.
  - A* = a − 1/2 + 6e > 0, and the coprime-label exponent z_r − 51/100 becomes z_r − 1/2 − 4e.
- **Lemma 19.1, §20 (11745, 11848: "δ ≥ 1/50").** Only uniformity constants, which now depend on δ_f.
- **(20.5) floor reserve.** E(h) ≤ C0 + (3/4)δ_f. This is the term that has to vanish for 13/15.
Conclusion: the floor a0 = 51/100 is a free parameter, and any a0 > 1/2 works with the same proofs. Together with §52
(principal region holomorphic for Re s > 2/3), the conditional chain is precise:
**DH at every scale on (1/2, 0.875] for the sextic family ⟹ θ = 13/15 in the OpenAI architecture.**
Source check of the endpoint (oai5.txt l.11570ff, OpenAI §20.2). The paper's exponent from Lemma 10.4 is
a − 7/8 + h(17/50 − 1/6) − a·l_y − ℓ/2 + qℓ + d(R + δ/2 − 17/50), with a = (1+δ)/2. This is exactly the formula in
floor/scales.py (with B in place of 7/8), so §50–§53 are evaluated on the paper's own endpoint. The paper's floor bound
(20.5), E(h) ≤ C0 + (3/4)δ0 = −7/1200, is the floor term of §50.

## 54. Local-to-global for LM: copies plus subtraction (PROVED reduction; HEURISTIC sizes)
By injectivity of (u, a) ↦ ua⁶ (paper Lemma 17.6), the all-elements moment at V = D^{1+c} splits by sixth-power-free scale:
  MS_all(V; D) = Σ_{U′ dyadic ≤ V} (V/U′)^{1/6}·MS_sf(U′; D)·(1 + o(1))   (coprimality averaged).
Copies keep only the term U′ = U and drop the others by positivity. That is lossless only if the scales near U′ ≈ D (and
the rows-dominant scales U′ ∈ (D, V]) carry no mass. They carry their diagonal, which is ≈ V^{1/6}D^{11/6}. Hence:
**If** (a) MS_all(V; D) and MS_sf(U′; D) for U′ ∈ (D, V] have asymptotics with power-saving error (rows-dominant side), and
(b) MS_sf(U′; D) ≥ (1 − D^{−η′})·diag for U′ ∈ [D^{1−c′}, D] (a lower bound just beyond the conductor), **then**
subtracting the top and near-conductor scales leaves
  MS_sf(U; D) ≤ c·UD + c·U^{1/6}D^{1+5(1−c′)/6} + U^{1/6}D^{11/6−η″},
i.e. slope α′ = (5/6)(1 − c′s/(s−1)) < 5/6 for every s. At s = 1.16 and c′ = 0.01 this is α′ = 0.773, a new record (MODEL; the frontier is α′ ≈ 0.814, §55).
Quantitatively, the exponent becomes e(t) = max(1, 1/6 + (5/6)(1 − c′)t). That is LM on [1, 1 + c′/(1 − c′)] and local slope
(5/6)(1 − c′) beyond, so the value is essentially that of LM on [1, 1 + c′] (§55) plus a small slope gain.
So a power-saving **asymptotic in a neighbourhood of the conductor** (two-sided in D/U) would improve the long slope
everywhere. By Theorem C both (a) and (b) are dual-LM statements: the rows-dominant asymptotic needs the dual
(Gauss-sum) mean square to equal its diagonal up to D^{2−η}. The paper's recursion supplies only ε-lossy upper bounds there.
The ζ-analogue of (b) is a power-saving variance of μ in intervals of length X^{small}, which is open (Matomäki–Radziwiłł
give o(H²)). This is a reformulation, not progress on the input: the near-conductor asymptotic is the smallest form of LM
that would move the record.
**The plateau constant is the coprimality factor (§45 addendum).** The columns are primary n coprime to 6, and ψ_u(n) = 0 when
(n, u) ≠ 1. So the true diagonal for row u is Σ_{(n,u)=1}μ²W², and its row average is the diagonal times
Π_{p ∤ 6}(1 − 1/(N(p)(N(p)+1))) = 0.9399 (product over prime ideals). The observed plateau 0.88–0.96 (mostly 0.92–0.95) matches
it. So the data say MS = (1 + o(1))·(true diagonal), with relative deviations of a few percent, at every tested s up to 76×
past the conductor. That is the LM asymptotic with the correct constant.
Source check for (a) (oai5.txt l.6991ff). The paper's canonical family estimate (Lemma 17.2) is
E = Z^{−V}Σ_f w(f)Σ_{q_k≪Z^M}|Z^{−N/2}Σ_n ᾱγ2νρ χ_n(k)χ_n(f)⁴ d(n)W|² ≪ Z^{N+V+ε} under 3M + 6z0 + c* ≤ 4F0. This is an ε-lossy
sharp-sieve-type upper bound with no main-term extraction, and its terminal step (reflected energy, Lemma 14.3) is also an
upper bound. So the paper's recursion does not supply the power-saving asymptotic (a); that would be new.

## 55. How much of the LM range is needed (MODEL, c61/lmrange.py; capc4 grid, current slot costs, baseline 0.87492)
LM (slope 0) only on s ∈ [1, s0], copies' slope 5/6 beyond:
| s0 | 1.00 | 1.05 | 1.10 | 1.16 | 1.20 | 1.30 | 1.50 |
|---|---|---|---|---|---|---|---|
| θ | 0.87492 | 0.87316 | 0.87140 | 0.86952 | 0.86864 | 0.86742 | 0.86667 |
LM just 5% past the conductor is already worth −0.0018 (dθ/ds0 ≈ −0.035 near s0 = 1). LM up to t* ≈ 1.16 gives 0.8695, about
65% of the way to 13/15. The table needs only the ε-lossy upper bound, so the smallest high-value target is
**Σ*_u |M_u(D)|² ≪ U^{1+ε}·D for D ∈ [U, U^{1+c}]** (LM just past the conductor), for some fixed c > 0. The stronger asymptotic
form, (1 + O(U^{−η}))·(coprime diagonal), would also improve the slope beyond U^{1+c} by §54. Even under GRH only the upper
bound is known: Gao–Zhao's quadratic asymptotic needs D ≤ U^{1−ε}. Partial DH from partial LM (same model; λ(σ) = R − (1−δ), worst amplitude class):
| LM on [1, s0] | DH (λ ≤ 0) holds for | max λ (at σ) |
|---|---|---|
| s0 = 1.00 (now) | none below 0.875 | 0.084 (0.716) |
| 1.05 | ≈ σ = 0.875 only | 0.068 (0.697) |
| 1.16 | σ > 0.781 | 0.038 (0.641) |
| 1.30 | σ > 0.659 | 0.013 (0.594) |
| 1.50 | all σ | 0 |
So LM up to t* ≈ 1.16 already gives a genuine density hypothesis for the sextic family on σ ∈ (0.781, 0.875].
By Theorem C the upper bound is an asymptotic-with-error
statement for the dual Gauss-sum family Σ_{h≤K}|B_h(D)|² slightly inside its
rows-dominant range (K = D²/U ∈ [D, D^{1+c/(1+c)}], just above the self-dual point K = D). The ζ-analogue is a power-saving variance of μ in intervals of length
X^{c/(1+c)}, which is open.

Partial slopes (same model, c61/lmrange.py): slope α′ on [1, s0], copies' 5/6 beyond.
| α′ \ s0 | 1.05 | 1.16 | 1.30 | 1.50 |
|---|---|---|---|---|
| 0 | 0.87316 | 0.86952 | 0.86742 | 0.86667 |
| 0.25 | 0.87370 | 0.87104 | 0.86844 | 0.86819 |
| 0.50 | 0.87418 | 0.87265 | 0.87080 | 0.87055 |
| 0.70 | 0.87467 | 0.87394 | 0.87372 | 0.87372 |
The model's baseline is 0.87492 (the record 0.8748389 adds refinements outside this model). Even the weakest entry
(α′ = 0.7 on [1, 1.05]) gains 0.00025 over the baseline, three times the 0.00008 gap to the record. If the gains are
roughly additive (MODEL, not checked in the record's own optimiser), slope 0.7 instead of 5/6 on a 5% window past the
conductor would already be a record.

**Record-moving frontier (MODEL, c61/lmfront.py, c61/lmeta.py).** The largest slope α′ on [1, s0] (5/6 beyond) that gains
≥ 0.00008, the model-baseline–record gap:
| s0 | 1.02 | 1.05 | 1.10 | 1.16 |
|---|---|---|---|---|
| α′ ≤ | 0.697 | 0.775 | 0.801 | 0.814 |
| (5/6 − α′)(s0 − 1) | 0.0027 | 0.0029 | 0.0032 | 0.0031 |
The last row is nearly constant. What counts is the integrated saving at the long-witness length t* ≈ 1.16. A uniform
power saving over copies, e(t) = max(1, 1 + (5/6)(t − 1) − η), gives θ = 0.87492, 0.87480, 0.87467, 0.87443 and 0.87402 for
η = 0.001, 0.003, 0.005, 0.01 and 0.02. Restricting the saving η = 0.003 to a window of lengths (c61/lmwin.py) gives:
- t ∈ [1.14, 1.18]: 0.87485, a gain of 0.00007, just short of the gap;
- t ∈ [1.10, 1.25]: 0.87480, the full gain.
So the smallest LM-type statement that would (in the model) move the record is
  Σ*_u |M_u(U^s)|² ≪ U^{1 + (5/6)(s−1) − 0.003} · U^s   for s ∈ [1.10, 1.25],
a saving of U^{0.003} over sixth-power copies on that window. By §54 and §60, no symmetry or subtraction argument gives
even this. It needs an asymptotic (main-term) input.

## 56. The dual coefficients are twisted CUBIC Gauss sums (NUMERICAL identity check, mmnum/hdcheck.py; corrects §29)
At all 2248 split primes of norm ≤ 2·10⁴:
- Hasse–Davenport holds exactly: g̃6(p) = χ̄_p(4)χ_p²(−1)·g̃3(p)²·g̃2(p), with 0 mismatches. Here g̃k is the normalised Gauss sum
  of χ_p^{6/k}.
- The quadratic Gauss sum g̃2(p) ∈ {1, i, −i} is determined by the primary generator π mod 4 (12 classes, none ambiguous).
- The cubic Gauss sum satisfies g̃3(p)³ = −π/|π| exactly (all phases ±π). The conjugate π̄ gives random phases.
- The cross terms of twisted multiplicativity match: (χ²)²·χ³ = χ⁷ = χ.
- Inert primes (23 primes q ≡ 2 mod 3, norm q² ≤ 4·10⁴): Hasse–Davenport defect 0, g̃2 = 1, and g̃3³ = −π/|π|.
- Composite n (mmnum/hdcomp.py): ε(n) = μ(n)g̃6(n)/((n/|n|)conj(g̃3(n))) is a 12th root of unity at all 17 976 split primes and
  all 25 557 products pq of norm ≤ 2·10⁵ (phases in multiples of π/6).
Hence, for squarefree primary n,
  μ(n)g̃6(n) = ε(n)·(n/|n|)·conj(g̃3(n)),
with ε a finite-order function (determined by residues mod 4·36 and quadratic-reciprocity signs). So the dual Möbius–Gauss
coefficients of Theorem C are **Grössencharacter-twisted cubic Gauss sums**, the objects of the cubic Kubota series of
K-type 1 (Heath-Brown–Patterson equidistribution). This is the "Möbius cancels after Poisson" of REPORT §2.
Correction to §29: a spectral route would run on the cubic metaplectic cover, where Patterson's theory is explicit. It does
not need the unknown sextic theta coefficients. The row characters χ̄_n(h) are still sextic in h, so the dual family is
mixed (sextic in h, cubic in the coefficients), and Theorem C's self-duality is unchanged.
The open input, stated in cubic terms, is a power-saving asymptotic for
  Σ_{h ≤ K} |Σ_{n~D} ε(n)(n/|n|)conj(g̃3(n))χ̄_n(h)W(N(n)/D)|²  at K = D²/U slightly above D,
i.e. an asymptotic sextic large sieve for weight-1 cubic Gauss sums just inside the rows-dominant range.
Literature remark (CITED, 16 Oct search). The closest general framework is the asymptotic large sieve of
Conrey–Iwaniec–Soundararajan (arXiv 1105.1176). It is built for the thick family of all primitive characters with q ≤ Q, and
it gets beyond-sieve asymptotics from shifted-convolution main terms of divisor-type coefficients. For Möbius coefficients
the corresponding input is a Chowla-type bound for Σ μ(n)μ(n+h) with power saving, which is open. Dunn–Radziwiłł
(2109.07463) correct the cubic large sieve by a main term, but only under GRH and with a tiny saving. No 2025–2026 result
gives an asymptotic sextic large sieve for Gauss-sum coefficients beyond the self-dual point.

## 57. Exact off-diagonal measurements and U = 10⁴ far beyond (NUMERICAL; mmnum/mmk3.c, lm_exact.py, lm_beyond.py)
mmk3 also returns each row's exact coprime diagonal td_u = Σ_{(n,u)=1}μ²W². So MS/Σtd is the exact relative size of the
off-diagonal, and with all rows there is no sampling. The ± is the row-fluctuation scale.
| U (all rows) | s = 0.9 | 1.0 | 1.05 | 1.1 | 1.2 | 1.5 | 1.8 | 2.0 |
|---|---|---|---|---|---|---|---|---|
| 10³ (3636) | 0.960 (s=0.8) | 0.986 | – | – | 0.977 | 0.963 | 1.005 | 0.950 |
| 3·10³ (10 872) | 1.006 (s=0.8) | 1.004 | – | – | 0.986 | 0.988 | 0.990 | 1.005 |
| 10⁴ (36 234) | 1.002 | 0.992 | 0.993 | 1.014 | | | | |
| 3·10⁴ (108 606) | 0.995 | 1.009 | 1.007 | 1.003 | | | | |
For s ≤ 1.1 the effective column count is only 1–5% of the rows, because W is narrow and rows count units, so these points
are rows-dominant in the effective sense (§44). They test side (a) of §54. The relative off-diagonal there shrinks
from ≈ 3% (U = 10³) to ≈ 1% (10⁴) to ≈ 0.6% (3·10⁴), roughly like U^{−1/2}: square-root-size fluctuation, no structured
main term. The random multiplicative control behaves the same.
The U = 3·10³ row is exact over all 10 872 rows and reaches cols_eff/rows = 1.05, 11.6 and 57 at s = 1.5, 1.8, 2.0. There the
off-diagonal is −1.2%, −1.0% and +0.5% of the exact diagonal, with row-fluctuation scale ±1.0%. The random multiplicative
control gives −1.3%, 0.0% and −0.9%. So 57 times past the effective conductor, the off-diagonal is at the
square-root-fluctuation level, with no sign of a term of the diagonal's order.
At and past the effective conductor for larger U (exact coprime diagonal; ± = row-fluctuation scale):
| U (rows) | s | cols_eff/rows | μ: MS/truediag | random mult. |
|---|---|---|---|---|
| 10⁴ (all 36 234) | 1.42 | 0.92 | 0.996 ± 0.006 | 0.998 |
| 10⁴ (all 36 234) | 1.65 | 7.6 | 1.012 ± 0.006 | 1.006 |
| 3·10⁴ (20 000 sampled) | 1.38 | 0.96 | 0.992 ± 0.008 | 1.000 |
| 3·10⁴ (all 108 606) | 1.38 | 0.96 | 1.0013 ± 0.0034 | 0.9991 |
| 3·10⁴ (30 000 sampled) | 1.50 | 3.3 | 0.9997 ± 0.0064 | 0.9988 |
| 3·10⁴ (15 000 sampled) | 1.62 | 11.4 | 1.0016 ± 0.0091 | 1.0010 |
| 10⁴ (2000 sampled) | 1.95 | 121 | 0.964 ± 0.025 | 1.003 |
| 3·10⁴ (20 000 sampled) | 1.50 | 3.3 | 0.994 ± 0.008 | 1.002 |
| 10⁵ (30 000 sampled) | 1.33 | 0.86 | 0.986 ± 0.006 | 1.000 |
| 10⁵ (30 000 sampled) | 1.45 | 3.4 | 1.004 ± 0.007 | 0.992 |
At the effective conductor (cols_eff/rows ≈ 1, the regime of §55's "smallest target"), the μ off-diagonal is −0.4% (U = 10⁴),
−0.8% (3·10⁴, sampled), +0.13% ± 0.34% (3·10⁴, all 108 606 rows) and −1.4% ± 0.6% (10⁵, 30 000 rows). At U = 10⁵,
3.4 times past the effective conductor, it is +0.4% ± 0.7%. No trend in U: the deviations scatter around 0 at
the 1–2σ level, as the random multiplicative control does. That is within one to two fluctuation units, and the same as the random multiplicative control.
U = 10⁴ far beyond the effective conductor (3000 sampled rows, normalised by the plain diagonal; plateau = 0.94):
| s | 1.5 | 1.65 | 1.8 |
|---|---|---|---|
| cols_eff/rows | 1.92 | 7.62 | 30.4 |
| μ | 0.933±.019 | 0.947±.019 | 0.967±.021 |
| random multiplicative | 0.931 | 0.940 | 0.964 |
| squarefree indicator | 0.658 | 0.506 | 0.260 |
At all three U tested so far beyond the effective conductor (10³, 3·10³, 10⁴), μ sits on the coprimality plateau, and the
smooth-type control collapses.
WMDS form of the dual (PROVED algebra, heuristic programme). For order 6, χ = χ^{−2}·χ³, so χ̄_n(h) = (h/n)_3(h/n)_2. Since
conj(g̃3(n))(h/n)_3 = conj(g̃3(h, n)) (the cubic Gauss sum with additive twist h, for (h, n) = 1),
  B_h(D) = Σ_n ε(n)(n/|n|)·conj(g̃3(h, n))·(h/n)_2·W(N(n)/D),
which is a partial sum of the h-th Fourier coefficient of a cubic Kubota–Eisenstein series of K-type 1, twisted by the
quadratic symbol (h/·)_2. The dual mean square Σ_h w(h/K)|B_h|² is therefore a double Dirichlet series in (h, n) of
Weyl-group-multiple-Dirichlet-series type (cubic Gauss sums with quadratic twists, cf. Brubaker–Bump–Friedberg A2/A3).
The quartic route's lemma (L³) was attacked with exactly this technology (Chinta–Gunnells G2 WMDS, tasks #26–#27) and
stopped at the convexity point Q^{3/2}. The sextic LM needs the analogous series beyond its convexity range with a
power-saving error. That is a concrete research programme, not a step available in this session.
Relation to the large-sieve extremisers (PROVED reading + HEURISTIC). Dunn–Radziwiłł's sequences attaining the cubic
(QN)^{2/3} term are β_b·g̃3(b) with β ≥ 0 (§11). Their excess is the residue of the K-type-0 cubic Kubota series (the cubic
theta function). LM's dual sequence b_n = ε(n)(n/|n|)conj(g̃3(n)) is the same object twisted by the K-type-1 Grössencharacter
n/|n|. That series has no theta pole, so the extremal main term is absent. This is why the numerics see no excess (§41).
It is also why only the true asymptotic, with error below K·D by the factor U/D (Theorem C), is at stake, not a
large-sieve constant. Unconditional bounds for twisted Kubota partial sums do exist (convexity of the K-type-1 series), but
they enter only through the (KD)^{2/3} term, which copies already beat by U^{1/6}D^{1/6}.

## 58. Function-field analogue of LM (CITED large-q limit; NUMERICAL small degree; ff/ff_lm.py, ff/brute.py)
Over F_q[t] with q ≡ 1 (mod 6), take rows u monic and sixth-power-free of degree d (U = q^d), and put
M_u(k) = Σ_{deg n = k} μ(n)(n/u)_6 (length D = q^k, s = k/d). Then M_u(k) is the T^k coefficient of 1/L(T, χ_u), and L is a
polynomial of degree ≤ d − 1. So M_u(k) = q^{k/2}h_k(Θ_u), where h_k is the complete symmetric polynomial in the unitarized
Frobenius eigenvalues. The script computes L from χ_u(P), deg P < d (Newton's identities). It is checked row by row against
brute-force Möbius sums.
- **Large q (CITED + representation theory).** If the χ-eigenspace monodromy of the family contains SU(N), Deligne–Katz
  equidistribution gives avg_u |M_u(k)|²/q^k → ∫|h_k|² = 1 for every fixed k, because Sym^k is irreducible. So LM holds in
  the q → ∞ limit at every length, with constant exactly 1. This is the same mechanism as Keating–Rudnick's variance of μ
  in short intervals and progressions (Algebra & Number Theory 10 (2016)). The coprimality plateau 0.94 of §45 becomes
  1 − O(1/q).
- **Small degree (NUMERICAL), and why it is no proxy.** For d ≤ 5 (q = 7, 13) every generic row has |h_k| ∈ {0, 1}: the
  eigenvalues are a common phase times roots of unity, i.e. the eigenspace has finite (CM) monodromy, so the large-q
  limit is not Haar there. At d = 6 (q = 7, all 117 642 rows), rows with distinct L-roots give MS/truediag = 0.85, 0.50, 1.07,
  0.66 and 0.85 at s = 1, 1.5, 2, 2.5, 3, with max|M|²/q^k ≤ 3.7. The full average is driven by CM rows with repeated
  Frobenius eigenvalues (e.g. u = t³ + c at d = 3, u = t⁴ + c at q = 13, d = 4), where |h_k|² = binom(k+m−1, m−1)² grows
  polynomially in k. These give MS/truediag up to 17.7 at s = 3. In Z[ω] the analogue would be rows whose L(s, ψ_u) has
  multiple zeros or tiny L′(ρ), which is the explicit-formula form of LM (LM_PROBLEM "Analogues"). So at accessible degree
  the function-field family measures finite-monodromy effects, not the number-field question.
- **Fixed q, d → ∞.** RH (Weil) gives only |h_k| ≤ binom(N+k−1, k), a loss U^{O(1/log q)}. Negative-moment bounds of
  Bui–Florea–Keating type (Trans. AMS 376 (2023), hyperelliptic ensemble) are the function-field tool that would give LM on
  average, and they use RH. Over Q(ω) the corresponding input is a family negative second moment
  Σ*_u ∫|W̃(σ+it)||L(σ+it, ψ_u)|^{−2}dt ≪ U^{1+ε} at σ = 1/2 + ε, plus control of the residues D^ρ/L′(ρ) at zeros with β > σ.
  Together these give LM with ε-loss: shift the Mellin contour of M_u(D) to Re z = σ. Under GRH both hold. Unconditionally
  neither is available.

## 59. Which σ-window of DH moves θ (MODEL, c61/dhrange.py; capc4 grid, current slot costs, slope 5/6, baseline 0.87492)
The DH count R ≤ 1 − δ is granted only for σ = (1+δ)/2 in a window, and the current counts are kept elsewhere:
| window | θ | | window | θ |
|---|---|---|---|---|
| [0.50, 0.875] | 0.86667 | | [0.50, 0.70] | 0.87492 |
| [0.60, 0.875] | 0.86922 | | [0.50, 0.75] | 0.87446 |
| [0.65, 0.875] | 0.87109 | | [0.50, 0.80] | 0.87255 |
| [0.70, 0.875] | 0.87335 | | [0.68, 0.78] | 0.87449 |
| [0.72, 0.875] | 0.87418 | | [0.70, 0.76] | 0.87482 |
| [0.75, 0.875] | 0.87492 | | [0.72, 0.75] | 0.87488 |
- DH above σ = 0.75 alone is worthless: the binding class sits at σ ≈ 0.733. Short windows around the binding point are also
  nearly worthless, because once that class improves, the binding moves to a neighbour on both sides.
- Each lower endpoint is worth about −0.005 per 0.1 in σ, and the upper end must reach about 0.8 or more.
- Partial LM is worth more than the DH window it implies. LM on [1, 1.16] gives DH only on σ > 0.781 (§55), which alone is
  worthless here, yet θ = 0.86952. The reason is that partial LM lowers λ at every σ (max λ 0.084 → 0.038), not just where
  it reaches 0.
So the DH formulation is the endpoint, but the right intermediate target is the LM slope, not a DH window.

## 60. Twisted copies give nothing (PROVED, elementary)
Natural attempt: beat the copy density 1/6 by also using w = u·a³ or w = u·a², whose characters
ψ_{ua³} = ψ_u·(a/·)_2 and ψ_{ua²} = ψ_u·(a/·)_3 have conductor only N(u)N(a). The all-elements moment at norm W then bounds
Σ_{u,a}|M_{ψ_u ξ_a}(D)|² for the twisted sums. Formally, if those were copies of ψ_u, the count would give slope 3/4 (cubes) or
2/3 (squares). They are not copies. ψ_u ξ_a = ψ_u on all n coprime to ua only when ξ_a is trivial, i.e. a = b² (resp. b³),
which gives back w = ub⁶. No inequality runs from Σ_a|M_{ψ_u ξ_a}|² to |M_u|²: the a-sum is a quadratic (cubic) large-sieve
average of the fixed sequence μψ_uW, of size ≍ A·D whatever M_u is. So sixth-power copies are the only exact symmetry of
the family. The slope 5/6 = 1 − 1/6 is the copy density, as the quadratic analogue's 1/2 = 1 − 1/2 is (§29), and Theorem C
says Poisson cannot add to it. Any improvement must come from an asymptotic, not from more symmetry.

## 61. Every Part II margin re-checked at the 13/15 parameters; DH is needed only for d ∈ [2/5, h+ζ] (PROVED relative to the paper; exact rational check, margins/part2_margins.py)
Paper §20 lists every margin of the seven-eighths proof: (20.3)–(20.7), (20.11), (20.12), §20.3, §20.5 and Prop. 20.3. Each was
re-evaluated in exact rational arithmetic at the Conditional-Theorem-D geometry: ℓ = ℓ_sup(λ) − λ/100 with
ℓ_sup = (1−2λ)/(5+6λ), x = (1−ℓ)/2, h = 1 − x + ℓ, l_y = 1 − ℓ − x, σ0 = B = 11/12 − ℓ/4, floor a0 = 1/2 + λ/8, ζ = λ/10⁴. The
requirement of (20.12) is E(d) < 0 for δ ≤ 2σ0 − 1, and E(d) < (1+δ)/2 − σ0 for larger δ ≤ 5/6, since β* ≤ 11/12 by Part I.
Worst value of each margin (λ = 0.001 with ℓ = ℓ_sup − λ/100, so θ = 0.866829 against the limit 0.866826; λ = 0.01, 0.02, 0.04 behave the same way; values from the run with
x = (1−ℓ)/2, which differ from the final interior-x run in the fifth decimal):
| item (paper ref.) | value | needs |
|---|---|---|
| E(d), DH_λ count, d ∈ [2/5, h+ζ] (20.4) | −1.2·10⁻⁵ | < 0 |
| E(d), paper's own no-prime count 76/75 − 2δ/3, d ∈ [1/100, 2/5] (20.11) | −0.0147 | < 0 |
| floor bin, R = 1 (20.5) | −0.0006 | < 0 |
| small rows h(17/50 − 1/6) − l_y/2 + 2d_min (20.6) | −0.0417 | < 0 |
| principal margins m_w = l_y/20, m_z = h/600 (20.3) | 0.0200, 0.0013 | > 0 |
| prime supply ℓ/(h+ζ) − 7/37 (§20.5) | 0.060 | > 0 |
| low side: completed-row excess, Gram b = 1 − ℓ − 2x (REPORT §3) | −0.0004, 0 | ≤ 0, ≥ 0 |
Consequences:
- **The hypothesis shrinks to d ∈ [2/5, h+ζ].** Below 2/5 the paper's own no-prime count (§20.5, t = 1) already works: at
  d = 2/5, E = −0.2587 + δ/3 ≤ −0.0142 for δ ≤ 11/15, and it stays below (1+δ)/2 − σ0 up to δ = 5/6. So neither the detected
  count nor any selected-prime machinery is needed below the top scales. This removes the mesh adjustment that §51's
  "existing counts below 1/2" variant would have required. The paper restricts selected primes to d ≥ 1/2, through
  w_i = ℓ_i/d ≤ 2ℓ/K.
- **The remaining changes are all region margins, each checked against the exact local Euler data:**
  (a) the D1 edge Re s ≥ 51/100 lowered to 1/2 + ε, for the floor (§52(ii), §53);
  (b) the D2 edge Re s ≥ 7/8 lowered to 0.86 for the principal row, together with Lemma 10.5's hypothesis σ0 ≥ 7/8 (§52(iii));
  (c) for the small rows, D1(3/8) replaced by D1(1/3). This is only a choice of ε0: Definition 10.1 covers every D1(ε0). The line (β*+e, 1/2, 17/50) with β* > 13/15 has Re(s+w) > 1.366.
      There the selected-prime exponents (−0.867, −2.04, −2.37, −1.54) are all ≤ 0, and the unselected normal-convergence
      exponents (−1.367, −2.907, −2.54, −3.24, −2.407) are all < −1, so (16.15) holds verbatim.
- C(s) = s − 5/6 + l_x/3 + ℓ/6 becomes s − 2/3 at the limit, so C(σ0) = 1/5 > 0. The prime normaliser A_T, κ_P = 7ℓ/(16K) and
  the large rows (z∞ → ∞) do not depend on θ.
So, relative to the paper and to these three table-level region checks, **Conditional Theorem D′ holds with the hypothesis
DH_λ on bins a ∈ (1/2 + λ/8, β*] at row scales d ∈ [2/5, h+ζ]**. Every other ingredient is the paper's own. REPORT §3's open
remark on small rows ("should be re-verified for each concrete (ℓ, x)") is closed: −0.0417 at the limit.
**Addendum (low side, Prop. 15.3).** The low bound is Z^{l_x/2 + b/12 + ε} with b = 1 − ℓ − 2x = l_y − l_x, and
C(s) − (l_x/2 + b/12) = s − (11/12 − ℓ/4) identically. This is exactly why B = 11/12 − ℓ/4, and the low margin is β* − B. The proof's
length inequalities, taken at the worst tuple d = ℓ, are X′ ≥ 1, Y′^{11/6}/P_a ≫ Z^{1/12} and M′ = 1 − 3ℓ > 0. At the limit
x = (1 − ℓ)/2 one has b = 0 and P_a = Z^{l_y − l_x} ≍ 1. To keep P_a ≥ 1 with room, take x at the midpoint of its allowed interval.
Then b > 0, and all lengths hold with margins ≥ 0.11 (part2_margins.py, λ = 0.001 … 0.04, final run: ALL OK).
Independent float cross-check (floor/scales25.py, optimiser of §51 with this split): DH on d ∈ [2/5, h] and the paper's count on
d ≤ 2/5 give θ = 0.869792 (floor 0.51), 0.866986 (floor 0.501) and 0.866699 (floor 0.5001). These are identical to DH at every
scale, and tend to 13/15.
Validation: the same exponent function reproduces the paper's own (20.3), (20.5), (20.6) and (20.11) exactly at the paper's geometry
(ℓ = 1/6, x = 17/48, σ0 = 7/8): −529/2400 + 25δ/96, −7/1200, −79/800, 23/960, 13/9600.

## 62. LM is stronger than DH: what the architecture needs from the zeros (HEURISTIC, explicit formula with |L′(ρ)| ≍ 1)
For a row with a simple zero ρ = β + iγ (bounded height) and no zeros to its right, M_u(D) ≈ Σ_ρ W̃(ρ)D^ρ/L′(ρ). Suppose
N(δ′) ≍ U^{R(δ′)} rows have a zero at δ′ = 2β − 1 and |L′(ρ)| ≍ 1. Then
  Σ*_u |M_u(U^s)|² ≈ U^s · max_{δ′} U^{R(δ′) + sδ′},
so LM_α on (1, 3/2] holds iff R(δ′) ≤ 1 − δ′ − (s − 1)(δ′ − α) for all δ′ and s. Consequences:
- **In a DH-saturated world** (R(δ′) = 1 − δ′ at every δ′ ≤ κ), the slope is α = κ ≈ 0.75. The closed form (§49) then gives DH
  only at δ ≥ 0.75. So the architecture cannot prove DH from LM in that world. Even if DH were true with equality, LM
  would fail.
- **What LM_α at slope α asks of the zeros:** strictly sub-DH counts above α, with R(δ′) ≤ 1 − δ′ − (1/2)(δ′ − α) at s = 3/2.
  Alternatively, |L′(ρ)| must be large on average at such zeros, or their residues must cancel. For example, the slope
  α = 0.466 needed at the binding bin allows at most U^{0.108} rows with zeros at σ = 0.875, against DH's U^{0.25}. At s = 1 the condition is
  exactly DH. LM at s = 1 is proved (the critical-length Möbius moment) while DH is open, so the |L′| ≍ 1 heuristic is not a
  theorem. The adversary's dip at r = 1 (§40, §47) is exactly the escape: the zeros' residues cancel at the conductor.
- **Single rows.** One row with a zero at β gives D^{2β}, which stays below U·D iff β ≤ (1 + 1/s)/2. With the proved
  bound β* ≤ 0.8748389 this holds for every s ≤ 1.334. So exceptional single rows cannot obstruct LM on the record-moving
  window s ∈ [1.10, 1.25] of §55. The ζ-analogue (power-saving variance of μ in short intervals) is different: there one
  zero near Re s = 1 already blocks it, and no zero-free strip is known. The family LM is not blocked that way, because
  the architecture itself provides the strip.
So LM is a top-down, sub-DH statement about the family's zeros, weighted by 1/|L′(ρ)|². It is sufficient for 13/15 but
morally stronger than the DH it is used to prove. A proof of DH that avoids LM would need a different long witness.

## 63. LM_α ⟹ DH on δ ≥ α directly, with no selected primes (PROVED relative to the paper; its endpoint step (20.7))
Take the detector (Prop. 8.3) with t = 3/2, as the paper does at δ = α. The Möbius witness then has r ∈ [1 − O(ε), 3/2], and
the no-slot amplified estimate (19.3) with LM slope α counts the bin rows by U^{1 − α + (α − δ)r + ε}, or U^{1 − δr} for r < 1.
For δ ≥ α the maximum is at r = 1 and equals U^{1−δ}. Bins with δ < α satisfy U^{1+ε} ≤ U^{1−δ+α} trivially. So
LM_α ⟹ DH_α, and with D′ (§61, DH_THEOREM.md §3a) LM_α ⟹ θ ≤ 11/12 − (1−2α)/(4(5+6α)) → 13/15 as α → 0.
- This needs neither the selected-prime moments nor Lemma 18.1, which is stated only for κ ∈ [3/4, 1] while κ may now drop
  to 11/15.
- For fixed α > 0 it is much weaker than the detected count (§49), because λ = α uniformly. It is the clean slot-free
  bridge in the limit α → 0.
- The paper's own endpoint (20.7) is the case δ = α = 5/6.
- Range note: the Cycle61 paper (c61paper.tex (1.3), (4.8)) imports β* ≤ 7/8 and works with κ = 2β* − 1 ∈ [7/10, 3/4]. So the
  sub-7/8 records, ours included, are not exposed to the [3/4, 1] statement range of Lemma 18.1. Near 13/15, κ ≥ 11/15 > 7/10
  stays inside that extended range.

## 64. Value of a partial slot-cost improvement, and three more closed doors (MODEL c61/cIscan.py; PROVED/HEURISTIC remarks)
Model θ (capc4 grid, record geometry) as a function of the inverse slot cost c_I, all else fixed:
| c_I | 2 | 1.95 | 1.9 | 1.8 | 1.6 | 1.4 | ≤ 1.2 |
|---|---|---|---|---|---|---|---|
| θ | 0.87492 | 0.87478 | 0.87471 | 0.87444 | 0.87373 | 0.87296 | 0.87289 |
So any power-saving relaxation of r + 2z ≤ 1 in the marked Möbius moment (to r + 1.95z ≤ 1) would pass the record in the
model. c_I ≤ 1.4 is worth about 0.002, twelve times all progress since 7/8. The input is below the conductor (NZ ≤ U), yet
still blocked:
- **Hölder with high slot moments (PROVED arithmetic).** Σ|MP|² ≤ (Σ|M|^{2a})^{1/a}(Σ|P|^{2b})^{1/b}, with the 4th moment of M from
  copies (loss 0.408 at r* = 0.745) and P^b sharp for bz ≤ 1/2. This gives exponent 1 + r + z + 0.408/b ≥ 1 + r + z + 0.816z.
  The loss 0.816 per unit slot length exceeds the slot gain 2q = 0.466 at the binding bin.
- **Marks after Poisson (PROVED structural).** Summing the twisted sums T(χ_{p1}χ̄_{p2}) over pairs turns the dual into the
  marked moment Σ_h |B_h(N)Q_h(Z)|² of the dual Gauss-sum family, with Q_h = Σ_p a(p)g̃6(p)χ̄_p(h). This is the same problem of
  length NZ in the dual, so it is self-dual again.
- **Copies in the dual, and symmetries forcing T(χ) = 0.** Copies fail because the twist χ_{p1}χ̄_{p2} is sextic, so χ(k⁶) = 1
  and the copy sums do not cancel. As for symmetries, none of u ↦ εu (changes |M_u| by a small-conductor twist), u ↦ ua⁶
  (χ-invariant) or u ↦ ū (relates different pairs) annihilates T(χ).
- **Twisted rows as copies (HEURISTIC, extends §60).** For a group H of characters mod q, |M_u|² ≤ Σ_{ξ∈H}|M_{u,ξ}|², and for
  quadratic H the twisted rows are the family members ud³ (d | q). But a zero of L(s, ψ_u) does not pass to ψ_uξ. So twists
  do not amplify large values; only exact copies do, which is why the copy efficiency 1/6 governs both α′ and c_P.

## 65. Session 22:28Z: anatomy of the slot cost and of the dip (PROVED readings + HEURISTIC structure)
**Where c_I = 2 comes from, precisely (paper Cor. 14.1, Lemma 14.2, Lemma 14.3).**
- Each active mark q (slot prime in the active branch of 1_{q|A} = 1 − χ_q(A)⁰) enters the reflection conductor
  c = c_F(h0)·r, r = Π(active primes).
- Cubic-theta reflection maps a sum of length X to one of length ≍ q_c²/X, so the mark enlarges the dual by q², i.e. +2z_a.
- Lemma 14.2 then pairs the dual variable n with the active-mark product P through χ_P(n)^{−2}. Its bracket
  {N + L + (NL)^{2/3}} saves at most L = q_P, i.e. z_a, by orthogonality of the L characters. Net: +z_a beyond cost 1.
- Re-deriving the marked dual through twisted multiplicativity, g̃6(mp) = g̃6(m)g̃6(p)(p/m)_3, gives a frequency-shifted sum
  B^{MP}_k = Σ_p a′(p)χ̄_p(k)B_{kp⁴}(N). Reflecting each shifted dual separately again puts p into the conductor, which
  enters squared, so this reproduces the paper's +2z_a. The cost of 2 is intrinsic to theta reflection with marks inside
  the reflected index, as §32 said.

**Dip structure at leading order.**
- For the fixed dyadic weight W (log-support log 4), W̃^{(k)}(ρ0)/k! decays like (log 2)^k/k!. So for a cluster of radius
  ≲ 0.2, the residue sum is W̃(ρ0)·m0(D) + O(radius·log 2) relative, with m0(D) = Σ_j D^{ρ_j}/L′(ρ_j).
- The adversary's dip at D = U is therefore a cancellation of the row's own "Mertens sum" m0. It does not depend on the
  choice among fixed-support weights, so using several dyadic weights cannot kill it.
- In the Hermite–Genocchi positivity of §47, the conductor term (1/2)log q in L′/L contributes a linear phase in height.
  That is a pure ξ-shift (functional-equation shift) and does not break positivity.
- What breaks positivity is the curvature from zeros just outside the cluster, which is §47's obstruction (ii).

**Band-limited reformulation (HEURISTIC).** In r = log D/log U, a cluster spread over heights X makes f an exponential sum of
type ≈ X log U.
- A tent dip of depth 0.03 log U over (r*, t*), of length 0.417, needs zero density of f of order 0.015 log U on the
  window. Pólya's bound X log U/π per unit r then forces X ≳ 0.11, consistent with §47's 0.2–0.5.
- So the right dichotomy is in height spread, not in zero count. Clusters narrower than ≈ 0.1 in height cannot dip;
  wide ones need a pair-type count.

## 66. Three record-level levers checked (MODEL, c61/thrad.py; exact record pipeline with B = B0 + 7/10000)
- **Binding constraint.** The record-style optimiser gives θ = 0.87483906 (h = 0.79578). The margins at the optimum are
  +7.9·10⁻⁹ (d = h, new count), +2.4·10⁻⁴ (d = h − 10⁻³, old count) and +0.069 (d = 1/2). So only the top-scale new-count
  exponent binds, through R(0.467, 1/2). Widening the old/new split or the d = 1/2 range cannot help.
- **Theorem 2.2 child density (§6.6 of SESSION3).** Replacing the Part I density ρ(δ) of the high-bin children (δ ≥ q0 = 0.68)
  by the record's own count R(δ, p) gains nothing. At δ = 0.68…0.74, max_p R(δ, p) = ρ(δ) to 10⁻⁴ for every κ ∈ [0.75, 5/6].
  The p = 0 class has no slot help and is exactly ρ.
- **Off-line numerator moments.** Xu–Sang's amplified first moment uses the fourth moment at Re w = 1/2, which gives
  (3R+1)/4. At Re w = 1/2 + η, copies give Σ|L|^{2k} ≪ U^{1+ε} for k ≤ 2/(1 − 2.4η). The Hölder gain in the high exponent is
  (1 − R)·d·0.6 ≈ 0.194 per unit η, while the Y-power loses l_y ≈ 0.46 per unit η. To the left of the line, the
  functional-equation factor costs d − l_y > 0. So Re w = 1/2 is optimal and higher moments do not pay.

## 67. Sensitivity of the record to the affine constants (MODEL, c61/thsens.py; exact record pipeline)
| change | θ |
|---|---|
| record (B = B0 + 7·10⁻⁴, Λ = 4.0343) | 0.87483906 |
| B + 10⁻³ / + 2·10⁻³ | 0.87483118 / 0.87482330 |
| Λ − 0.05 / − 0.2 / − 0.5 | 0.87482414 / 0.87481583 / 0.87480172 |
- dθ/dB ≈ −7.9·10⁻³. The Λ gain saturates near 4·10⁻⁵ because the Theorem 2.2 window then binds.
- The whole affine-constant lever is worth at most about 4·10⁻⁵. That is about thirty times the last record step, but not a
  drastic change.
- Both remaining refinements are already used or closed. Using amplitude classes in Theorem 2.2's high-bin children is
  timbalud's marked high-zero counting (B0 + 7·10⁻⁴). Using the record's R for the children gives nothing (§66).
- Plain sixth moment via copies at the binding m* = 0.417 gives count 1/6 + 3m(5/6 − δ) = 0.626, above A_P = 0.5935. So the
  integer-moment route on the plain leg is also closed, consistent with Theorem A.

## 68. Session 23:00Z: the z-contour, the threshold q0 of Theorem 2.2, and where c_I = 2 sits in Lemma 14.2
**The Re z = 17/50 contour is not a lever (PROVED, paper (20.4)).** In E(d) = a − 7/8 + h(z_r − 1/6) − a·l_y − (1−a)ℓ − (δ/2 − q)ℓ
+ d(R + δ/2 − z_r), the coefficient of z_r is h − d. It vanishes at the binding top scale d = h and is positive below it, so
the smallest admissible value z_r = 17/50 (the lower edge of region D1, Lemma 7.1 (7.14)) is already optimal.

**Theorem 2.2's constants are a one-parameter family (PROVED reading of Cycle61 (2.10)–(2.23)).** The low/high bin split q0 is a
free parameter. Given q0:
- K_a = (2q0/9)/(1 − 4q0/3), the low reserve (2.17) with zero g-coefficient;
- ρ0 = ρ(q0), the Lemma 2.1 density;
- c = K_a(5 − 6ρ0);
- B0(q0) = (1 − ρ0 + c)/(q0 + c), the high reserve (2.23) at zero;
- Λ(q0) = (q0 + 9c/2)/(q0 + c), the identity (2.22).

At q0 = 17/25 this reproduces B0 = 14200/14361 and Λ = 19312/4787 exactly. The side conditions of (2.20) (4ρ − 1 − δ ≤ 0 etc.)
hold for q0 ∈ [0.65, 0.74]; Lemma 2.1 needs a ≥ 0.825, i.e. q0 ≥ 0.65. The capacity (B − 2m)/Λ at the binding m = 0.417 is
0.03840 / 0.03844 / 0.03842 / 0.03837 / 0.03813 at q0 = 0.65 / 0.66 / 0.67 / 0.68 / 0.70.

Record pipeline (c61/thq0.py, timbalud's +7·10⁻⁴ kept fixed):

| q0 | 0.66 | 0.67 | 0.68 (record) | 0.70 |
|---|---|---|---|---|
| θ | 0.87483792 | 0.87483779 | 0.87483907 | 0.87484557 |

So q0 ≈ 0.67 gains about 1.3·10⁻⁶ (MODEL). Timbalud's high-bin saving enters B as s/(q0 + c), so it grows as q0 decreases;
c61/thq2.py rescales it and allows separate q0 for the detector and numerator uses. A real version needs Cycle61 Theorem 2.2
and timbalud's marked high ledger (constants 519/4550, 224/1755, 13q/400, 1121M/1250 are tuned at q0 = 17/25) re-proved
at the new q0. This is a squeeze of order 10⁻⁶, not a structural gain.

**Where c_I = 2 sits in Lemma 14.2 (PROVED reading, paper (14.7)).** In the binding regime (K ≥ NB, N ≥ L²) both large sieves in
(14.7) are at their diagonals:
- the quadratic sieve gives Σ_k|Σ c_n χ_k(n)|² ≈ K‖c‖² for any c;
- the cubic sieve gives the mean square of A(n) = Σ_P a_P q_P^{−1/2} χ_P(n)^{−2} over n as N·Σ|a_P|²/L.

So (14.7) is sharp for the quantity it bounds, and no structure of the dual coefficients β(n, b) can be used after it. A
pointwise Heath-Brown–Patterson bound for twisted theta sums (Σ_n β(n)χ(n) with a small-modulus χ) would have to enter
before the quadratic sieve over k. But β is twisted by χ_k³χ_P^{−2}, whose modulus kP is large, so no pointwise saving is
available. In addition, the twisted theta series has the Patterson pole at 5/6, a bias of size N^{5/6}, which dominates
the square-root scale once L > N^{1/3}. The loss is therefore upstream, in Corollary 14.1: the mark q | A enters through
1 − χ_q(A)⁰, and twisted multiplicativity makes the q-part a cubic character of conductor q in the complementary index.
Its reflected length is then multiplied by q². This confirms §32 from the other side.

**Why the sixth moments are out of reach of multiple Dirichlet series (CITED framework + PROVED classification).** The moment
Σ_u Π_{i≤r} L(s_i, ψ_u^{±1}) u^{−w} has the star graph K_{1,r} as its Weyl-group diagram (Bump–Friedberg–Hoffstein, Diaconu–Goldfeld–Hoffstein):
- r ≤ 3 is finite type (A_2, A_3, D_4);
- r = 4 is affine D̃_4, the critical fourth moment, which the paper reaches by theta reflection;
- r ≥ 5 is indefinite.

The plain sixth moment (§6, §34) and the numerator sixth moment (QRH catalogue P30-14, which moves the record once its loss is
< (1 − R)/2 ≈ 0.20) are r = 6, with no meromorphic continuation theory. The Möbius inputs (LM, MLV, Möbius fourth moment) have
no Dirichlet-series structure at all, since 1/L has no multiple Dirichlet series. Quadratic analogue: the sixth moment of
quadratic Dirichlet L-functions is open; the cubic moment (D_4) and the fourth moment (D̃_4) are known.

**Addendum to §68: timbalud's ledger ties its saving to q0 (MODEL, c61/thBL.py; Lean source read).** Timbalud's B0 + 7·10⁻⁴ comes from
the marked high ledger:
- the reserve needed is 73/20000;
- the saving is 37/10000 ≤ κ·removed + δ·frozen + (224/1755)·min(live, 13q/400), with supply ≥ 33M/1000;
- κ = 519/4550 = q0 − 3ρ0/2 comes from CoupledRadicalSaving.first_common_uniform, which uses the bound −(δ − 3ρ/2)c.

The cheapest channel is "removed", so the saving is s(q0) = min(37/10000, 0.033·(q0 − 3ρ(q0)/2)). It falls below 37/10000 as soon as
q0 < 0.6785. With ΔB = (s − 74·10⁻⁶)/(q0 + c), the capacity at m = 0.417 peaks at q0 ≈ 0.669 (0.038580 against 0.038542), and the pipeline gives:

| q0 | B | Λ | θ |
|---|---|---|---|
| 0.68 (record) | 0.989479 | 4.0343 | 0.87483908 |
| 0.67 | 0.986977 | 3.9652 | 0.87483840 |
| 0.668 | 0.986441 | 3.9513 | 0.87483851 |

The q0 lever is therefore worth about 7·10⁻⁷. Realising it would need the Theorem 2.2 cells, the Cycle63 ledger constants and the
record certificate re-proved, so it is not pursued. The 3·10⁻⁶ figure above assumed a q0-independent saving and is superseded.

## 69. Möbius-weighted cubic Gauss sums carry no Patterson bias (NUMERICAL, mugauss/mug2.py; CITED Dunn–Radziwiłł arXiv 2109.07463)
**Setup.** Over Z[ω], with c ≡ 1 mod 3 squarefree, take g(c) = Σ_{d mod c}(d/c)₃ e(Tr(d/c)) and g̃ = g/|c| (Dunn–Radziwiłł's
convention). The code checks g̃(ϖ)³ = −ϖ/|ϖ| at every prime to 10⁻¹³. Composite c are built by twisted multiplicativity
g(c·ϖ) = (c/ϖ)₃² g(c)g(ϖ), using cubic reciprocity; inert primes are computed directly over F_{q²}.

**Results to X = 2·10⁵** (partial sums at X/2^k, 70 554 squarefree c; mugauss/mug_2e5.log):
- Σ_{N(c)≤X} g̃(c) / X^{5/6} = 0.355, 0.356, 0.358, 0.350, 0.361, … at X = 2·10⁵, 10⁵, … This is the theta (Kubota–Patterson) bias.
- Σ_{N(ϖ)≤X} g̃(ϖ) / (X^{5/6}/log X) = 1.10, 1.15, 1.19, 1.15, 1.03 (18 020 primes). Dunn–Radziwiłł's constant under GRH is 2(2π)^{2/3}/(5Γ(2/3)) ≈ 1.006, so this is
  consistent up to lower-order terms.
- Σ_{N(c)≤X} μ(c) g̃(c) = +222.1, −19.7, −141.9, +1.0, +59.2, … at X = 2·10⁵, 10⁵, …. That is |S|/X^{1/2} ∈ [0.006, 0.63] on
  every scale, against 0.355·X^{5/6} ≈ 9300 for a theta-size bias at 2·10⁵. There is no X^{5/6} term.

**What it means.**
- Hasse's formula g(c)³ = μ(c)c²c̄ gives the exact identity μ(c) = g̃(c)³·λ̄(c), with λ(c) = c/|c|, for every c ≡ 1 mod 3. So
  Möbius against Gauss sums is the sum of g̃⁴λ̄. Numerically it behaves like a square-root-cancelling sequence, not like the
  Gauss-twisted extremal sequences β_b g̃(b) (β ≥ 0) that make Heath-Brown's cubic large sieve sharp (Dunn–Radziwiłł Thm 1.4).
- Hence the dual sequences of LM (§56: Möbius times twisted cubic Gauss sums) lie outside the class where the (KN)^{2/3} term is
  attained, consistent with the flat LM numerics (§57).
- Pointwise square-root cancellation of the dual sums T(k) = Σ_n μ(n)g̃(k, n) is still not enough for LM (PROVED, exponent count).
  - LM needs Σ_{k≤K}|T(k)|² = K·N + O(N²) with K = N²/U > N, i.e. the dual mean square equal to its diagonal up to relative error U/N.
  - Pointwise |T(k)| ≪ N^{1/2+ε} only gives K·N^{1+ε}, which is the trivial bound N² for the original sum.
  - The needed asymptotic is an off-diagonal statement over k, and Poisson in k returns the original family (§3, Theorem B).

So the Patterson bias is not the obstruction for LM. The obstruction is the self-dual asymptotic large sieve.
Cross-check with the earlier sextic run (mmnum/gbias.py, gbias_1.5e5.log):
- Squarefree sextic Gauss sums grow like 0.17·X^{2/3}, the sextic theta pole at 1/2 + 1/6.
- The Möbius-weighted sums stay at |S|/X^{1/2} ≤ 0.28.
- Prime sums grow faster than X^{1/2}.
The cubic and sextic pictures agree. In the record window itself, mmnum/lmbig_1e6.log gives LM at U = 10⁶, s = 1.10, as
0.9435 ± 0.017 of the naive diagonal, i.e. 1.00 of the true (coprimality-weighted) diagonal 0.94.
**Why μ removes the bias (HEURISTIC), and what LM's dual exactly needs (PROVED reduction).**
- If g̃ were multiplicative on squarefree c, Σ μ(c)g̃(c)N(c)^{−s} would be 1/ψ(s), with ψ(s) = Σ g̃(c)N(c)^{−s} Patterson's series. It would
  then vanish at the theta pole s = 5/6 where ψ has its residue; −ψ′/ψ has a pole there, which is the prime bias. The numerics match
  this picture.
- In the sextic family, μ(n)g̃₆(n) = ε(n)·(n/|n|)·conj(g̃₃(n)) (§56; ε of finite order). So the LM dual is literally
  Σ_{N(h)≤K} |Σ_{n~N} c_n χ̄_n(h)|², with c_n = ε(n)(n/|n|)conj(g̃₃(n)), sextic row characters χ̄_n(h) and K = N²/U > N.
- LM at length N holds iff this equals K·Σ|c_n|² + O(N^{1+ε}·Σ|c_n|²), i.e. a sharp cubic large sieve with exact leading term for the
  single explicit bias-free sequence c_n. Heath-Brown's sieve carries the (KN)^{2/3} term, which exceeds N for every
  K > N^{1/2}. (OpenAI's sextic sieve has the same K + D + (KD)^{2/3} shape.)
- Dunn–Radziwiłł's corrected sieve (Prop. 9.1/9.2) applies to this kind of sequence, but saves only (KN)^{−ρ} with ρ < 10⁻⁴. That
  reaches K ≤ N^{0.5003}, far from K > N.
- So LM is equivalent to "the theta bias is the only obstruction to the sharp cubic large sieve" for this one sequence, in the
  rows-dominant range K ∈ (N, N^{1.5}]. de Faveri–Dunn–Hoffstein (arXiv 2607.07911) prove the bias IS an obstruction; whether it is the
  only one is open. Poisson in k returns the original family (self-duality, §3).
**Spectral (Rankin–Selberg) route to the dual asymptotic: where it breaks (PROVED structural).**
- Writing χ̄_n(h) = (h/n)₂·(h/n)₃ gives T(h) = Σ_n ε(n)λ(n)(h/n)₂·conj(g̃₃(h, n))·W.
- Without the factor (h/n)₂, T(h) would be the smoothed h-th Fourier coefficient of the K-type-1 cubic metaplectic Eisenstein
  series (Kubota). Σ_h |T(h)|² would then come from a regularised Rankin–Selberg unfolding ⟨|E_ℓ(·, s)|², E(·, w)⟩: main term from
  the pole at w = 1, no theta pole since ℓ = 1, and power-saving error. That is exactly the asymptotic large sieve LM needs.
- The quadratic factor (h/n)₂ is a character of n to modulus h. It changes with h, so T(h) is not the coefficient sequence of one
  automorphic form, and the unfolding is unavailable. Hasse–Davenport moves the factor into a sextic Gauss sum, but that brings
  back μ(n)g̃₆(h, n), the original LM summand.
- OpenAI handle the two parts separately: Heath-Brown's quadratic large sieve over k (with its (KN)^ε loss on the leading
  term) and the cubic theta reflection. The asymptotic version needs both at once, with exact leading term K‖c‖² and error
  N^{1+ε}‖c‖².
- For the quadratic part alone this is plausible: quadratic Gauss sums are explicit, so Poisson in k is self-improving. For the
  mixed quadratic–cubic family it is open.

## 70. Copies are the positivity shadow of an exact identity (PROVED reduction; re-derives §54 independently, adds the λ_min form of (b))
**The identity.** Every row u′ of norm ~U′ factors uniquely as u′ = v·k⁶ with v sixth-power-free, and ψ_{vk⁶}(n) = ψ_v(n)·1_{(n,k)=1}.
Hence, exactly,
  A(U′; N) := Σ_{u′~U′} |M_{u′}(N)|² = Σ_k F^{(k)}(U′/N(k)⁶; N),  F^{(k)}(V; N) := Σ*_{v~V} |M_v^{(k)}(N)|²,
where M^{(k)} is the Möbius sum restricted to (n, k) = 1.

**Copies.** The copies bound (slope 5/6) is this identity plus positivity: drop every block except k ~ K, then use the paper's
A ≪ U′^{1+ε}N at U′ = N.

**What an asymptotic would add.** Suppose:
- (i) A(U′; N) = D·U′N + O(U′^{1−η}N) for U′ ≥ N^{1+δ};
- (ii) two-sided asymptotics F^{(k)}(V; N) = c_k·VN + O(V^{1−η}N) at the intermediate scales V = U′/N(k)⁶ ∈ (U, U′].

Subtracting the blocks k < K and dropping k > K by positivity gives K·F(U) ≤ c·UKN + O(U′^{1−η}N). That is LM at scale U with
saving ≈ U′^{η} over copies, and §71 shows such savings move the record.
- The diagonal main terms match exactly: Σ_k c_k U′N(k)^{−6} = D.
- Positivity alone, i.e. lower bound 0 for the blocks k < K, reproduces copies exactly.
- The lower half of (ii) is the obstacle. A lower bound F ≥ diag − o(diag) at supercritical scales means the dual Gram matrix
  satisfies λ_min ≥ K − o(K) on the Möbius–Gauss vector. That is again an asymptotic large sieve in the rows-dominant range
  (§69), so the induction over scales cannot start from upper bounds.

## 71. Partial LM on the exact record pipeline (MODEL, c61/thLM.py)
**Setup.** The pipeline is the record's: detector-affine plus timbalud B = B0 + 7·10⁻⁴, all three row constraints, and the geometry
checks. A saving of U^η over copies on a window [s1, s2] lowers the long leg to long(t) − η for t ∈ [s1, s2]. R is recomputed exactly
as min_t max(short, long_η), piece by piece.

| hypothesis | θ | change |
|---|---|---|
| none | 0.8748391 | record (reproduced) |
| η = 0.001 on [1.10, 1.25] | 0.8747953 | −4.4·10⁻⁵ |
| η = 0.003 on [1.10, 1.25] | 0.8747080 | −1.31·10⁻⁴ |
| η = 0.01 on [1.10, 1.25] | 0.8744044 | −4.35·10⁻⁴ |
| η = 0.03 on [1.05, 1.30] | 0.8735620 | −1.28·10⁻³ |

- On the record pipeline dθ/dη ≈ −0.044, linear in η.
- Every positive saving over the sixth-power-copies bound on [1.10, 1.25] gives a new record. Even U^{0.001} is worth about forty times
  the last record step. The earlier "η ≥ 0.003" threshold (§55) came from the grid model, whose baseline sits 8·10⁻⁵ above
  the record.
- So the cheapest record-moving open input is any power saving over copies for the long Möbius mean square in that window.

**§71 addendum: the inverse slot cost on the exact record pipeline (MODEL, c61/thcI.py, fastRwc.py).**
Here A_I = 1 − 2q/c_I − (δ − 2q/c_I)·r, all else as in the record.

| c_I | 2 | 1.95 | 1.9 | 1.8 |
|---|---|---|---|---|
| θ | 0.8748391 | 0.8747323 | 0.8746166 | ≤ 0.8743642 |

At c_I = 1.8 the optimum h hit the search bound 0.80, so the true value is lower. dθ/dc_I ≈ 2.2·10⁻³ near c_I = 2.
**§71 addendum: large values cannot supply even U^{0.001} on the long leg (PROVED exponent count).** The long-leg count is over rows with
|M_u(N)|² ≥ N^{1+δ}, N = U^{t*} = U^{1.162}, δ = 0.467. Halász–Montgomery gives |B|·N^{1+δ} ≤ N(N + |B|·G), with
G = max_{u≠u′}|Σ_{n≤N} μ²(n)ψ_uψ̄_{u′}(n)|, a character sum of length U^{1.162} and conductor ≍ U². A gain needs G < N^δ = U^{0.543}.
- Burgess (r = 2) gives G ≤ U^{0.956}.
- Even Lindelöf for the pair family gives G ≈ N^{1/2} = U^{0.581} > U^{0.543}.

So no large-values argument through the pair Gram matrix helps at σ = 0.733 < 3/4, consistent with §3. The long-leg saving has
to come from the mean square (LM) itself.

## 72. de Faveri's optimal n-th order large sieve (arXiv 2610.04045, 2 Oct 2026): checked, does not move the record (CITED + PROVED exponent count)
**The theorem (Thm 1.1).** For n ≥ 3, over K ⊇ μ_n:
  Θ_n(A, B) ≪ (AB)^ε (A + B + A^{1−1/n}B^{2/n} + A^{2/n}B^{1−1/n}).
This is optimal for n = 3, 4 by de Faveri–Dunn–Hoffstein. For n = 6 it replaces BGL's (AB)^{2/3} by A^{5/6}B^{1/3} + A^{1/3}B^{5/6},
and for n = 3 it equals Heath-Brown's (AB)^{2/3}.

**Where large sieves enter the record chain (paper text).**
- Part II (§12–§20) uses only Goldmakher–Louvel's quadratic sieve (sharp) and Heath-Brown's cubic sieve (14.6), in Lemma 14.2 and
  near line 5780.
- The sextic sieve, Lemma 9.1 with (KD)^{2/3}, is used only in Part I's Prop. 9.2, the 11/12 step. The record chain supersedes
  that step with β* ≤ 7/8.
- Even there, the adversary's witness at r = 1/2 still gives 1 − δ/2. de Faveri improves only the r = 1 branch, from 4/3 − δ to
  7/6 − δ.

**Direct use on the record's legs (sextic family, rows U, length D = U^λ).** Σ_u|Σ_n c_nψ_u(n)|² ≤ (U + D + U^{5/6}D^{1/3} + U^{1/3}D^{5/6})‖c‖².
This is sharp only for λ ≤ 1/2.
- Marked inverse leg: the count is 5/6 + r/3 − δr + z(1/3 − 2q), optimised at z = 1 − r. That gives 7/6 − δ at p = ½, i.e. 0.700
  against A_I = 0.5935.
- Plain leg: 2m = 0.834 > 1/2, so the sieve is not sharp there.
- Long leg: U^{1/3+5t/6}, worse than copies' U^{1/6+5t/6}.
- Sixth plain moment at m = 0.417: loss 0.375, against the needed < 0.09.

So de Faveri's sieve improves Part I-type and generic zero-density bounds. His Cor. 1.5 gives exponent 0.752 at σ = 0.733 for
n = 6, against the record's 0.5935. It does not touch the theta-reflection inputs that bind at the record.

## 73. Two conditional records, certified in exact rational arithmetic (PROVED relative to the record's proof; c61/certLM.py, certcI2.py, sidechecks_rad.py)
**Setup.** The record's row-inequality certificate (handoff cert.py) is rerun with timbalud's B = B0 + 7/10000 − 10⁻⁵ and one open
input changed. The certificate is adaptive box subdivision of [0, 3/4] × [0, 1/2] in exact rationals. It uses the monotonicity of R
in (δ, p), which persists under both changes, and the bound R ≤ max(short(t), long(t)) for any t ∈ [1, 3/2]. All side conditions
(3.3)–(5.15) of the detector-affine handoff pass exactly at both points.

**(A) Partial LM.**
- Hypothesis: for s ∈ [11/10, 5/4], every long Möbius mean square the detector uses (paper (1.4) / Cycle61 long envelope) satisfies
  Σ*_u |M_u(U^s)|² ≪ U^{1 + (5/6)(s−1) − 1/1000 + ε}·U^s. That is U^{0.001} better than sixth-power copies.
- Then at θ = 87480/100000, h = 7953/10000 all three row inequalities hold:

| row | boxes | min slack |
|---|---|---|
| new, d = h | 121 487 | 2.7·10⁻¹⁰ |
| old, d₋ | 30 429 | 2.7·10⁻⁹ |
| old, d = 1/2 | 29 | 9.7·10⁻⁴ |

- So θ ≤ 0.87480 < 0.8748389, relative to the record's proof.

**(B) Inverse slot cost 39/20.**
- Hypothesis: the marked inverse moment holds whenever r + (39/20)z ≤ 1, in place of r + 2z ≤ 1.
- The capacity is capped by the supply ℓ/d₊ (0.2104; the cap starts at r < 0.5897 and never binds at the optimum).
- Then at θ = 87474/100000, h = 7971/10000:

| row | boxes | min slack |
|---|---|---|
| new, d = h | 100 543 | 3.5·10⁻⁹ |
| old, d₋ | 31 579 | 6.9·10⁻⁸ |
| old, d = 1/2 | 29 | 1.5·10⁻³ |

- So θ ≤ 0.87474.

**Scope.** These are conditional theorems relative to the record's proof structure: the detector count R, the high-side row inequalities,
the affine-in-d and junction steps, and timbalud's ledger, which is unaffected because B and Λ are unchanged. The hypotheses are open
(§54, §69, §32). Neither is a frontier move. They fix the smallest record-moving forms of the two cheapest inputs. The η = 0 sanity run of the modified script at the exact record point (θ = 8748389/10⁷, h = 7951/10000) reaches depth 40
at the binding corner (δ, p) ≈ (0.47013, ½) with an excess of 2·10⁻¹⁰. The record's model margin there is about 8·10⁻⁹, below
what corner-bound boxes resolve; the record itself was certified with Bernstein polynomials (task 13). A sanity run with slack
(θ = 0.874845, h = 0.7957, η = 0) certifies the new row (108 529 boxes, min slack 1.5·10⁻⁹; c61/certLM_rec2.log), so the modified
script reproduces the record pipeline at η = 0.
**§72 addendum: de Faveri's mechanism is the copies mechanism on the sieve side (CITED reading of 2610.04045 §1.3).**
- His A^{1−1/n}B^{2/n} term comes from writing the dual Poisson variable as a = d·e^n, with e ≤ (B²/A)^{1/n}. The e^n part acts trivially
  on χ_b, so the e-sum is a sum of exact copies and is bounded trivially.
- In the LM range N ∈ (U, U²), his bound gives Σ_u|Σ_{n~N} c_nψ_u(n)|² ≪ U^{1/3}N^{5/6}·‖c‖². This has the same long slope 5/6 = 1 − 1/n
  as sixth-power copies, with intercept 1/3 instead of the 1/6 that copies get from the paper's sharp critical Möbius moment.
- Conjecture (1.8) is proved optimal for n = 3, 4 by Gauss-sum bias (dFDH) and expected for all n. So slope 5/6 is the
  general-sequence optimum, and any improvement of the long leg must use the arithmetic of μ (cf. §48). Section 69 says which property
  is needed: the dual sequence ε(n)λ(n)conj(g̃₃(n)) carries no bias.

## 74. The frontier moved upstream: Gao's directional affine, 0.874836 (CITED: GTimG/QuasiRiemannTracker PR #62, merged 16 Oct 18:30Z; paper read as data)
**What happened.** Upstream merged "Directional marked affine source at 0.874836" (θ = 218709/250000, Zhongpai Gao). It builds on our
detector-affine (PR #56) and timbalud's radical count (PR #50), so it contains everything in our unsubmitted local 0.8748389 (9acccf3).
**The live frontier is therefore 0.874836**, and our local value is superseded by 2.9·10⁻⁶.

**Gao's new ingredients.**
- **Stronger Theorem 2.2 region.** λ = 19793/5000 = 3.9586 and B = 987047/10⁶, so the capacity at m = 0.417 rises from 0.03854 to 0.03866.
- **Lower threshold.** q = 839/1250 = 0.6712, the same direction as the q0 lever of §68.
- **LOW branch.** It is paid with the incumbent κ0 = 2(0.874840) − 1 (6κ0 in place of 9/2).
- **Width coupling in HIGH.** The HIGH payment keeps the total selected width S(b) tied to the actual slot width b (Lemma 1, conserved
  width). This replaces the fixed supply 33/1000 and the single-rate κ = q − 3ρ/2. It uses the two-piece common-factor rate
  γ(ν) = ν − max(4ρ − 1, (2 + 22ρ)/21), from the valuation coupling 4L_c ≤ 5(c + d).
- **Where the affine count is used.** Only on [.451, .49] × [.49, .5] and [.445, .451] × [.4975, .5], with h = 3181/4000.

**Effect on tonight's work.**
- The q0 lever of §68 and the "raise the window supply" idea are subsumed by Gao's construction.
- The conditional prices of §71 and §73 transfer with the same slopes (dθ/dη ≈ −0.044, dθ/dc_I ≈ 2.2·10⁻³), now measured from 0.874836.
- The remaining squeezes inside Gao's construction are of order 10⁻⁷ or less: re-running LOW with the new incumbent, and tightening the
  1/20000 and 1/200000 loss allocations. A drastic move still needs one of the open inputs of §71/§73.
**§74 addendum: Gao's HIGH comparisons are tight at both ends (PROVED arithmetic, exact rationals; Gao paper (8)–(10)).**
At ν = q = 839/1250 the four comparisons (10) have slack:
- b = 0.03401: 5.1·10⁻⁶ (γS, binding), 9.5·10⁻⁴, 5.4·10⁻⁴, 1.1·10⁻³;
- b = 0.04: 3.1·10⁻⁴, 1.4·10⁻³, 6.2·10⁻⁶ (rτ_sF, binding), 5.1·10⁻⁴.

The binding rate pieces are γ(ν) = ν − (4ρ − 1) (the c ≤ d case of (5)) at the small-width end, and timbalud's τ_s = 13/400 at the
large-width end. Since ∂H/∂B ≈ −4.48, a ΔB of 10⁻³ (θ about −8·10⁻⁶) would need about 4.5·10⁻³ more saving at both ends at once.
For example γ would have to rise from 0.118 to about 0.25. My pipeline with Gao's (λ, B) gives the model optimum 0.87483584, only
1.6·10⁻⁷ below his certified 0.874836. So the affine route is exhausted at the 10⁻⁷–10⁻⁶ level.
**§74 addendum: open-input prices measured from Gao's frontier (MODEL, my pipeline with λ = 19793/5000, B = 987047/10⁶; c61/gao_*.log).**

| input | θ |
|---|---|
| none (Gao's constants) | 0.87483584 (certified upstream 0.874836) |
| LM U^{0.001} below copies on [1.10, 1.25] | 0.8747923 |
| same, U^{0.003} | 0.8747054 |
| same, U^{0.01} | 0.8744030 |
| c_I = 1.95 | 0.8747293 |
| affine ceiling (detector capacity at the window cap 0.04; numerator uncapped, so an upper estimate of the gain) | ≥ 0.8748017 |

The slopes are unchanged from our point: dθ/dη ≈ −0.044, dθ/dc_I ≈ 2.1·10⁻³. The whole remaining affine/Theorem 2.2 game is worth at most
about 3.4·10⁻⁵, while a single U^{0.001} saving in LM is worth 4.4·10⁻⁵.

## 75. The zero-side bootstrap for the long leg is an exact fixed point (PROVED exponent identity)
**The bound.** By the explicit formula (§34), Σ_u|M_u(U^s)|² ≲ Σ_δ U^{R(δ)}·U^{s(1+δ)}·⟨|L′(ρ)|^{−2}⟩ plus the near-line bulk.
Grant the strongest averaged input on L′, the Gonek–Hejhal-type ⟨|L′(ρ)|^{−2}⟩ ≪ U^ε, which is open even for ζ. LM-partial at s
then needs R(δ) + sδ ≤ 1 + (5/6)(s − 1) − η for every δ.

**The identity.** The record's own count satisfies R(δ) = long(t*) = 1/6 + (5/6 − δ)t* at its crossing t*(δ). Hence
R(δ) + t*δ = 1 + (5/6)(t* − 1) exactly. At the binding class (δ = 0.467, t* = 1.162): R = 1/6 + 0.36633·1.162 = 0.59235, tδ = 0.54265, and the
sum is 1.13500 = 1 + (5/6)(0.162).

**Consequence.** Feeding the record's zero counts back into the explicit formula returns exactly the copies exponent at s = t*.
The long leg is self-consistent, and no information about L′ at zeros, however strong, improves it without a sub-R zero count at
δ ≈ 0.467 (cf. §62). This is the zero-side face of the self-duality (§3, Theorem B).
**§72 addendum: de Faveri–Dunn–Hoffstein (arXiv 2607.07911) gives lower bounds only (CITED reading).**
- They prove Ξ₃ ≫ (AB)^{−ε}(A + B + (AB)^{2/3}) unconditionally and the quartic analogue, using Patterson's cubic theta coefficients
  and Rankin–Selberg for quartic metaplectic Eisenstein series to show that theta coefficients are large on average.
- Their extremal sequences are the biased Gauss-twisted ones. They give no upper bound or asymptotic for bias-free sequences.
- So the input LM needs has not appeared anywhere: the ε-free sextic sieve in the rows-dominant range for c_n = ε(n)λ(n)conj(g̃₃(n)),
  with the row-dependent quadratic twist.
- The nearest available pieces:
  - Lindelöf on average over indices for Kubota series (dFDH);
  - the Rankin–Selberg asymptotic for the pure cubic K-type-1 case (no quadratic twist);
  - Heath-Brown's quadratic sieve with its (KN)^ε loss on the leading term.
- Each misses one of the three requirements: the exact leading constant, the mixed twist, or an error of size N^{1+ε}.

## 76. The dual sieve of §69, measured directly (NUMERICAL, dualsieve/ds2.py)

Setup. Columns are squarefree composite n built from split primes of distinct rational norm, with N(n) ∈ [P, 2P]. Gauss sums
g6(n) are computed directly over Z/N(n) with additive character e(d·Tr(n)/N(n)). Rows are primary h with N(h) ≤ K. The script
accumulates T(h) = Σ_n c_n Π_i conj((h/π_i)_6) and reports S/(rows·Σ|c|²) at each K. (An earlier version, ds.py, used prime-only
columns. Those carry the sextic theta bias of primes, so they were replaced.)

| sequence | P = 1000 (264 cols, K = 500 … 32000) | P = 3000 (778 cols, K = 1500 … 96000) |
|---|---|---|
| μ·g6 (the LM dual sequence) | 0.938 0.943 0.953 0.936 0.955 0.943 0.944 | 0.993 0.954 0.952 0.941 0.946 0.938 0.940 |
| g6 (no Möbius) | 0.920 1.005 0.986 1.023 1.038 1.064 1.088 | 1.061 1.048 1.052 1.067 1.082 1.089 1.088 |
| random unimodular | 1.082 0.980 0.945 0.931 0.957 0.954 0.955 | 0.963 0.930 0.929 0.929 0.940 0.942 0.946 |

Reading.
- The LM dual sequence μg6 = ε(n)(n/|n|)conj(g̃3(n)) sits on the coprimality plateau 0.94 (predicted Π(1 − 1/(N(p)(N(p)+1))) =
  0.9399), from the conductor scale K ≈ P up to K = 32P. It is indistinguishable from the random control. That is the shape LM
  predicts: exact leading term K‖c‖²·(local factor), with the off-diagonal inside the noise.
- The unweighted Gauss sequence g6 rises above the plateau, by about +0.15 at K = 32P at both sizes. That is the
  theta-bias excess: the large sieve is not sharp for it, in line with Dunn–Radziwiłł's lower bound for biased sequences.
- So Möbius removes exactly the obstruction that makes the general sextic sieve non-sharp. This does not prove anything, but it
  locates the difficulty: a proof has to use the absence of a pole in Σ_n Ψ(n)conj(g̃3(h, n))(h/n)₂N(n)^{−s} (angular twist
  ℓ = 1), uniformly in the row h.

## 77. Where LM sits: the thick-family analogue, and why first moments do not substitute (PROVED exponent comparisons; CITED classical results)
**Thick analogue.** Replace the U sextic characters of conductor ≍ U by the ≍ Q² primitive Dirichlet characters of conductor ≍ Q.
LM becomes Σ*_{cond χ ≍ Q} |Σ_{n≤x} μ(n)χ(n)|² ≪ Q²·x^{1+ε} for x > Q² (beyond the family size). In residue classes this is
Σ_{q≍Q} Σ_a |Σ_{n≤x, n≡a(q)} μ(n)|² ≪ Q·x^{1+ε} for Q ≤ x^{1/2−c}, a power-saving Barban–Davenport–Halberstam/Montgomery–Hooley bound
below the square-root barrier. Unconditionally the asymptotic is known only for Q ≥ x(log x)^{−A} (Montgomery, Hooley). Below that, the
error is x²(log x)^{−A} (Siegel–Walfisz plus large sieve). Under GRH the bound is immediate.
- Restricting to primitive characters of conductor ≥ x^δ removes the Siegel obstruction but does not help with L². Gallagher's sieve
  gives (Q + x/Q₁)x for the variance, i.e. x² in character form when Q₁ ≍ Q ≤ x^{1/2}. Vaughan plus the large sieve gives Type II
  pieces at LM quality only when both factors are ≤ Q, so only for x ≤ Q². This is Theorem D again.
- So LM is the sextic thin-family case of "GRH on average beyond the family size". The ζ (Matomäki–Radziwiłł) and thick-Dirichlet
  analogues are both open with power saving. The one structural difference in our family is that Poisson in u turns μ into Gauss sums
  (§56: μ(n)g̃6(n) = ε(n mod 4)(n/|n|)conj(g̃3(n))). In the thick family the dual of μ is a Chowla sum. Our dual is an asymptotic
  sieve for twisted cubic Kubota coefficients, blocked only by the row-dependent twist (h/n)₂.

**Where μ goes (NUMERICAL, g2per/g2per.py; corrects a first draft of this paragraph).** For squarefree primary n built from
1, 2 or 3 split primes (3522 elements, norms ≤ 6·10⁴), both the normalised quadratic Gauss sum g̃2(n) and
ε(n) = μ(n)g̃6(n)/((n/|n|)conj(g̃3(n))) are functions of n mod 4 alone: 12 classes, with no ambiguity mod 4, 12 or 36.
- So the dual coefficients are exactly c_n = ε(n mod 4)·(n/|n|)·conj(g̃3(n)): a periodic factor, a Hecke character and the
  conjugate cubic Gauss sum. No Möbius factor and no quasi-multiplicative sign survive. (An earlier draft said μ survives as
  the sign g̃2(n). That is false, because g̃2 is periodic mod 4, as ε_n is over Q.)
- μ is absorbed in the relation between the two Gauss-sum orders: μ(n) = ε(n)·(n/|n|)·conj(g̃3(n))·conj(g̃6(n)) on squarefree n. At a
  prime this is g̃3(π)³ = −π/|π| combined with Hasse–Davenport. The Jacobi sum J̃(π) = J(χ_π³, χ_π⁴)/√p is −(Weil's Hecke value) in
  this normalisation, which is consistent.
- So the only non-automorphic feature of the dual is the row character. The columns carry a cubic Gauss sum, but the row character
  conj((h/n)_6) = (h/n)₂·(h/n)₃ is sextic. Its cubic part is absorbed as conj(g̃3(h, n)), and its quadratic part (h/n)₂ is the
  obstruction.

**First moments do not substitute for LM (exponent check).** A count of large values from an L¹ bound Σ_u|M_u(D)| ≪ U·D^{1/2+ε}
(square-root cancellation in L¹) gives U·D^{−δ/2} rows above D^{(1+δ)/2}. Copies in L² give U^{1+(5/6)(s−1)}D^{−δ}. L¹ wins iff
δ < (5/3)(s − 1)/s, which is 0.232 at s = t* = 1.162. The binding class has δ = 0.467, so even a perfect L¹ input is worse than copies there.
Higher moments are worse too (k = 2 against k = 1 wins only for δ > 5/6; Theorem A).

## 78. Theorem 2.2's constants also come from copy efficiency 1/6 (PROVED reading of Cycle61 (2.10)–(2.23), dl_c61/c61paper.tex)
The proof amplifies with one common width a = K_a·g + α, where g = (n₁ + n₂ + (9/2)z − 1)₊ is the excess over the old capacity. Then:
- **Low bins (δ ≤ q0).** Slots are deleted greedily until the old 9/2 capacity holds. The pointwise prime bound charges δ ≤ q0 per unit of
  deleted width, i.e. rate r = q0/(9/2) = 2q0/9 per unit of excess. The amplifier enlarges the excess by 6a (K_main in (2.14)–(2.15)), so
  the reserve is a(1 − 6r) − r·g (2.17). Nonnegativity forces K_a = r/(1 − 6r) = (2q0/9)/(1 − 4q0/3).
- **The 6a is the copies mechanism.** The amplifier's extracted alternatives are (w, w_o) = (λ_p, λ_p), (6λ_p, 0), (7λ_p, λ_p): the
  amplifier averages a row over u ↦ u·p⁶, which fixes every column character coprime to p. Its radius therefore grows by 6a for a gain
  of a: efficiency 1/6 again.
- **High bins (δ ≥ q0).** Lemma 2.1's density ρ(δ) gives the reserve 1 − ρ − δA − ā(5 − 6ρ) (2.21). The identity (2.22) cancels the b
  dependence and leaves B ≤ B0 = (1 − ρ0 + c)/(q0 + c) with c = K_a(5 − 6ρ0).

**Consequence.** If the amplifier had copy efficiency e instead of 1/6, then K_a = r/(1 − r/e), and B0 and Λ would improve. Neither
shortcut works:
- Lemma 2.1's ρ is already the best count on [q0, 3/4] (§66: max_p R(δ, p) = ρ(δ) to 10⁻⁴).
- Twisted copies do not amplify, because a zero of L(s, ψ_u) does not pass to ψ_uξ (§64).

So every cost on the record's critical path traces to one of two sources:
- the copy efficiency 1/6 behind α′ = 5/6, c_P = 6κ and K_a's factor 1 − 6r;
- theta reflection squaring the conductor of a mark, behind c_I = 2.

With these two mechanisms fixed, the affine route is exhausted at the 10⁻⁷–10⁻⁶ level (§74).
**§78 addendum: the copy group is fixed by the theta function, not chosen (CITED, OpenAI paper overview, oai5.txt l.238–271).** The base sum
averages smoothed cubic-theta coefficients over the completed indices c·n³. Poisson in the averaging variable has nonzero frequencies
u·a⁶ with u sixth-power-free, so every row comes with its sixth-power multiples built in. A family with better copy efficiency (cubic
twists, efficiency 1/3; quadratic, 1/2) would need a different theta series in the base sum. The quartic version over Q(i) was worked
through in earlier sessions (QUARTIC notes, (L³) stopped at a proved barrier). The classical theta function loses the sparse cubic
support c·n³ that gives Part I its 2/3. So, inside this architecture, 1/6 is structural.

## 79. Convexity in the dual beats copies on one half and fails on the other (PROVED exponent bookkeeping; HEURISTIC entire-ness of the twisted Kubota series)
Write OffDiag(K) = Σ_{n≠n′} c_n c̄_{n′} Σ_h χ̄_nχ_{n′}(h)W(N(h)/K) = (1/2πi)∫ K^w W̃(w) Σ_{n≠n′} c_n c̄_{n′} L(w, χ̄_nχ_{n′}) dw, move to Re w = 1/2,
take absolute values, and insert the balanced approximate functional equation (conductor ≍ D², split at D).
- **Direct half.** Σ_{N(h)≲D} N(h)^{−1/2}(|B_h|² + diag_h). Poisson in h with the paper's subcritical Möbius moment (r ≤ 1) gives
  Σ_{h~H}|B_h|² ≪ HD + D^{2+ε} for H ≤ D, so this half is ≪ D^{3/2+ε}. The MS error is then (U/D)·K^{1/2}·D^{3/2} = U^{1/2}D^{3/2}.
  That is strictly better than copies U^{1/6}D^{11/6}, by the factor (U/D)^{1/3}.
- **Root-number half.** The root number of χ̄_nχ_{n′} is conj(g̃6(n))g̃6(n′) times a reciprocity sign, and c_n·conj(g̃6(n)) = μ(n). So
  this half is S₂ = Σ_{N(h)≲D} N(h)^{−1/2}|M_h(D)|², after splitting into classes mod 36 to factor the sign. That is the original
  Möbius family at rows H ≤ D and length D, weighted towards small rows.
- **What S₂ needs.** Beating copies by η at length s needs S₂ ≤ D^{11/6 − (1/3+η)/s}, which is D^{1.5465} at s = t* = 1.162. Row band H
  needs A(H; D) ≤ H^{1/2}D^{1.5465}.
  - Copies suffice for H ≥ D^{0.86}.
  - Below that the saving needed over copies is 0.287s′ − 1/3 in exponent of H (D = H^{s′}). It grows without bound as H → 1, where
    the requirement is the pointwise bound |M_h(D)| ≪ D^{0.773}.
  - The proved strip gives D^{0.8748} only.

So convexity in the dual does not make a self-improving loop: LM-partial at s is reduced to stronger LM at every s′ ≥ s, and finally to
pointwise Möbius cancellation for small conductors, which is GRH-type. The direct half shows that the copies loss comes entirely
from the root-number (Möbius) half of the functional equation.
**§74 addendum 2: slack inside Gao's own certificate (CITED, Gao paper (10), (13); MODEL comparison).** The weighted-numerator comparison (13)
is certified as E ≤ ϑ − 10⁻⁷ by Bernstein polynomials. The HIGH comparisons (10) carry ≥ 10⁻⁶ after the 1/20000 and 1/200000 loss
allocations. My pipeline with his (λ, B) has its optimum at 0.87483584. So his construction leaves about 2·10⁻⁷ (0.874836 → ≈ 0.8748358)
if the 10⁻⁷ margin is reduced and θ re-chosen. Realising it would mean regenerating every θ-dependent certificate, because ℓ = 11/3 − 4θ,
x, y and κ∗ all move (Target836 has ≈ 300 Lean files). This is a tracker-level squeeze of someone else's construction, not a structural
gain, and it is not pursued.

## 80. The plain slot cost 6κ is the dual copies amplifier (PROVED reading, OpenAI paper §18.5–18.6, oai5.txt l.10020–10140, 10490–10560)
- **Mechanism.** The induction of Lemma 18.1 amplifies in the dual rows by h ↦ h·p⁶, p in a pool of norm Z^{ℓ*}. The paper's words: "the new row
  hp⁶ has norm O(Z^{K+6ℓ*}) … summing over h gives the factor |P|⁻¹ = Z^{−ℓ*}". This gains ℓ* and extends the rows by 6ℓ* = 2σ. The
  extracted errors have widths (w, w_o) = (iℓ_p, e_iℓ_p) for i = 1, 6, 7.
- **Where 6κ comes from.** Each child then carries an excess F_act ≤ 6(w + ℓ) + small, paid by deleting slots. Deleting width d costs κd in the squared
  exponent and lowers A + (6κ − 1)z by 6κd, so the induction region is n₁ + n₂ + 6κz ≤ M, i.e. c_P = 6κ. With an amplifier of row cost
  c·ℓ instead of 6ℓ, the region would be n₁ + n₂ + cκz ≤ M.
- **Why only p⁶.** G(u, hp⁶) = G(u, h) for p ∤ u, exactly. For k < 6, G(u, hp^k) = χ̄_u(p)^k·G(u, h) carries a phase that depends on the
  column. Averaging |Σ_u B(u)χ̄_u(p)^kG(u, h)|² over p leaves Σ_p χ̄_uχ_v(p), a prime character sum of length Z^{ℓ*} and conductor Z^{2a₀},
  far outside Pólya–Vinogradov range, so it needs GRH-type input. Unit multiples h ↦ εh are exact copies per residue class of u mod 36 at no row
  cost, but they give only a factor 6, not a power.
- **Collet's Theorem 2.2** (§78) lowers the effective weight to Λ ≈ 3.96 in a window by combining the same copies amplifier with Lemma 2.1's density
  in high bins. It does not change the amplifier.

So the critical-length targets with substantial value, c_Pκ → 2 (≈ −1.5·10⁻³ in the grid model) and c_I → 1 (≈ −2·10⁻³), are blocked in the
same two places as LM: replication only by sixth powers, and theta reflection squaring a mark's conductor. §31's numerics say their true cost is
1, so an unconditional proof needs an amplifier that is not a sixth-power copy. The only candidate shifts h ↦ hp^k bring column phases χ_u(p)^k
whose average is a short prime character sum.

## 81. The HIGH-side sixth moment needs the central block (PROVED exponent bookkeeping; CITED Xu–Sang arXiv 2610.15291 §5, Prop. 5.1)
Xu–Sang bound Σ_{u∈B}|L(1/2 + ϑ + it, ψ̄*_u)| block by block: AFE blocks S_m of length U^m carry D^{−ϑ}. For m ≤ 1/2, two copies of S_m plus prime
slots of length z_p = v(1 − 2m) fit the plain fourth moment (2m + 6κz_p ≤ 1). Hölder gives (3R+1)/4 − qz_p/2 − ϑm = (3R+1)/4 − ϑ/2 for every m
once ϑ = qv (their (23)).
- A sixth moment with loss λ improves the HIGH exponent by d(1 − R − 2λ)/12. That is up to d(1 − R)/12 ≈ 0.027 (θ about −7·10⁻³) at λ = 0, and it
  pays iff λ < (1 − R)/2 ≈ 0.20.
- Block by block, the sixth moment gives 5R/6 + (1 + λ(m))/6 − ϑm, with λ(m) = (5/6)max(0, 3m − 1) from copies (granting a three-plain-factor
  critical lemma, which the paper does not have; Lemma 18.1 has two plain factors).
  - Short blocks m ≤ 1/3 come out at 0.66 − ϑm, against 0.694 − ϑ/2. That is about equal at ϑ ≈ 0.06.
  - The central block m = 1/2 has no slot room and λ = 0.417, giving 0.729 − ϑ/2 > 0.694 − ϑ/2. Fractional moments between 4 and 6 do not
    help either (dK/dk = (2R − 0.332)/16 > 0 at k = 2).
  - So the gain needs the sixth moment of the central block of length U^{1/2}, i.e. the sixth moment of L(1/2, ψ_u) itself (r = 6, §68).
- Unbalancing the approximate functional equation does not help: the dual blocks have length up to U^{1 − m₀} and weight U^{−ϑ(1−m)}, which is
  worse than ϑ/2.
- The shift ϑ = q/(6κ) is itself capped by the slot cost (v ≤ 1/c_P), so the HIGH side also inherits the copy efficiency through ϑ.

## 82. The central block's sixth moment shows no loss (NUMERICAL, mmnum/sixth.py, sixth2.py; all sixth-power-free rows N(u) ≤ U)
S_u(X) = X^{−1/2} Σ_n χ_n(u) W(N(n)/X), the centrally normalised plain block, over all ideals n (not only squarefree). The control is a random completely
multiplicative ±1 twist of the same block.

| U (rows) | X | m₄/m₂² plain / control | m₆/m₂³ plain (all rows) | m₆/m₂³ plain, Θ rows removed | m₆/m₂³ control |
|---|---|---|---|---|---|
| 10⁴ (36 234) | √U | 1.905 / 1.914 | 5.241 | – | 5.276 |
| 10⁵ (362 214) | √U | 2.041 / 2.050 | 6.462 | 6.327 | 6.426 |
| 10⁶ (3.6·10⁶) | √U | 2.085 / 2.084 | 7.236 | 6.755 | 6.790 |
| 3·10⁶ (1.09·10⁷) | √U | 2.110 / 2.111 | 7.808 | 6.984 | 7.040 |
| 10⁶ | √U/4 | 1.999 / 1.999 | 5.985 | – | 5.980 |

- The only outlier is the row u = 1 (principal character): |S|²/m₂ = 37, 120, 208 at U = 10⁵, 10⁶, 3·10⁶. It is the fixed-family (Θ) row that
  Lemma 18.1 and Proposition 5.1 already exclude. Removing the 216 rows with {2,3}-smooth norm puts the plain block exactly on the control.
  The slow rise of both columns is the shared arithmetic factor of three-fold products, not a power.
- So input 4 (FOUR_INPUTS.md), the sixth moment of the central approximate-functional-equation block that the HIGH side needs, has no visible loss
  (λ ≈ 0) up to U = 3·10⁶. Copies certify only λ = 0.417, and λ < 0.20 is what moves the frontier.
- The block of length U^{1/2} has cube length U^{3/2}, i.e. 1.73·10⁹ multiplied out at U = 3·10⁶ against 1.09·10⁷ rows: far past the
  conductor, like LM, and as flat as LM.
**§82 addendum: what unconditional bounds give for the central block (PROVED exponent bookkeeping).** The fourth moment (Lemma 18.1, any lengths)
together with a Burgess-type pointwise bound for a character sum of length √(conductor) (|S|² ≪ U^{3/8} after normalisation) gives λ ≤ 3/8 by
Hölder. That is slightly better than copies (5/12), but the HIGH side needs λ < 0.20. Three plain factors (1/2, 1/2, z) cannot be reflected
below total length 1 + z, so even a short third factor is past the conductor, by z. The central block leaves no room for slot amplification
(2·(1/2) + c·z ≤ 1 forces z = 0). Halász-type pretentiousness constrains only o(X) cancellation, not square-root-scale large values.

## 83. The plain slot cost reduces to twisted fourth moments with a partial saving (PROVED reduction + exponent bookkeeping; HEURISTIC Θ estimate)
Write the marked plain moment as Σ_u|S₁S₂Q|² = Σ_{m,m′} q_m q̄_{m′} T(m, m′), where T(m, m′) = Σ_u |S₁S₂(u)|² ψ_u(m)ψ̄_u(m′), m and m′ run over slot products of norm
≤ U^z, and Q is normalised so that mean|Q|² ≍ 1.
- **Diagonal.** m = m′ gives ‖q‖²·Σ_u|S₁S₂|² ≪ U^{1+ε} for any lengths (Lemma 18.1, unmarked).
- **Reduction.** If T(m, m′) ≪ U^{1+ε}N(mm′)^{−θ′} for m ≠ m′ (with n₁ + n₂ + z ≤ 1), then Σ|S₁S₂Q|² ≪ U^{1+(1−2θ′)z+ε}. The plain-leg count becomes
  A_P = 1 − 2δm − (2q − 1 + 2θ′)z with z ≤ 1 − 2m (and the supply ℓ/d). Lemma 18.1 gives A_P = 1 − 2δm − 2qz with z ≤ (1 − 2m)/(6κ). The new
  bound wins iff 2q − 1 + 2θ′ > 2q/(6κ), i.e. θ′ > 0.318 at the binding class (q = δ/2 = 0.233), or θ′ > 0.325 against Theorem 2.2's Λ = 3.96.
  θ′ = 1/2 (cost 1) would move A_P at m* = 0.417 from 0.5935 to about 0.534.
- **The Θ rows do not block it (HEURISTIC estimate).** After Poisson, the marked dual is B^{(Q)}_h = Σ_m q_m g̃6(h, m)·σ·B_{hm⁴} (twisted multiplicativity,
  with g̃6(k, n′)(m/n′)₃ = g̃6(km⁴, n′)). Shifted Θ rows km⁴ are generic, so the bias is removed. The copies of a Θ row k with aligned slots
  contribute about (U/D)·K^{1/6}|B_k|²·Z^z ≈ U^{5/6}D·Z^z, i.e. U^{z−1/6} relative to the diagonal: harmless for z < 1/6. This agrees with §31
  (no excess even for a(p) = 1).
- **Why it is still open.** Poisson in u gives T(m, m′) = (U/D)N(mm′)^{−1/2}·Σ_{h≤D²N(mm′)/U} χ̄(h)·B_h[a·ψ(m)]·conj(B_h[a·ψ(m′)]). Its dual row range
  exceeds the column length by U^z, so it is past the conductor by exactly the slot length. Cauchy plus the sieve gives θ′ = 0, the trivial bound.
  Lemma 18.1 cannot reflect the slots, because prime-supported factors are not reflectable; that is why its proof uses the dual copies amplifier (§80).
- **So input 2 in its sharpest form** is a near-conductor (by z ≈ 0.04–0.17) twisted fourth moment for plain products with saving θ′ > 0.32 in the twist
  modulus. This is weaker than LM in two ways: the coefficients are plain products, which reflect and whose Θ main terms are explicit, and only
  two thirds of the square-root saving is needed.
**§83 addendum: the twisted fourth moments save at least N(mm′)^{−1/2} numerically (NUMERICAL, mmnum/twist4.py; all sixth-power-free rows).**
T(m, m′) for distinct prime slots m, m′, with S₁ = S₂ the plain block at length U^n (2n + z ≤ 1). θ_eff = log(|T|/Σ|S₁S₂|²)/log(N(mm′)^{−1}):

| U | n | slot norm ≈ U^z | pairs | median |T|/Σ|S₁S₂|² | median N(mm′)^{−1/2} | θ_eff |
|---|---|---|---|---|---|---|
| 10⁶ | 0.35 | U^{0.30} | 15 | 0.0013 | 0.0210 | 0.86 (row-noise floor 0.0005) |
| 3·10⁶ | 0.35 | U^{0.30} | 28 | 0.0005 | 0.0144 | 0.90 (noise floor 0.0003) |
| 10⁶ | 0.30 | U^{0.40} | 276 | 0.0006 | 0.0055 | 0.71 (noise floor 0.0005) |
| 3·10⁶ | 0.40 | U^{0.20} | 6 | 0.0367 | 0.0636 | 0.60 |

The binding geometry (slot norm U^{0.04–0.08}) is out of reach at these U, since slot primes would have norm below 10. In every feasible range the
twisted moment sits at or below N(mm′)^{−1/2}, with θ_eff ≥ 0.6 ≥ 0.32. So the partial saving that §83 needs is visible in the data, while
Poisson plus the sieve proves only θ′ = 0.
**§83 addendum 2: the anatomy of the 6, and cost 2 as a level-of-distribution statement (PROVED reduction).** Take single slot primes m ≠ m′ of norm Z
and the plain product P = S₁S₂ of length D.
- **Upper bound.** Poisson in u gives T(m, m′) = (U/(DZ))·Σ_{h≲(DZ)²/U} B^{(m)}_h·conj(B^{(m′)}_h)·σ, with B^{(m)}_h = g̃6(h, m)·Σ_α σ_α B^α_{hm⁴} (twisted
  multiplicativity, plus g̃6(h, n)(n/m)₃ = σ·g̃6(hm⁴, n), classes α mod 36). Cauchy then reduces T to Σ_{h≤H}|B_{hm⁴}|², the unmarked dual on a
  sublattice of index Z⁴.
- **Dropping the sublattice gives 6.** Positivity and the Poisson-dual form of Lemma 18.1 (with V̂ ≥ 0 weights) give Σ_{h′≤HZ⁴}|B_{h′}|² ≪ (D + HZ⁴)‖c‖²U^ε.
  So T ≪ U/Z + DZ⁵, which is θ′ = 1/2 exactly when n₁ + n₂ + 6z ≤ 1. The paper's 6κ is this (κ ≤ 1 from the pointwise cap).
- **Equidistribution would give 2.** If |B_{h′}|² were equidistributed over h′ ≡ 0 mod m⁴, then T ≪ U/Z + DZ, i.e. cost 2 (n₁ + n₂ + 2z ≤ 1). That is c_Pκ ≈ 1.5
  against 4.5 now.
- **The same question on the original side.** ψ_u(n)·(n/m)₃ = σ·ψ_{um²}(n) by reciprocity, so B_{hm⁴} is the dual of the plain product at rows u·m².
  Cost 2 is equivalent to Σ_{u≡0 (m²), N(u)≤U}|S₁S₂(u)|² ≈ Z^{−2}Σ_{N(u)≤U}|S₁S₂(u)|² on average over slot primes. That is a level of distribution U^{2z}
  (≈ U^{0.08–0.17}) for the row function |S₁S₂|².
- **Why this stops.** The large sieve over characters mod m² turns it into the eighth moment Σ_u|S₁S₂|⁴. Its cube-type length is U^{1.668}, past the
  conductor, and copies lose U^{0.56} there. This is Proposition E's route again (DH_NOTE §7).
**§83 addendum 3: final form of input 2, a fourth moment over mildly thin subfamilies (PROVED reduction).**
- **The subfamily.** The rows u = vm² (v ≤ U/Z²) form the family {ψ_{vm²}} = {σ·ψ_v·(·/m)₃}. Its size is U/Z², but its conductor is ≍ N(v)·N(m) ≍ U/Z: it
  is thin by a factor Z = U^z.
- **What Lemma 18.1 gives.** With moving residue-symbol factors of norm Z^q, the bound is Z^{m+q}. That is the conductor-sized bound U/Z, a loss of Z
  against the natural U/Z², and it is exactly what yields θ′ = 0, and hence 6κ after the amplifier.
- **What is needed.** Cost 2 needs Σ_{v≤U/Z²}|S₁S₂(vm²)|² ≪ (U/Z²)^{1+ε} (natural size) for n₁ + n₂ + 2z ≤ 1. Beating 6κ needs only a loss below Z^{0.36}
  (θ′ > 0.32), on average over the slot primes m.
- **Why it is the same wall.** Poisson in v has dual rows D²Z²/U ≤ D (columns-dominant), but its coefficients c_n(n/m)₃g̃6(n) are the unmarked dual at rows
  h·m⁴. Positivity there returns cost 6, and the general sieve needs K ≤ D^{1/2}. So input 2 is a thin-family Lindelöf-on-average statement, the same type as
  the quartic (L³) problem (thinness Q² there). Here the thinness is only U^{0.04–0.08}, which makes it the mildest instance of that type in the whole
  architecture.

## 84. Why every attack on input 2 meets self-duality (PROVED structure)
For distinct slot primes m ≠ m′, the twisted moment T(m, m′) = Σ_u|S₁S₂(u)|²ψ_u(m)ψ̄_u(m′) has no main term in either representation.
- **Primal.** The diagonal (Poisson frequency h = 0) vanishes, because u ↦ ψ_u(m)ψ̄_u(m′) is a nonprincipal character mod mm′.
- **Dual.** The row-diagonal of Σ_h B^{(m)}_h conj(B^{(m′)}_h) is Σ_n|c_n|²·Σ_{h≤H, (h,n)=1} χ_mχ̄_{m′}(h)·(units). That is O(Z‖c‖²) by orthogonality of a
  character of modulus mm′, not H‖c‖².
- **The identity.** Cauchy–Schwarz (which gives 6, §83) throws exactly this cancellation away. Without it, the dual off-diagonal is T itself: Poisson in h
  returns T with coefficient (U/(DZ))·(H/(DZ)) = 1.

So any θ′ > 0 must break self-duality, as for LM (§3, Theorem B). The data (§83 addendum) put |T| at the row-noise level √(rows) × typical, far
below N(mm′)^{−1/2}.
- **Multiple Dirichlet series do not help.** The plain-product dual coefficients Σ_{d,e} g(h, de)N(d)^{−s₁}N(e)^{−s₂} are partial sums of the sextic A₂
  Weyl group multiple Dirichlet series (twisted multiplicativity gives the A₂ cocycle (d/e)₆(e/d)₆). Their mean square over the index h is Poisson-dual to
  the plain fourth moment, i.e. affine D̃₄ (§68). That type continues only to a half-space, which supports the paper's upper bound but not a twisted
  asymptotic with power saving.
- **What it would take.** The Dirichlet-character analogue of input 2 (twisted fourth moments) is proved through Kuznetsov's formula, Kloosterman sums and
  Weil's bound. The sextic version would need the Kuznetsov formula on the cubic/sextic cover with metaplectic Kloosterman sums, plus the
  twisted-multiplicativity bookkeeping of §83.
**§83 addendum 4: the plain slot weight priced on the current record pipeline (MODEL, c61/cp/price.py; Gao's constants, baseline 0.874836
reproduced).** The slot weight c (now 6κ ≈ 4.498) is replaced in both the plain leg and the HIGH capacity C_O. For c < 4.5 the affine branch is replaced by
the same weight everywhere (B = 1, no window):

| c | 4.498 (no affine) | 3.96 | 3.0 | 2.0 | 1.0 |
|---|---|---|---|---|---|
| θ | 0.874887 | 0.874732 | 0.874313 | 0.873382 | 0.870521 |

So input 2 at θ′ = 1/2 (weight 2) is worth −1.45·10⁻³ against the frontier, and natural cost 1 is worth −4.3·10⁻³. Even weight 3.96 everywhere with
B = 1, i.e. Collet's Λ without his window and without B₀ < 1, would be −1.0·10⁻⁴.
- **The window is not the issue; B₀ < 1 is (MODEL, cp/price2.py).** With Gao's (B, Λ) = (0.987047, 3.9586), widening Theorem 2.2's window from z ∈ [0.034, 0.04],
  n ∈ [0.41, 0.44] all the way to z ∈ [0.001, 0.2], n ∈ [0, 0.5] leaves θ at 0.8748359. The −1.0·10⁻⁴ for weight 3.96 above comes entirely from B = 1.
  And B₀(q0) = 1 − (ρ(q0) − (1 − q0))/(q0 + c): its deficit is the high-bin density's excess over DH at δ = q0 ≈ 0.67 (0.061, divided by q0 + c ≈ 5.07).
  That density is Lemma 2.1's ρ, set by the copies long leg t* = 5/(5 − 2δ) (§66). So pushing B₀ toward 1 is input 1 (LM) again, at δ ≈ 0.67–0.75.

## 85. Why the inverse slot cost is 2 but the plain one is 6κ (PROVED reading, OpenAI Lemmas 17.1, 17.2, 18.1)
- **Lemma 17.1 (marked inverse).** Σ_u|M_u(Z^r)Q_u|² ≪ Z^{m+ε} for r + 2z ≤ m − c₁ and 2r + 8z ≤ 3m − c₂. It is proved through the canonical family (Lemma 17.2): rows k
  and squarefree columns n with coefficients ᾱγ₂νρ(n)·d(n) and row character χ_n(k)χ_n(f)⁴ averaged over f. This is the Poisson dual of μ·(sextic Gauss)
  = ε(n mod 4)(n/|n|)conj(g̃3(n)) (§77). The cubic structure can be reflected, and marks cost 2 (conductor q², saving q).
- **The canonical family's restriction.** Its column coefficient d must be a product of short prime slots of total length ≤ z₀: "no additional residual
  coefficient b(n), even one independent of k, f".
- **Lemma 18.1 (marked plain).** Poisson in u turns the plain product S₁S₂ into a_n·g̃6(n). On squarefree n this is a_nμ(n)·ε(n mod 4)(n/|n|)conj(g̃3(n)) with
  a_nμ(n) = Σ_{n=ab} μ(a)μ(b): the dual of a product of two long inverse polynomials of lengths U^{n₁}, U^{n₂}. That is not admissible in Lemma 17.2,
  so the paper works on the plain side with reflection plus the dual copies amplifier, at cost 6κ (§80, §83).
- **Consequence.** Bringing the plain leg to cost 2 by the inverse route would need a canonical-family estimate with a long second factor in the
  coefficients.
  - The marked Möbius-times-Möbius moment Σ|M₁M₂Q|² is itself at critical length (r₁ + r₂ + z ≤ 1). Its dual coefficients are divisor-weighted cubic
    Gauss sums: (μ∗μ)(n)g̃6(n) = 2^{ω(n)}ε(n mod 4)(n/|n|)conj(g̃3(n)) on squarefree n. These are positive weights, outside Lemma 17.2's admissible
    coefficients.
  - Through Proposition E, it needs the Möbius fourth moment, which is past the conductor for r > 1/2 (task #6).
  - So the asymmetry is explained, but the inverse cost cannot be transferred to the plain leg.
**§83 addendum 5: the shift identity is exact (NUMERICAL, g2per/shift4.py).** For 58 random pairs of split primes n ≠ m (norms ≤ 200) and 2–4 random rows h each:
g̃6(h, nm) / (g̃6(h, m)·g̃6(hm⁴, n)) is a sixth root of unity, independent of h, to 2·10⁻¹⁵. This confirms B^{(m)}_h = g̃6(h, m)·Σ_α σ_α B^α_{hm⁴}, the identity behind
the 2 + 4 anatomy of the plain slot cost.

## 86. A cubic-character base probe would lose the Möbius function (PROVED reading, OpenAI §6.1, oai5.txt l.2684–2780)
The corrected row is T_{m,s,η}(t) = Σ_{c sf, n} γ₂(c)α(cn³)Ψ_{m,s,η}(cn³)G(cn³)N(c)^{−1/2−t}N(n)^{−1−3t}, with Ψ ∋ χ_a(m) = (m/a)₆. Its sextic symbol splits as
χ_{cn³}(m) = χ_c(m)·(m/n)₂.
- Poisson in m turns χ_c(m) into the sextic Gauss sum g̃₆(c). Through μ(c)g̃₆(c) = ε(c mod 4)(c/|c|)conj(g̃₃(c)) (§77), the γ₂α-weighted cubic-theta coefficients
  become the Möbius function: this is the source of 1/L_F(s, η) in the Poisson rows (paper §7.2).
- A cubic symbol (m/a)₃ would produce g̃₃(c) instead. Since g̃₃(c)³ = μ(c)(c/|c|), a cubic Gauss sum is a cube root of μ times a Hecke character, not μ
  itself, so the reciprocal would be lost.
- So the sextic family is forced by the detector's need for 1/L, and with it the Poisson frequencies u·a⁶ and the copy efficiency 1/6 (§78).
  Replacing the family is not an available lever.

## 87. The LOW lever: price, source of b/12, and why it reduces to §12.4 (MODEL c61/cp/pricelow.py; PROVED reading OpenAI Prop. 15.2–15.3, oai5.txt l.6293–6500)
**Price (MODEL, record pipeline with Gao's constants).** A saving Δ in the LOW constraint (θ ≥ 11/12 − ℓ/4 − Δ, i.e. ℓ = 11/3 − 4θ − 4Δ):

| Δ | 0 | 0.0005 | 0.001 | 0.003 | 0.01 |
|---|---|---|---|---|---|
| θ | 0.8748359 | 0.8744495 | 0.8740631 | 0.8725176 | 0.8671084 |

dθ/dΔ ≈ −0.77, the steepest slope of any lever found. At the record b = y − x = 2h − 3ℓ − 1 ≈ 0.0885, so the whole b/12 is worth Δ ≈ 0.0074.
**Source of b/12 (PROVED reading).** Prop. 15.2's Gram bound is (Q/Y′)(1 + P_a^{1/6} + P_a²/Y′):
- the diagonal k = 0 gives Q/Y′;
- the nonexceptional frequencies give (Q/Y′)P_a²/Y′ by 4-dimensional lattice Poisson;
- the exceptional frequencies (k) = a₁⁶e, counted O(K_R^{1/6}q_C^ε) and bounded with |L_C| ≤ q_C and the trivial point count, give (Q/Y′)P_a^{1/6}. The C = d = 1, R = 1 dyad
  dominates.
Prop. 15.3 has ly − d − 11b/6 ≥ 1/12, so the exceptional term dominates. Cauchy–Schwarz with Lemma 15.1 (Σ|B_m|² ≪ Z^{M′}) gives (X′)^{1/2}P_a^{1/12}, the b/12.
**Why it is not an independent lever (REPORT §12.2–§12.4, re-checked).** The exceptional frequencies form a fixed translate k₀·a⁶ with nonzero complete
mean (fixed sign), so the term is a genuine energy: the twisted sextic theta coefficients r_σ(m) of A_m's polar part. Their pairing with B is
J = κS(Z) + (row terms), with S(Z) the η-Möbius sum. So any bound better than Cauchy–Schwarz on the theta-paired part bounds S(Z), which is the zero-free
statement itself. The LOW side is a theta-dual large sieve for S(Z) with alignment loss Q^{1/12}. Removing b/12 gives no input beyond the theorem, and the price
above is not reachable through the LOW side. What remains is the row energy below its threshold, i.e. the cap M + ℓ ≤ 1 (REPORT §13.2, Möbius
fourth moment at r ≈ 0.7, §14).
**HIGH-side constants checked in the same pass (PROVED reading, Lemma 7.1 (7.14), (20.4)).** 13h/75 + d(1/4 − 17/50) = z_r(h − d) − h/6 + d/4 with
z_r = 17/50, the lower edge of region D1. Its coefficient vanishes at the binding d = h, so z_r is not a lever (as in §68).
**Input 4 against the ζ analogue (CITED).** The best unconditional sixth moment of ζ is ∫_T^{2T}|ζ(1/2+it)|⁶dt ≪ T^{5/4+ε} (Hölder between the fourth
and Heath-Brown's twelfth moment), a loss of 1/4. Input 4 needs a loss λ < 0.20 for the sextic family's central values. So it is beyond the ζ state of the art,
as LM is (§36). Even a twelfth moment ≪ U^{2+ε} for the family would give only λ = 1/4 by the same Hölder step.
**Detector count against DH at the binding (MODEL, fastRw.R with Gao's B, Λ).** R(δ, 1/2) − (1 − δ) = 0.056, 0.058, 0.0595, 0.060, 0.060, 0.058, 0.044
at δ = 0.40, 0.43, 0.467, 0.50, 0.55, 0.60, 0.70. Every substantial lever acts on this roughly uniform excess of about 0.06 over DH. It equals (5/6 − δ)(t* − 1), with the
crossing at t* ≈ 1.16. The count 0.592 at σ = 0.7335 is already below the Montgomery-type 3(1−σ)/(2−σ) = 0.631 and the Guth–Maynard-type (15/13)·2(1−σ) = 0.615.

## 88. A large sieve over the slot twists: the slot cost is free given M4, and the unconditional square-root version loses (PROVED inequality; MODEL exponent count)
Write w_u = |M_u(N)|², Q_u = Z^{−1/2}Σ_p a_pψ_u(p) with |a_p| ≤ 1, and T(p, p′) = Σ_u w_u ψ_u(p)ψ̄_u(p′). By reciprocity, ψ_u(p)ψ̄_u(p′) is a character of u modulo pp′,
so T is a twisted Möbius second moment with twist conductor ≤ Z². Then:
- Σ_u w_u|Q_u|² = Z^{−1}Σ_{p,p′}a_pā_{p′}T(p, p′) ≤ Σ_u w_u + (Σ_{p≠p′}|T(p, p′)|²)^{1/2}.
- The classical large sieve for the ≤ Z² characters of modulus ≤ Z² in the variable u gives Σ_{p≠p′}|T|² ≤ (Z⁴ + U)Σ_u w_u².
- Hence, for Z⁴ ≤ U: **Σ_u|M_u(N)Q_u|² ≤ Σ_u|M_u(N)|² + (U·Σ_u|M_u(N)|⁴)^{1/2}.**
Consequences:
- Given the Möbius fourth moment at length N (Σ|M_u|⁴ ≪ U^{1+ε}), slots of any length z ≤ 1/4 are free (c_I = 0 for z ≤ 1/4). This sharpens §33 (cost 1 from M4).
  M4 at r ≈ 0.75 is a super-DH statement (REPORT §14.6), so this changes no unconditional value.
- With a fourth-moment loss L₄ the bound is U^{1+L₄/2}. It beats Lemma 17.1 at the binding r = 0.75 (δ = 0.467, q = 0.2335) only if L₄ < 0.117. The best
  available is L₄ = 0.75 (the trivial |M|² ≤ N), or 0.417 via copies, which would need M4 at critical length in the inflated family (known only for r ≤ 1/3).
- Restricted to the saturated rows B₁ (|B₁| ≤ U^{1−rδ}) and with Q^k, Z^k = U^{1/4}, the same sieve counts rows with large slots by U^{1 − rδ/2 − q/2} = U^{0.708},
  against Lemma 17.1's U^{1 − rδ − q(1−r)} = U^{0.592}. The square root in the Hilbert–Schmidt step costs more than the reflection loses.
- The same inequality for the plain leg needs the plain eighth moment at total length 2(n₁ + n₂) ≈ 1.67, again supercritical.
So the slot costs sit exactly between the critical second moment (which Lemmas 17.1 and 18.1 handle with cost 2 and 6κ) and the supercritical fourth moment (which makes
them free). That is the "common root" of FOUR_INPUTS in one inequality.

## 89. Joint prices of the open inputs on the record pipeline, with a correction to this session's quick price scripts (MODEL, c61/cp/combo_wa.py, combo_ws.py)
**Correction.** c61/cp/price.py, price2.py and pricelow.py (written this session) did not enforce the geometry and LOW-validity conditions that
thoptw_core.margin, thLM.py and the certificates do (ℓ < x < y < 2x, y < 1/2 < h < 2y, ℓ < 1/5, h > 4ℓ, 4ℓ − y > 0; LOW form (5.5) needs
ℓ < x ≤ y). c61/cp/combo_wa.py enforces them, searches h ∈ (0.60, 0.82), and for plain weights c < Λ uses that weight everywhere (B = 1, wide window).
Re-run with Gao's (B, Λ), baseline 0.874836 reproduced:
- unchanged (the old optima were valid): plain weight 3 / 2 / 1 → 0.874313 / 0.873382 / 0.870524 (§83 addendum 4); LOW saving Δ = 0.001 / 0.003 / 0.01 →
  0.874063 / 0.872518 / 0.867108 (§87); LM slope α = 0.80 / 0.75 / 0.70 → 0.874589 / 0.874175 / 0.873662;
- changed: c_I = 1 → 0.870352. The unconstrained 0.870184 sat at x > y.
The certified conditional records of §73 and the partial-LM table of §71 were produced by scripts that enforce all conditions, so they are unaffected.
**Joint prices (Gao's constants, cap M + ℓ = 1, no slot-supply cap):**

| c_I | LM slope α | plain weight | θ |
|---|---|---|---|
| 2 | 5/6 | 6κ | 0.874836 (frontier) |
| 1.8 | 5/6 | 6κ | 0.874353 |
| 1.5 | 5/6 | 6κ | 0.873179 |
| 1.0 | 5/6 | 6κ | 0.870352 |
| 1.8 | 0.80 | 4 | 0.874107 |
| 1.5 | 0.75 | 3 | 0.872384 |
| 1.0 | 0.70 | 2 | 0.869446 |
| 1.0 | 0.50 | 1 | 0.868545 |

The c_I rows ignore the inverse slot supply ℓ/h, which certcI2.py enforces; see the supply-capped re-run below. The LOW-Δ prices are not reachable (§87).
**§88 addendum: an unconditional corollary for a single plain block (PROVED from §88 and the plain fourth moment).** With w_u = |S_u(X)|² and
Σ_u|S_u(X)|⁴ ≪ U^{1+ε} for X ≤ U^{1/2} (plain × plain at total length 1, Lemma 18.1 without slots), §88 gives Σ_u|S_u(X)Q_u|² ≪ U^{1+ε} for every slot polynomial
of length Z ≤ U^{1/4}. Lemma 18.1 with one plain factor allows only m + 6κz ≤ 1 (z ≤ 0.11 at m = 1/2). This does not reach the record, which uses only two-copy legs:
- the plain leg Σ|S²Q|² would need the eighth moment;
- the HIGH numerator with this single-block bound and Hölder exponent 2 gives (R + 1)/2 − q/4 ≈ 0.737, against Xu–Sang's (3R + 1)/4 − ϑ/2 ≈ 0.69.
Similarly, M4 for r ≤ 1/3 (REPORT §14.1) gives Σ|M_u(U^r)Q_u|² ≪ U for Z ≤ U^{1/4}, which is weaker there than Lemma 17.1 (z ≤ (1 − r)/2).
**§89 addendum: slot supply, affine window and floor (MODEL, c61/cp/combo_ws.py, combo_cert.py, combo_cert2.py).** Three more conditions from the certificates
(certcI2.py, sidechecks_rad.py) were missing from the quick model:
- the inverse slot length is capped by the supply ℓ/h;
- Gao's/Collet's affine HIGH capacity C_N at d = h is valid only when the Theorem 2.2 window (3.3) holds there (otherwise C_O is used);
- the floor bin θ > 51/100 + h/3 + 51ℓ/100.
With all three:
- baselines: Gao's constants 0.874838 (h 0.7961); radical constants 0.874841 (h 0.7958), the record 0.8748389 within the model's 10⁻⁵ margin;
- c_I = 1.8 / 1.75 / 1.5 / 1.2 / 1.0 → 0.874353 (supply cap only) / 0.874218 (radical) / 0.873204 / 0.872762 / 0.872762. The same values hold on both constant sets, since
  at h ≈ 0.77–0.81 the window fails and the affine capacity drops out. Below c_I ≈ 1.2 the supply binds;
- LM slope 0.75 → 0.874186 (Gao, window-aware), against 0.874175 without the window check.
**§88 addendum 2: the Möbius fourth moment the inequality needs is flat numerically (NUMERICAL, mmnum/m4.py, existing runs).** At U = 3·10⁵ (1.09·10⁶ rows),
M4/(2·M2²) for M_u(U^r) is 1.051 / 1.070 / 1.081 / 1.090 at r = 0.7 / 0.8 / 0.9 / 1.0. A random completely multiplicative control gives 1.052 / 1.065 / 1.077 / 1.090,
and an independent random model gives 1.00–1.01. The same pattern holds at U = 10⁴, 3·10⁴ and 10⁵. So M4 at r ≈ 0.75, which by §88 would make inverse slots free and by
REPORT §9 gives DH, shows no power loss, only the arithmetic factor of a multiplicative model.
**§89 final table: all conditions enforced (window-aware, slot supply, floor; Gao's constants; c61/cp/combo_cert2.py, wa_*.log).** These supersede the joint table above.

| c_I | LM slope α | plain weight | θ |
|---|---|---|---|
| 2 | 5/6 | 6κ | 0.874838 (frontier 0.874836) |
| 1.8 | 5/6 | 6κ | 0.874355 |
| 1.5 | 5/6 | 6κ | 0.873204 |
| 1.2 (and anything below: the supply binds) | 5/6 | 6κ | 0.872762 |
| 2 | 0.75 | 6κ | 0.874186 |
| 2 | 0.70 | 6κ | 0.873655 |
| 2 | 5/6 | 3 / 2 / 1 | 0.874313 / 0.873382 / 0.870524 |
| 1.5 | 0.75 | 3 | 0.872385 |
| 1.0 | 0.70 | 2 | 0.870775 |

So even three simultaneous critical-length improvements (c_I = 1, slope 0.70, plain weight 2) reach only about −4·10⁻³. The architecture's large moves (to 13/15) come
only from supercritical inputs: full LM, DH for the family, M4 at r ≈ 0.75, or the central sixth moment with λ < 0.20 (up to −7·10⁻³).

## 90. Conditional records at substantial levels, certified in exact rationals (PROVED relative to the radical record's proof; c61/certG.py, sidechecks_rad.py)
Same standard as §73: all three row inequalities G_d(δ, p) ≤ θ − v_d on [0, 3/4] × [0, 1/2], certified by exact-rational box subdivision with the R
upper bound from monotonicity, together with the geometry, LOW-validity, supply and reserve side conditions. New in certG.py:
- the inverse cost c_I and the long-leg slope α are parameters;
- with OLDTOP=1 the top row d = h uses the paper's old capacity min(d, 2y − d)/(6κ). This is needed when the Theorem 2.2 window (3.3) fails at d = h, which
  happens at every geometry below.
With old capacity throughout, the old exponent is affine on [1/2, h] (y < 1/2) and convex below 1/2, so the rows d = h, d₋ and 1/2 suffice. The junction (4.5) is then
not needed. Cycle61's reserves (4.6) and (4.7) need only strict positivity: Cycle61 §4 (4.8)–(4.10) absorbs the κ-feedback by monotonicity, and the d₊ extension costs
5·10⁻¹⁰. The values 1/10 and 1/100 printed there are the record's own margins.
- **LM slope 4/5 ⟹ θ ≤ 87462/100000 = 0.87462** (h = 789/1000; this is RD with c₀ = 1/24).
  - Rows: 113 055 / 20 597 / 31 boxes, min slack 4.7·10⁻¹⁰.
  - Side conditions: only (3.3), which is not used.
  - That is 2.2·10⁻⁴ below the frontier.
- **LM slope 3/4 ⟹ θ ≤ 87425/100000 = 0.87425** (h = 773/1000).
  - Hypothesis: Σ*_u|M_u(U^s)|² ≪ U^{1+(3/4)(s−1)+ε}·U^s for s ∈ [1, 3/2]; copies give 5/6.
  - Rows: d = h certified (18 691 boxes, min slack 1.5·10⁻⁷); d₋ (11 485, 9.8·10⁻⁸); d = 1/2 (35, 1.5·10⁻⁴).
  - Side conditions: all pass except (3.3), which is not used; (4.6) = 0.100005, floor (4.7) = 0.0201.
  - That is 5.9·10⁻⁴ below the frontier 0.874836.
- **Inverse slot cost c_I = 6/5 ⟹ θ ≤ 8728/10000 = 0.8728** (h = 7707/10000).
  - Hypothesis: Lemma 17.1 with r + (6/5)z ≤ m − c₁. The second condition, rescaled as 2r + (6c_I − 4)z ≤ 3m − c₂, does not bind under the supply cap.
  - The slot length is capped by the supply ℓ/d₊ = 0.2277 (rk = 0.727).
  - Rows: d = h certified (39 699 boxes, min slack 1.1·10⁻⁸); d₋ (21 091, 1.2·10⁻⁷); d = 1/2 (45, 4.3·10⁻⁴).
  - Side conditions: (3.3) and (4.5) are not used; (4.6) = 0.0927; floor (4.7) = 0.0164.
  - That is 2.0·10⁻³ below the frontier. Because the supply saturates (§89), any c_I ≤ 6/5 gives the same value from this proof.
- **LM slope 7/10 ⟹ θ ≤ 87368/100000 = 0.87368** (h = 758/1000).
  - Rows: d = h certified (25 107 boxes, min slack 2.1·10⁻⁸); d₋ (11 371, 8.6·10⁻⁸); d = 1/2 (41, 2.9·10⁻³).
  - Side conditions: (3.3) and (4.5) are not used; (4.6) = 0.0925; floor = 0.0233. The geometry has y − x = 1.6·10⁻⁴ > 0, so LOW (5.5) is valid.
  - That is 1.16·10⁻³ below the frontier.
- **LM slope 3/5 ⟹ θ ≤ 87255/100000 = 0.87255** (h = 7648/10000; RD with c₀ = 7/18).
  - Rows: 28 491 / 11 641 / 37 boxes, min slack 2.0·10⁻⁹.
  - Side conditions: (4.6) = 0.0897; floor = 0.0176; y − x = 2·10⁻⁴; κ = 0.7451.
- **LM slope 1/2 ⟹ θ ≤ 87148/100000 = 0.87148** (h = 7713/10000; RD with c₀ = 2/3).
  - Rows: 27 675 / 11 669 / 35 boxes, min slack 8.1·10⁻⁹.
  - Side conditions: (4.6) = 0.087; floor = 0.0122; y − x = 3.6·10⁻⁴; κ = 0.743. The model is 0.871452 (slope 0.6 models at 0.872529).
  - That is 3.4·10⁻³ below the frontier.
- **Inverse slot cost c_I = 3/2 ⟹ θ ≤ 87325/100000 = 0.87325** (h = 8095/10000; supply ℓ/d₊ = 0.2145, rk = 0.678).
  - Rows: d = h certified (226 197 boxes, min slack 1.9·10⁻¹⁰); d₋ (86 521, 5.0·10⁻⁹); d = 1/2 (29, 7.1·10⁻⁴).
  - Side conditions: (3.3) is not used; (4.6) = 0.107; floor = 0.0048. The model optimum is 0.873204.
  - That is 1.59·10⁻³ below the frontier. A run at θ = 0.8734 also passed (62 985 / 43 447 / 29 boxes).
- **c_I = 6/5 together with LM slope 7/10 ⟹ θ ≤ 8716/10000 = 0.8716** (h = 7705/10000).
  - Rows: 48 759 / 21 157 / 43 boxes, min slack 6.0·10⁻⁹.
  - Side conditions: (4.6) = 0.0873; floor = 0.0128; y − x = 2·10⁻⁴; κ = 0.7432.
- **c_I ≤ 6/5, LM slope 1/2 and plain weight 1 together ⟹ θ ≤ 8692/10000 = 0.8692** (h = 7853/10000).
  - Rows: 14 995 / 12 267 / 19 boxes, min slack 1.6·10⁻⁷.
  - Side conditions: (4.6) = 0.081; floor = 0.0006; y − x = 0.001; κ = 0.7384. The model is 0.869112.
  - Here the bin floor 51/100 nearly binds. Going lower needs the floor lowered, which §52–§53 show is a margin.
- **Plain slot weight 2 ⟹ θ ≤ 8734/10000 = 0.8734** (h = 8097/10000; certG2.py, CPW = 2; plain slots capped by the supply).
  - Hypothesis: Lemma 18.1 with n₁ + n₂ + 2z ≤ 1 in place of 6κz, i.e. input 2 at θ′ = 1/2. The HIGH capacity becomes min(d, 2y − d)/2.
  - Rows: d = h certified (148 965 boxes, min slack 8.6·10⁻¹⁰); d₋ (45 905, 5.6·10⁻⁸); d = 1/2 (19, 7.9·10⁻³).
  - Side conditions: (3.3) is not used; (4.6) = 0.108; floor = 0.0052. The model is 0.873382.
  - That is 1.44·10⁻³ below the frontier.
- **All three critical-length inputs together (c_I = 6/5, LM slope 7/10, plain weight 2) ⟹ θ ≤ 8709/10000 = 0.8709** (h = 7753/10000; c61/certG2.py with CPW = 2;
  plain slots also capped by the supply).
  - Rows: d = h certified (28 687 boxes, min slack 1.0·10⁻⁷); d₋ (19 583, 2.5·10⁻⁷); d = 1/2 (37, 3.7·10⁻⁴).
  - Side conditions: (3.3) and (4.5) are not used; (4.6) = 0.0857; floor = 0.0091; κ = 0.7418 ∈ [7/10, 3/4] for Cycle61 (4.8); x = 0.4078 < y = 0.4092.
  - That is 3.9·10⁻³ below the frontier.
  - Tightened: **θ ≤ 0.8708** (63 659 / 28 417 / 39 boxes, min slack 9.0·10⁻⁹). All geometry and LOW conditions hold, with y − x = 2·10⁻⁴; (4.6) = 0.0852; floor = 0.0088.

## 91. Correction to the input-4 price: the slot term ϑ/2 (MODEL, exponent bookkeeping; refines §81)
Per block S_m (m ≤ 1/2) the HIGH numerator has two routes:
- the fourth route (Xu–Sang, with slots): A = (3R+1)/4 − ϑ/2, the same for every m;
- the sixth route: B_m = 5R/6 + (1 + λ_m)/6 − ϑm, with D^{−ϑ} damping and no slots.
The numerator exponent is max_m min(A, B_m), so the sixth route helps only if B_m < A for every block. That means
  λ_m < (1 − R)/2 − 3ϑ(1 − 2m)   for all m ∈ [0, 1/2].
At R ≈ 0.59 and ϑ ≈ 0.05:
- tiny blocks need λ₀ < 0.055 (true: λ = 0);
- m = 1/6 needs < 0.105 (true unconditionally: de Faveri's sieve, 2610.04045, gives λ(m) = 0 for m ≤ 1/6);
- m = 1/3 needs < 0.155 (de Faveri gives 1/6; copies would give 0 but need a three-plain-factor critical lemma);
- m = 1/2 needs < 0.205 (the known bound is 3/8).
With λ = 0 on every block the exponent becomes B₀ = (5R + 1)/6, a gain of (1 − R)/12 − ϑ/2 ≈ 0.009 per unit d, not (1 − R)/12 ≈ 0.034.
So §81's "up to −7·10⁻³" ignored the slot term. A rough corrected value is about −2·10⁻³ for the full λ = 0 statement on all blocks with m ∈ (1/6, 1/2].
Slot amplification inside the sixth route (S_m³ with slots) would recover part of the difference; that is not modelled here.
**§90 sanity check: the certificates are sharp (c61/certG_sanity_*.log).** Below the model optimum the same scripts fail, at the binding class p = 1/2:
- c_I = 6/5 at θ = 0.8726 (model 0.872762) fails at δ ≈ 0.552;
- LM slope 3/4 at θ = 0.8740 (model 0.874186) fails at δ ≈ 0.521.
So the certified values are tight to within the model's 10⁻⁵–10⁻⁴ resolution, not artefacts of a loose check.
**§88 addendum 3: general form and where it could act (PROVED; MODEL).** For any P_u = Σ_{N(n)≤L}b_nψ_u(n) with L⁴ ≤ U and any weights w_u ≥ 0:
  Σ_u w_u|P_u|² ≤ ‖b‖²·(Σ_u w_u + U^ε(U·Σ_u w_u²)^{1/2}),
by the multiplicative large sieve over the characters ψ_•(n)ψ̄_•(n′) of modulus ≤ L² in the row variable. It acts only where the weight's fourth moment is known:
- w = |S_m|² with m ≤ 1/2 (two-plain lemma);
- w = |M(U^r)|² with r ≤ 1/3;
- w = |S_m|⁴ with m ≤ 1/8 (eighth moment from de Faveri's sieve, 2610.04045). This frees slots up to U^{1/4} in the Xu–Sang fourth-moment route for blocks m ≤ 1/8, against
  (1 − 2m)/(6κ) ≈ 0.17.
None of these binds. The HIGH numerator is a max over blocks and is set by the central block m = 1/2, and the detector legs sit at m ≈ 0.42 and r ≈ 0.75. So §88 moves no unconditional value.
**§88 addendum 4: the twist sieve cannot replace the copies amplifier of Lemma 18.1 (PROVED structure).** The amplifier h ↦ hp⁶ works because the rows hp⁶ are exact
copies (G(u, hp⁶) = G(u, h)), which gives a lower bound Σ_p|row(hp⁶)|² = |P|·|row(h)|². For k < 6 the twisted rows carry χ̄_u(p)^k. A large sieve over the
twists (§88) bounds Σ_p|Σ_u b(u)χ̄_u(p)^kG(u, h)|² from above, but the induction needs a lower bound in terms of the untwisted row. Iwaniec-type
amplification with amplifier Σ_p x_pχ_{h′}(p) turns this into §88 with weights w_{h′} = |row(h′)|², and so needs the dual fourth moment Σ_{h′}|row(h′)|⁴. That moment is
supercritical. So c_P = 6κ stays tied to sixth-power replication, as in §80.
**§76 addendum (17 Oct 04:23Z): the RD plateau at a larger scale (NUMERICAL, dualsieve/ds2.py, P = 10⁴: 2574 squarefree composite columns, 128 997 rows).**
Σ_{h≤K}|Σ_n c_nχ̄_n(h)|² / (K‖c‖²):
- LM dual μg̃₆: 0.924 / 0.931 / 0.938 / 0.949 / 0.945 / 0.943 at K = P·(1, 2, 4, 8, 16, 32);
- random control: 0.950 / 0.937 / 0.947 / 0.943 / 0.944 / 0.945;
- unweighted g̃₆: 1.108 / 1.149 / 1.137 / 1.147 / 1.135 / 1.122 (theta bias).
The bias-free sequence stays on the coprimality plateau across the rows-dominant range K ≤ 32D, which is the range RD_{c₀} needs (PLAN_RD.md).
At P = 3·10⁴ (7 796 columns, 386 939 rows), K = P·(1, 2, 4, 8, 16, 32):
- LM dual: 0.941 / 0.937 / 0.939 / 0.943 / 0.945 / 0.945;
- random: 0.942 / 0.944 / 0.946 / 0.944 / 0.945 / 0.945;
- unweighted g̃₆: 1.199 / 1.176 / 1.193 / 1.198 / 1.180 / 1.174.
Across P = 3·10³, 10⁴ and 3·10⁴ the bias-free sequence matches the random control to within 1%, while the unweighted one carries a growing theta excess (1.05 → 1.11 → 1.18).

## 92. Averaging the slot primes before positivity: plain slot cost 5 with no zero-free hypothesis (PROVED reduction, at the level of §83 addendum 2; exponent bookkeeping)
**Statement.** Let P = S₁S₂ be the plain product at total length D = U^{n₁+n₂}, and let Q = Σ_m q_mψ_u(m) be a slot over primes (or squarefree products of slot
primes) of norm ≍ Z = U^z, with Σ|q_m|² ≍ 1. Then
  **Σ_u W(u/U)|P(u)Q(u)|² ≪ U^ε(U + D·Z⁵),** i.e. no loss on the region n₁ + n₂ + 5z ≤ 1, for every κ (no zero-free input).
For comparison: the per-pair route (§83 addendum 2) gives n₁ + n₂ + 6z ≤ 1. Lemma 18.1 gives 6z at κ = 1 and 6κz under β* ≤ (1 + κ)/2 (κ ∈ [7/10, 1] in Cycle61).
**Proof.**
1. Expand: Σ_u|PQ|² = Σ_{m,m′}q_m q̄_{m′}T(m, m′), with T(m, m′) = Σ_u W|P|²ψ_u(m)ψ̄_u(m′). The diagonal m = m′ is ≍ Σ_u|P|² ≪ U^{1+ε} (zero-slot Lemma 18.1).
2. For m ≠ m′, Poisson in u and CRT give T(m, m′) = (U/(DZ))·Σ_j Σ_h V_j(h)·B^{(m),j}_h·conj(B^{(m′),j}_h). Here B^{(m)}_h = g̃6(h, m)·B_{hm⁴}, where
   B_k is the dual of P at row k (the identity G(nm, h) ∝ g̃6(h, m)g̃6(hm⁴, n), exact up to a sixth root of unity: §83 addenda 2 and 5). The finitely many j come from expanding
   the reciprocity units, which depend only on classes mod 36, in characters of that finite group. The joint kernel in N(nm), N(n′m′) is separated by one Mellin integral, at cost U^ε.
   Terms with m | n are plain products again (m divides a plain variable), and the zero-slot lemma bounds them by U + D.
3. Let X(m, m′) be the same expression for all pairs, m = m′ included. Then Σ_{m≠m′}q q̄T = Σ_{all}q q̄X − Σ_m|q_m|²X(m, m).
4. Apply Cauchy–Schwarz over m inside the h-sum, not pair by pair: |Σ_{all}| ≪ (U/(DZ))·(Σ_m|q_m|²)·Σ_m Σ_{h≲H}|B_{hm⁴}|², with H = (DZ)²/U. The same bound,
   times max|q_m|² ≍ 1/Z, covers the subtracted diagonal.
5. **The averaging step.** Σ_m Σ_{h≲H}|B_{hm⁴}|² ≤ max_{h′}#{m : m⁴ | h′}·Σ_{h′≲HZ⁴}|B_{h′}|². The multiplicity is ≤ ω(h′) (or τ(h′) for slot products), so it is ≪ U^ε.
   This is the paper's own amplifier device (OpenAI §18.5: "a fixed output row has only boundedly many representations as hp⁶"), applied to the slot primes.
6. The zero-slot Lemma 18.1 in Poisson-dual form (no restriction on lengths) gives Σ_{h′≲X}|B_{h′}|² ≪ U^ε(D + X)‖c‖² for every X.
7. With X = HZ⁴ = D²Z⁶/U, the off-diagonal is ≪ U^ε·(U/(DZ))·(D + D²Z⁶/U) = U^ε(U/Z + DZ⁵). ∎
**Anatomy.** The 6 of §83 is 2 (conductor) + 4 (the sublattice m⁴). Averaging over m before dropping the sublattice condition recovers one of the Z⁴ factors, because the
sublattices m⁴O for different m overlap boundedly. Bounding each pair separately and then summing |q_m q_{m′}| over Z² pairs loses exactly that factor Z.
**Value.** 5 < 6κ iff κ > 5/6, i.e. β* > 11/12. At the record κ₀ = 2θ − 1 ≈ 0.7497 the paper's 6κ₀ = 4.498 is better, so **this does not move θ by itself**.
It strictly improves the hypothesis-free (κ = 1) case of Lemma 18.1, from 6 to 5. Inserted into the κ-induction (deletion at rate κ per unit width), it would give 5κ₀ ≈ 3.75.
That would need the averaging inside the induction of §18.5–18.8, where the slots are column variables and the 6 comes from the p⁶ amplifier's extracted
terms (18.28)–(18.29). This is not done here; see the price below.
**Price of 5κ₀ (MODEL, c61/cp/combo_cert3.py rad 2 5/6 3.7485; wa_rad_w375.log):** θ = 0.874661 (h = 0.7968), only 1.75·10⁻⁴ below the frontier. Since CPW < Λ = 3.9586, Theorem 2.2's window is replaced by the global plain capacity. The plain leg stops binding and the inverse leg (c_I = 2) takes over, so even a ratio-5 amplifier would be worth little. The substantial levers are c_I and LM, both inside RD (PLAN_RD.md).
**§92 addendum: combining with deletion and with the slot lemma does not give 5κ (PROVED exponent bookkeeping).**
- **Conversion plus deletion.** Delete width d pointwise (cost κd; Lemma 18.1's (18.12)) and convert the remaining slot by §92. The excess is
  (A + 4z − M − 5d)₊ + κd, with A = n₁ + n₂ + z. The optimum d = (A + 4z − M)/5 leaves κ(A + 4z − M)/5 = κ(z − z₁) on the plain leg, with z₁ = (1 − 2m)/5.
  The paper's lemma with deletion gives (A + (6κ − 1)z − M)/6 = κ(z − z₀), with z₀ = (1 − 2m)/(6κ). Both have marginal cost κ > 2q beyond their threshold, and z₁ < z₀ when 6κ < 5.
- **Partial conversion (convert z_a, keep z_b as live slots).** The marked dual sieve at the extended rows X = D²Z_b²Z_a⁶/U holds only off the supercritical
  range (Lemma 18.1 case 2). That forces D·Z_a⁶·Z_b^{6κ} ≤ U, which is worse than both pure routes.
- **Reason.** The averaging gain of §92 lives in the supercritical dual range X > D. Only the zero-slot assertion ("no restriction on lengths") reaches it, so every live slot left
  in place gives the gain back. The paper's 6 is the p⁶ amplifier's ratio (row extension 6ℓ* against pool saving ℓ*, (18.38)). Getting 5κ needs an amplifier of ratio 5,
  and §80 shows why shifts h ↦ hp^k with k < 6 carry column phases.
- **Shift amplifiers of order k < 6 (h ↦ hp^k) checked (HEURISTIC, 06:20Z).** The column phase is χ_n(p)^k. Splitting columns by its value costs only a constant for one p, but the
  amplifier needs the whole pool and must be large at every column: |Σ_{p∈P} a_pχ_n(p)^k| ≫ |P| for each n. Only k ≡ 0 mod 6 guarantees this. For k < 6 the pool sum has size |P|^{1/2}
  on average over n (large sieve for the order-6/gcd(k,6) character in n) and is ≥ |P|^{1−ε} only on a sparse set of columns. So it gives no uniform amplification, and restricting to
  that sparse set loses more than the shorter row extension kℓ* gains. The ratio-6 amplifier is the only one available, consistent with §80.

## 93. The classical-strength ladder: c_I = 7/4, 5/3, 107/66 certified, and what they mean (PROVED relative to the record's proof; analogy labelled)
- **Certificate (c61/certG.py, window-aware top row; certG_cI175b.log).** Lemma 17.1 with r + (7/4)z ≤ m − c₁ ⟹ θ ≤ 87425/100000 at h = 806/1000.
  - Rows: new d = h, 63 519 boxes, min slack 3.3·10⁻⁹; old d₋, 35 621, 4.9·10⁻⁸; old d = 1/2, 29, 3.1·10⁻³.
  - Side conditions (sidechecks_rad.py): the window (3.3) holds, so the affine capacity is legitimate at d = h. The junction (4.5), (4.6) > 1/10, (3.7) and (3.8) pass (S_b(t₀) = 1.8·10⁻⁵), and the case
    conditions pass. (4.7) = 0.0091 > 0, and only positivity is needed. The earlier trial at h = 0.8108 failed (3.8). At θ = 0.87425 the admissible strip is narrow: h = 0.806 and 0.807 pass, 0.805 fails the window (3.3), and 0.808 fails the case condition 41/100 + N + (9/2)zA_max < 1 (scan of sidechecks_rad.py).
  - The model optimum is 0.874219 (cp/fixh.py), so the certificate is within 3·10⁻⁵. Sanity: θ = 0.8742 fails at (δ, p) = (0.4755, 1/2) (certG_sanity_cI175_low.log).
- **Meaning (PROVED exponent identity).** At the binding r = 3/4 the marked dual length is J = N²Z²/U = U^{2r+2z−1}. So c_I = 2 is exactly self-dual (J = N). c_I = 7/4 needs J = N^{1+c₀} with
  c₀ = (1 − r)/(7r) ≈ 0.048, the marked version of RD at about the same strength as the unmarked RD_{1/24} that gives LM slope 4/5.
- **Analogy (labelled).** In the large-sieve form of FOUR_INPUTS input 3, the progressions u mod pp′ have barrier q ≤ U/N, i.e. c_I = 2. A gain of 8/7 in the level exponent, the factor by
  which Bombieri–Friedlander–Iwaniec beat the √x barrier for well-factorable weights, gives exactly c_I = 7/4. So **a BFI-strength theorem for the twisted Möbius mean square in progressions to
  well-factorable slot moduli would put the frontier at 0.87425**, 5.9·10⁻⁴ below Gao's 0.874836. That is about 150 times the last upstream step (0.874840 → 0.874836).
- **Both inputs at classical strength (certG_bfi_a080.log; OLDTOP, since only the window (3.3) fails at h = 0.795).** c_I = 7/4 together with LM slope 4/5 ⟹ **θ ≤ 87405/100000**.
  Rows: 61 909 / 24 501 / 31 boxes, min slack 7.0·10⁻⁹. The float cross-check gives −8.8·10⁻⁵, and the model optimum is 0.874030 (cp/fixh.py). That is 7.9·10⁻⁴ below Gao's frontier,
  from RD at c₀ ≈ 1/24 for the marked and unmarked canonical families.
- **Lemma 17.1's second condition at these points (PROVED arithmetic).** Under the supply cap z ≤ ℓ/d₊ =: L, the largest value of 2r + 8z on the used (r, z) curve is 2 + (8 − 2c_I)L, reached at the kink
  r = 1 − c_I L. For 2″, 2‴ and c_I = 5/3 at h = 0.805 it is 2.947 / 2.965 / 2.989 < 3, so the paper's second condition 2r + 8z ≤ 3m − c₂ holds unchanged and the hypothesis is
  the first condition alone. For c_I = 5/3 with slope 4/5 at h = 0.80 it is 3.0002. There the hypothesis includes the rescaled second condition 2r + (6c_I − 4)z ≤ 3m − c₂
  (maximum 2.572), the convention already used for rows 5 and 7.
- **Maynard strength (CITED: Maynard, arXiv 2006.07088, level x^{3/5−ε} for well-factorable weights, extending BFI's x^{4/7−ε}).** The factor 6/5 over the barrier gives
  2z ≤ (6/5)(1 − r), i.e. c_I = 5/3, which is the marked RD with c₀ = (1 − r)/(5r) ≈ 0.067.
  - Alone (certG_maynard.log, OLDTOP; only the window (3.3) fails at h = 0.805): **θ ≤ 0.8740**. Rows 67 837 / 33 529 / 29, min slack 1.5·10⁻⁸. Float cross-check −1.8·10⁻⁴; model optimum 0.873958.
  - With LM slope 4/5 (certG_maynard_a080.log, OLDTOP; only the window (3.3) fails at h = 0.80): **θ ≤ 0.8738**. Rows 58 957 / 29 005 / 29, min slack 1.5·10⁻⁸. Float cross-check −1.7·10⁻⁴;
    model optimum 0.873760. That is 1.04·10⁻³ below the frontier.
- **Lichtman strength (CITED: Lichtman, arXiv 2309.08522, level x^{66/107} for triply well-factorable weights, the current record).** Factor 132/107 gives c_I = 107/66 ≈ 1.621.
  Our slot moduli are products of many short primes, so they are highly factorable.
  - Alone (certG_lichtman.log, OLDTOP, h = 0.805): **θ ≤ 0.87382**. Rows 81 247 / 38 591 / 29, min slack 1.2·10⁻⁹. Model 0.873779; θ = 0.8737 fails at δ = 0.511.
  - With LM slope 4/5 (certG_lichtman_a080.log, h = 0.80): **θ ≤ 0.87362**. Rows 78 781 / 33 635 / 29, min slack 1.5·10⁻⁸. Model 0.873590; θ = 0.8735 fails at δ = 0.508.
  - Here the paper's second condition with coefficient 8 would reach 3.01–3.02 at the supply kink, so the hypothesis includes the rescaled 2r + (6c_I − 4)z ≤ 3m − c₂.
  - Side conditions: only the window (3.3) fails, so the top row uses the old capacity; floor reserves 0.0081 and 0.0091.
- **Summary of the classical-strength ladder (analogy).** The sextic analogue of each of these classical level-of-distribution theorems would move the frontier as follows:
  - BFI 0.87425;
  - Maynard 0.8740;
  - Lichtman 0.87382;
  - with RD_{1/24} (LM slope 4/5) added: 0.87405 / 0.8738 / 0.87362.
  That is 5.9·10⁻⁴ to 1.2·10⁻³ below 0.874836. A cubic/metaplectic version of any of these theorems was not found in the literature (17 Oct search).
- **c_I = 8/5 (Lichtman's level 5/8 under Selberg's eigenvalue conjecture; MODEL only, cp/fixh.py):** 0.873690 at h = 0.805, only 9·10⁻⁵ below the Lichtman row. Not certified.

## 94. Why the BFI dispersion does not transfer directly: the smooth variable carries the Patterson pole (PROVED reformulation + HEURISTIC obstruction)
- **The bilinear form (PROVED, from §88 and §83 addendum 5).** The marked off-diagonal is (U/(NZ²))·Σ_{j≤J}|Σ_p a_p g̃(j, p)B(jp⁴)|², with J = (NZ)²/U and B(k) = Σ_n c_nχ̄_n(k), c_n = μ(n)g̃₆(n)
  (the RD dual sequence). Expanding B(jp⁴) gives the dual Type II form R(j) = Σ_{n,p}c_n a′_p (n/p)₃ χ̄_{np}(j) of FOUR_INPUTS input 3.
- **What BFI needs.** Dispersion needs a long smooth variable to apply Poisson to after Cauchy–Schwarz. Here the only candidate is the smooth factor m in a Heath-Brown decomposition of μ(n) (Type I: n = dm).
  But then c_n contains g̃₆(m)·(d/m)-cross terms, so the m-sum is a smooth sum of cubic Gauss sums twisted by χ̄_m(j)(m/p)₃.
- **Obstruction (HEURISTIC, consistent with §69 and Dunn–Radziwiłł).** Smooth sums Σ_m g̃(m)χ(m)W(m/M) have the Patterson main term ≍ M^{5/6} (the cubic theta residue). Each Type I piece therefore
  carries a theta main term, far above the square-root scale. Only the full Möbius combination cancels it. In §69's data the smooth sum is Σ_{N(c)≤X}g̃(c) = 0.355·X^{5/6} (theta bias), while |Σμg̃| ≤ 0.63·X^{1/2} on every scale up to 2·10⁵. So a sextic BFI needs, before any Kloosterman step,
  a proof that the theta main terms cancel across the combinatorial decomposition of μ. That is the "Möbius orthogonality to theta" of PLAN_RD step 2. The classical setting has no analogue of
  this: the Type I pieces of the primes have smooth main terms that cancel trivially against x/φ(q).
**§94 addendum: the Type I pieces do carry the theta bias (NUMERICAL, mugauss_t1/typeI.py, X = 6·10⁴, 6 055 split/inert primes, cube check 2·10⁻¹³).**
S_D(X) = Σ_{N(c)≤X} g̃(c)·λ_D(c), with λ_D(c) = Σ_{d|c, N(d)≤D} μ(d) (= Σ_{N(d)≤D} μ(d) Σ_m g̃(dm), the Type I sum with a smooth variable m). Normalised by X^{5/6}:

| X | D = 1 | D = 10 | D = 100 | D = 10³ | D = 10⁴ | |Σμg̃|/X^{1/2} |
|---|---|---|---|---|---|---|
| 60 000 | 0.355 | 0.192 | 0.043 | 0.009 | 0.018 | 0.52 |
| 30 000 | 0.356 | 0.203 | 0.051 | 0.017 | −0.002 | 0.23 |
| 15 000 | 0.353 | 0.192 | 0.023 | −0.015 | −0.008 | 0.41 |
| 7 500 | 0.356 | 0.187 | 0.023 | −0.048 | 0.001 | 1.10 |
| 3 750 | 0.344 | 0.189 | 0.044 | −0.027 | 0.001 | 0.23 |

- Each truncation level has a stable theta coefficient c(D) with S_D ≈ c(D)·X^{5/6}: c(1) = 0.355, c(10) ≈ 0.19, c(100) ≈ 0.03, and it is at noise level from D = 10³.
- Only the full Möbius combination is at the square-root scale.
- So in any Heath-Brown/Vaughan decomposition with Type I level D, the Type I part carries c(D)·X^{5/6}. This must cancel against the Type II part before any dispersion or Kloosterman step can give
  square-root-scale bounds. That is the concrete form of the "Möbius orthogonality to theta" obstruction.
**§94 addendum 2: the theta coefficient is a multiplicative density ρ(p) = 1/(N(p) + 1) (NUMERICAL fit and direct test; presumably derivable from the residues of the twisted Kubota series), and Möbius kills it through the pole of ζ_{Q(ω)} (HEURISTIC).**
- **Fit.** The model c(D)/c(1) = Σ_{N(d)≤D, d squarefree} μ(d)·Π_{p|d} 1/(N(p) + 1) predicts 0.550 / 0.132 / 0.022 at D = 10 / 100 / 10³. The data give 0.541 / 0.113 / 0.014 (X = 6·10⁴; the D = 100, 10³ values are noisy).
  The alternatives N(p)^{−1}, N(p)^{−5/6} and N(p)^{−2/3} fit worse or wrongly: they give 0.46/0.075/0.006, 0.29/−0.36/−0.52 and 0.06/−1.03/−1.41. So the theta residue of Σ_m g̃(dm)N(m)^{−s} at s = 5/6 is, to this
  accuracy, the untwisted residue times Π_{p|d}(N(p) + 1)^{−1}, real and positive. Gauss-sum phases do not enter the local factor.
- **Consequence (HEURISTIC).** With a multiplicative theta density ρ, the theta coefficient of Σ_c a(c)g̃(c) for a(c) = Σ_{d|c}α(d) is c(1)·Σ_d α(d)ρ(d). For a = μ, α = μ∗μ, and the coefficient is
  c(1)·Π_p(1 − 2ρ(p)) ≈ c(1)·Π_p(1 − 2/N(p)) = 0. Equivalently Σ_d μ(d)ρ(d)N(d)^{−w} ≈ 1/ζ_{Q(ω)}(1 + w), which vanishes at w = 0. So **the absence of the Patterson bias for Möbius
  (§69) is the pole of the Dedekind zeta at 1**, the same mechanism as Σμ(n) = o(x).
- **At X = 2·10⁵ (mugauss_t1/typeI2.py, 18 020 primes, cube check 8·10⁻¹³).** c(D)/c(1) against the model Σ_{N(d)≤D}μ(d)Π_{p|d}(N(p)+1)^{−1}:

| D | 10 | 30 | 100 | 300 | 10³ | 3·10³ |
|---|---|---|---|---|---|---|
| data (X = 2·10⁵) | 0.545 | 0.313 | 0.120 | 0.047 | 0.022 | 0.018 |
| model | 0.550 | 0.319 | 0.132 | 0.056 | 0.022 | 0.014 |

- **Direct per-prime test (mugauss_t1/perprime.py, X = 6·10⁴).** The fraction of Σ_{N(c)≤X}g̃(c) (= 0.3551·X^{5/6}) carried by the multiples of p, against 1/(N(p) + 1):

| N(p) | 4 | 7 | 7 | 13 | 13 | 19 | 19 | 25 |
|---|---|---|---|---|---|---|---|---|
| fraction | 0.2004 | 0.1298 | 0.1298 | 0.0732 | 0.0732 | 0.0486 | 0.0486 | 0.0395 |
| 1/(N(p)+1) | 0.2000 | 0.1250 | 0.1250 | 0.0714 | 0.0714 | 0.0500 | 0.0500 | 0.0385 |

  Imaginary parts are ≤ 0.003. The agreement is within the expected X^{−1/3} ≈ 0.025 relative noise.
- **Multiplicativity (mugauss_t1/perpair.py, X = 10⁵).** Single primes at X = 10⁵ give 0.1972 / 0.1254 / 0.0740 / 0.0482 against 0.2000 / 0.1250 / 0.0714 / 0.0500 (N(p) = 4, 7, 13, 19).
  For pairs, multiples of p₁p₂ carry 0.0293 / 0.0141 / 0.0085 / 0.0142 / 0.0088–0.0097 / 0.0047 of the mass at N = (4,7) / (4,13) / (4,19) / (7,7) / (7,13) / (13,13). The product model
  Π(N(p)+1)^{−1} gives 0.0250 / 0.0143 / 0.0100 / 0.0156 / 0.0089 / 0.0051. All deviations are below the noise level ≈ (X/N(p₁p₂))^{1/2}/(0.356·X^{5/6}) ≈ 0.01.
- **Interpretation.** ρ(p) = 1/(N(p) + 1) is exactly the proportion of squarefree c divisible by p. So the theta measure has the local densities of the squarefree counting measure, and the
  density argument shows that the main term vanishes in the same way Σμ(c) = o(X) does. **It gives no rate.** As with the prime number theorem, a power saving needs a pole-free strip for
  Σ_c μ(c)g̃(c)χ(c)N(c)^{−s} left of 5/6, and this series has no Euler product (twisted multiplicativity). The data (|Σμg̃| ≤ 0.63·X^{1/2}) point to square-root cancellation, i.e. GRH-like
  behaviour of this one sequence.
- **What a power-saving "Möbius orthogonality to theta" therefore needs.** The theta main terms of the decomposition cancel through 1/ζ_K(1 + w). A power saving requires the error terms of the
  twisted series Σ_m g̃(dm)N(m)^{−s} uniformly in d, which is the second-variable (row) analysis of the metaplectic double Dirichlet series. That makes PLAN_RD step 2 concrete:
  a uniform-in-d power-saving error term for twisted Kubota series past their theta pole.

## 95. Where RD's first step sits against the literature: the Type III wall of Dunn–Radziwiłł (CITED reading, dl_dr/dr.txt §1.5; PROVED exponent remarks)
- **ρ(p) from Patterson (sketch).** With g̃(pm) = g̃(p)g̃(m)(m/p)₃² and cubic reciprocity, the multiples of p contribute g̃(p)N(p)^{−s}·ψ(s; p), where ψ(s; p) = Σ_{p∤m} g̃(p, m)N(m)^{−s} is the Kubota series
  with numerator p. Its residue at the theta pole is ∝ τ(p) ∝ conj(g̃(p))N(p)^{−1/6} (theta coefficient at a prime). So the multiples of p carry |g̃(p)|²N(p)^{−5/6−1/6}·(1 + 1/N(p))^{−1} = 1/(N(p) + 1)
  of the theta mass, which is §94 addendum 2's density. The derivation is standard in shape (Heath-Brown–Patterson), and the numerics agree.
- **Θ rows are harmless for RD (PROVED exponent count).** In the RD dual sum, the copies of a small row k₀ contribute |B(k₀)|²(K/N(k₀))^{1/6}. With Heath-Brown's B(1) = Σμg̃₆ ≪ D^{5/6+ε} this is
  ≪ D^{5/3}K^{1/6}. That is below the error budget D² iff K ≤ D², which holds throughout the RD range K ≤ D^{1+c₀}. So the theta pole is not the RD obstruction; the generic rows are.
- **The power-saving cancellation is open even pointwise (CITED).** Dunn–Radziwiłł (§1.5) explain that the sums of g̃ over A·B·C with A = B = C = X^{1/3} ("Type III") are self-dual under Voronoi,
  and that the only admissible Cauchy–Schwarz leads to a longer off-diagonal: "bounding these sums with a power-saving presents a real challenge that we do not know how to solve". They avoid it with GRH.
  A power saving below X^{5/6} for Σμ(c)g̃(c)χ(c), which is the pointwise shadow of "Möbius orthogonality to theta", needs exactly those Type III sums. So PLAN_RD step 2 contains a problem that
  the leading experts mark as open, and Heath-Brown's cubic large sieve is sharp for the Type II part at A = X^{1/2}. This is the cleanest statement of why RD, and with it LM and c_I < 2,
  is not reachable by the dispersion/large-sieve toolkit. The self-duality at A = B = C = X^{1/3} is the same self-duality as Theorem B of DH_NOTE at the level of the whole family.
- **DR's corrected cubic sieve and a bootstrap check (CITED (1.13)–(1.14) of dr.txt; PROVED exponent count).**
  - DR subtract the theta correlation, (2π)^{2/3}g̃(a)/(3Γ(2/3)N(a)^{1/6})·Σ_b β_b/N(b)^{1/6}, from Σ_a μ²(a)|Σ_b β_b g̃(b)(b/a)₃|². The remainder is (A/B)·Σ_{0<|h|≤B²/A, h ≠ cube}|Σ_b β_b(b/h)₃|²
    + O((AB)^{2/3−δ}B). Their main term in (1.14) is ∝ Σ α_aβ_bμ²(ab)/N(ab)^{1/6}, the same squarefree-density structure as §94 addendum 2.
  - To beat Heath-Brown, the dual term needs cancellation of Σ_b β_b(b/h)₃ against cubic Hecke characters. DR take it from GRH.
  - A self-improving loop would feed in the tracker's own zero-free region Re s > θ for these characters. That gives |Σ_{N(b)≤B} μ(b)(b/h)₃| ≪ B^{θ+ε} and a dual term ≪ B^{1+2θ+ε} (A ≍ B), which beats
    Heath-Brown's B^{7/3} only if θ < 2/3. At θ ≈ 0.875 the bootstrap gives nothing, so it cannot start from the current frontier.

## 96. A density-estimate route past Heath-Brown's cubic sieve for Möbius sequences (PROVED exponent count, relative to DR (1.13); new direction)
- **Setup.** In DR's corrected sieve (1.13), take β = μ (rough part; smooth parts by the fundamental lemma) and A ≍ B. The theta correction vanishes in the limit by §94 addendum 2, so the size of the left
  side is governed by the dual term (A/B)·Σ_{h≲B}|Σ_{N(b)≤B} μ(b)(b/h)₃|², against Heath-Brown's (AB)^{2/3}B = B^{7/3}.
- **Dual term from zero counts.** Under a zero-free strip Re s > θ for the cubic Hecke L-functions L(s, (·/h)₃), 1/L ≪ (N(h)(|t| + 2))^ε to the right of the strip (Littlewood-type).
  Then |Σ_{N(b)≤B}μ(b)(b/h)₃| ≪ B^{β_h+ε}, where β_h is the largest real part of a zero of L(s, (·/h)₃) up to height B^{O(1)}, and hence
  Σ_h|Σ_b μ(b)(b/h)₃|² ≪ max_σ N_cub(σ; H ≍ B)·B^{2σ+ε}, with N_cub the number of h with β_h ≥ σ.
- **Threshold.** This beats B^{7/3} iff N_cub(σ; H) ≪ H^{7/3 − 2σ − η} for every σ ∈ (2/3, θ]. Above θ there are no zeros, from the tracker's own theorem.
- **Self-consistency check (PROVED).** Density estimates derived from Heath-Brown's cubic sieve by Montgomery's method give N_cub(σ) ≪ H^{3/2−σ} for σ ≤ 5/6 (detector length H^{1/2}) and H^{4−4σ}
  for σ ≥ 5/6 (length H²). Both meet the threshold H^{7/3−2σ} exactly at σ = 5/6 and nowhere else. So the loop reproduces the exponent 2/3 and gives no gain; it can't start without new input.
- **What breaks it.** Any zero-density estimate for the cubic family with N_cub(5/6; H) ≪ H^{2/3−η} (DH would give H^{1/3}), uniform in a short T-aspect, gives a cubic large sieve beyond (AB)^{2/3}
  for μ-supported sequences, via DR (1.13). The OpenAI-type detector beats Montgomery-type counts for the sextic family at the binding (§87: 0.592 against 0.631 at σ = 0.7335).
  **The same architecture run on the cubic family {(·/h)₃} at σ = 5/6 is the natural source of such an estimate.** That would be the first step of a self-improving loop: density (cubic) ⟹ corrected
  cubic sieve for μ ⟹ Type II progress for Σμg̃ ⟹ RD-type input ⟹ better density. Each arrow is a research step; only the exponent bookkeeping above is done.
- **Remark: a threshold at θ = 5/6 (PROVED exponent count).** Heath-Brown-derived densities stay strictly below the threshold for every σ ≠ 5/6. So if the zero-free region for the cubic Hecke family
  reached Re s > θ with θ < 5/6, the maximum over σ ∈ (2/3, θ] would be strictly below B^{7/3}, and the loop would start from Heath-Brown alone. With θ = 5/6 − γ the cubic sieve for μ-sequences improves
  to B^{7/3 − c·γ} (c > 0 from the slopes 1 and 2 of the two density branches). The current frontier 0.8748, and even 13/15, is above 5/6, so this does not fire yet. It marks 5/6 as a
  point where the problem would become self-improving.
- **Where the loop stops (PROVED exponent count).** The sextic detector itself beats the Montgomery-type densities with room to spare at the relevant points. At σ = 5/6 (δ = 2σ − 1 = 2/3)
  its count is R ≈ 0.38 against the cubic threshold 2/3. For the general sextic sieve (de Faveri's A^{5/6}B^{1/3}), the self-consistency point is σ = 2/3, where R(1/3) ≈ 0.72 against the threshold 5/6.
  So density input is not what is missing for the critical (A ≍ B) sieve. But the arrow "critical sieve for μ ⟹ RD" fails: RD is the rows-dominant statement, and its dual is LM at s > 1 (§70).
  The zero-count route to LM gives slope α ≥ (2σ − 1) + E(σ)/(s − 1), with E = R − (2 − 2σ) ≈ 0.06 the detector's excess over DH. Near s = 1.1–1.25 that is > 1, and even with E = 0 it is 2θ − 1 = 3/4.
  So density estimates, however good, cannot give RD (α = 0 just past s = 1). They give at most the slope 2θ − 1, which sits below the copies' 5/6 only once DH-quality counts are available.
  §96 therefore improves the critical cubic sieve for Möbius sequences (a statement of independent interest) but does not move θ by itself.
- **The detector's density curve (MODEL, fastRs at Gao's constants, θ = 0.874836, h = 0.7961; R_max = max_p R(δ, p), σ = (1 + δ)/2).**

| σ | 0.600 | 0.650 | 0.667 | 0.700 | 0.734 | 0.750 | 0.800 | 0.833 | 0.850 | 0.875 |
|---|---|---|---|---|---|---|---|---|---|---|
| R_max | 0.855 | 0.773 | 0.744 | 0.683 | 0.617 | 0.583 | 0.474 | 0.394 | 0.352 | 0.286 |
| excess over DH 2(1 − σ) | 0.055 | 0.073 | 0.077 | 0.083 | 0.084 | 0.083 | 0.074 | 0.061 | 0.052 | 0.036 |
| Montgomery-type (de Faveri sieve) | 0.900 | 0.850 | 0.833 | 0.750 | 0.666 | 0.625 | 0.500 | 0.417 | 0.375 | 0.312 |

  The excess peaks at 0.084 near σ = 0.73, the binding point. max_σ(R_max + 2σ) = 2.084: as a zero-count bound for the critical sextic Möbius moment this is B^{2.084}, worse than the paper's sharp
  critical moment B^{2+ε} (Lemma 17.1 at z = 0). This confirms that the density side and the moment side meet at the fixed point of §75.

## 97. Balanced Type II sums with Möbius structure cancel to square-root size (NUMERICAL, mugauss_t1/typeII.py, X = 10⁵, 9 602 primes)
With w(c) = #{d | c : X^{0.4} < N(d) ≤ X^{0.6}} (balanced factorisations, the range where Heath-Brown's cubic sieve is sharp at X^{5/6} and where Dunn–Radziwiłł's Type II bound stops):

| N(c) ≤ | Σ g̃(c)w(c) / X^{5/6} | |Σ μ(c)g̃(c)w(c)| / X^{1/2} |
|---|---|---|
| 10⁵ | 0.269 | 0.21 |
| 5·10⁴ | 0.259 | 0.46 |
| 2.5·10⁴ | 0.238 | 0.26 |
| 1.25·10⁴ | 0.253 | 0.62 |

The theta-correlated sum sits at the X^{5/6} barrier, as it must (it is the Heath-Brown extremal shape). The Möbius-structured sum Σ_{d∈range} μ(d)Σ_m μ(m)g̃(dm) is at square-root size, about X^{1/3}
below the barrier. So the obstruction to beating X^{5/6} in the balanced range is the theta-correlated shape, and Möbius-type bilinear forms show full cancellation. This is numerical
evidence for the missing Type II/III step of PLAN_RD step 2: the truth is far below the barrier, and what is missing is a proof.
**§96 addendum: literature (CITED, 17 Oct search).** No zero-density estimate for the cubic Hecke family over Q(ω) with an explicit exponent beating the Montgomery/Heath-Brown type was found.
The nearest results are about one-level densities: David–Güloğlu under GRH, with support past (−1, 1) using Heath-Brown's work, and a more restricted unconditional version (as cited in arXiv 2511.08783).
Thorner–Zaman give unconditional log-free density estimates for GL_n families via a new large sieve, but not with the exponent needed here. The input N_cub(5/6; H) ≪ H^{2/3−η} of §96 therefore appears to be new.

## 98. The c_I = 2 barrier is the dual diagonal, not a missing upper bound (PROVED exponent count from the paper's LM-copies bound)
- **The rows-dominant large-sieve upper bound holds for the Möbius dual.** Unnormalised, with B(k) = Σ_{n≤N} μ(n)g̃(k, n): Poisson gives Σ_{k≤K}|B(k)|² = K·N + (N/U′)(Σ_{u≤U′}|M_u(N)|² − U′N), with
  K = N²/U′ and N = U′^s. The paper's copies bound Σ_u|M_u(N)|² ≪ U′^{1+(5/6)(s−1)+ε}N, valid for s ∈ [1, 3/2], gives an excess N²(N/U′)^{5/6} = N²(K/N)^{5/6}.
  This is ≤ KN exactly when K ≥ N. **So Σ_{k≤K}|B(k)|² ≪ U^ε·(N + K)·N for all K ≤ N^{4/3}**, an orthogonal large-sieve upper bound in the rows-dominant range, unconditionally.
  The LM loss (N/U′)^{5/6} is invisible against the dual's large diagonal K·N. In fact the trivial large-sieve shape Σ_u|M_u|² ≪ (U′ + N)N (slope 1) already gives an excess N³/U′ = NK, so the dual upper bound (N + K)N is just the standard large sieve shape seen from the dual side, and copies are not needed for it. RD is hard only because it demands the error O(N²), far below that diagonal.
- **Consequence for the inverse slot.** With the per-slot-prime bound S₁ = max_m Σ_{h≤H}|B_{hm⁴}|² ≪ D + H, the slot-pair route of §92 gives Off ≪ (U/D)(D + H) = U + DZ², i.e. cost 2. This
  holds in the rows-dominant range too, so extending the dual range buys nothing. The averaged (union) route gives U/Z + DZ⁵, i.e. cost 5. Beating 2 needs
  Σ_h B_{hm⁴}·conj(B_{hm′⁴}) = o(diagonal) for m ≠ m′, an orthogonality between distinct sublattice shifts of the Möbius dual. That is the marked RD, and no upper bound of large-sieve type
  can give it, since Cauchy–Schwarz in h erases exactly that cancellation.
- **Same anatomy on the plain side.** For plain products the per-m bound S₁ ≪ D + H is itself missing; the general sextic sieve loses (D²/V)^{1/6} ≈ U^{0.13} on the thin subfamilies v·m², against the
  allowed Z^{0.36} ≈ U^{0.014–0.029}. If it held, cost 2 (row 6, 0.8734) would follow by the same two lines. So input 2 is an upper-bound problem (thin-family large sieve), while input 3 is an
  orthogonality problem (marked RD). These are different kinds of statement, and §92's averaging addresses neither.
- **Input 2 by Hölder over the slot primes (PROVED exponent count; negative).** Σ_mΣ_v|P(vm²)|² = Σ_u|P(u)|²ω_Z(u), with ω_Z(u) = #{m ~ Z prime : m² | u} of mean 1/Z. Hölder with exponents (p, p′) and the plain
  eighth moment (copies loss U^{0.56}) gives U^{1+0.56ε}Z^{−ε} at p = 1 + ε, against the needed U·Z^{−0.64}. This fails for every ε. The additive-character expansion of 1_{m²|u} returns
  the shifted-sublattice dual sums of §83, whose sum over all residues is the identity, so positivity again gives no saving. The thin-family large sieve (§98) is the only form of input 2.

## 99. Bootstrapping the incumbent region into the HIGH numerator: a second threshold at θ₀ ≈ 0.755 (PROVED exponent count)
- **Pointwise bound from the incumbent.** For the AFE block S_m(u) = Σ_{n~U^m}ψ_u(n)n^{−1/2−ϑ−it}, Perron with L(1/2 + ϑ + it + w, ψ_u) and a shift to Re(1/2 + ϑ + w) = θ₀ + ε gives
  |S_m(u)| ≪ U^{m(θ₀ − 1/2 − ϑ) + ε}. Here θ₀ is the incumbent zero-free abscissa, valid for these Hecke characters by the tracker's theorem, and L ≪ U^ε right of it (Littlewood-type).
  At m = 1/2, θ₀ = 0.874836 and ϑ ≈ 0.05: |S_{1/2}|² ≪ U^{0.325}. This beats the Burgess cap U^{3/8} of §82 and is the same mechanism as the paper's κ = 2β* − 1 for prime slots.
- **Value.** The sixth route needs λ = max|S_{1/2}|² exponent < (1 − R)/2 ≈ 0.205, i.e. θ₀ < 1/2 + ϑ + 0.205 ≈ 0.755. Hölder between the fourth route and the pointwise cap helps
  only if max|S_{1/2}| ≪ U^{(1−R)/4} ≈ U^{0.102}, the same threshold. So at the current θ₀ the bootstrap improves input 4's known λ from 3/8 to 0.325 but does not reach 0.205.
- **Structure.** With §96's 5/6 for the cubic-sieve loop, this is a second self-improvement threshold, θ₀ ≈ 0.755 for the HIGH numerator. Both lie far below the frontier, and between
  them and 0.8748 the architecture's own pointwise information is already used, through κ.

## 100. "Marks outside the conductor" (§32's suggestion) is BFI dispersion on the row function |M_u|², and it needs a shifted Möbius fourth moment (PROVED reduction, exponent count)
- §32 noted that c_I = 2 comes from the slot primes sitting inside the theta-reflection conductor, and that removing it needs the marks detected after Poisson. Concretely, keep ψ_u(p)ψ̄_u(p′) as a
  periodic weight in the smooth row variable u (period q = pp′) and expand it in additive characters mod q. Each frequency b/q shifts the Poisson dual window to h ≈ b·nn′/q. That couples n and n′, so the
  dual no longer factors, and this is exactly the situation of the dispersion method.
- The data are: moduli q = pp′, well-factorable products of slot primes with product weights a_pā_{p′}; and the sequence w_u = |M_u(N)|² = Σ_{n,n′}μ(n)μ(n′)χ_{nn′}(u). As a function of u this is
  a sum of characters of conductor N² ≫ U. So it is not a Type I sequence in u (conductor above the length), unlike the Type I pieces of Λ in BFI.
- BFI's step is Cauchy–Schwarz over one modulus factor q₁ = p and opening the square. That produces Σ_{p′}Σ_{u ≡ u′ (p′)}w_u w_{u′}, i.e. shifted correlations Σ_u w_u w_{u+p′k} of the Möbius mean-square
  weights. This is a shifted Möbius fourth moment, the analogue of the additive divisor problem with |M_u|² in place of τ. Unshifted it is M4 at r ≈ 3/4, which is supercritical (§13.2, §33). So the
  dispersion route to c_I < 2 needs M4-type information. That matches §88, where M4 makes slots free outright: the two routes meet at the same input.
- **M4 from zero counts (PROVED exponent count, with §96's density curve).** Rows with |M_u(N)| ≈ N^{σ′} need a zero with real part ≥ σ′ (Perron plus the Littlewood-type bound), so they number ≤ U^{R(2σ′−1)}. Hence
  Σ_u|M_u(U^r)|⁴ ≪ max_{σ′∈[1/2,θ₀]} U^{R(2σ′−1) + 4σ′r}. At r = 0.745 the excess over the target 1 + 2r is λ = 0.15 / 0.22 / 0.31 / 0.37 / 0.40 at σ′ = 0.60 / 0.65 / 0.73 / 0.80 / θ₀.
  So λ_zc ≈ 0.40, the same as copies (0.408) and far from the needed 0.114 (§33). Reaching 0.114 would need θ₀ ≲ 0.575. M4 meets the density side at the same fixed point as LM.
**§98 addendum: the thin-family plain fourth moment is at natural size (NUMERICAL, mmnum/thin.py, U = 10⁶, plain block X = U^{0.35}, 3.6·10⁶ sixth-power-free rows).**
The mean of |S₁S₂(u)|² over rows with m² | u, divided by the mean over all rows:

| N(m) | 7 | 13 | 19 | 25 | 31 | 37 | 43 | 61 | 67 | 73 | 79 | 97 | 103 | 109 | 121 | 127 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| ratio | 0.67 | 0.66 | 1.01 | 1.02 | 1.04 | 1.02 | 1.01 | 1.14 | 0.91 | 1.05 | 0.92 | 0.72 | 0.87 | 1.08 | 1.06 | 1.11 |
| random-subset sd | 0.01 | 0.01 | 0.02 | 0.03 | 0.04 | 0.05 | 0.05 | 0.06 | 0.07 | 0.08 | 0.09 | 0.12 | 0.11 | 0.14 | 0.13 | 0.14 |

The twist rule of Lemma 18.1 corresponds to a ratio ≈ N(m) and positivity to N(m)². The data sit at ratio ≈ 1 (a local factor below 1 for N(m) = 7, 13), i.e. loss Z⁰. The target needs only
loss < Z^{0.36}. So input 2's thin-family large sieve is true with room to spare in the data. What is missing is a proof that does not reflect through the thin family's conductor.
At the binding plain length X = U^{0.417} (same rows): ratio 0.73 / 0.86 / 0.87 / 0.99 / 0.86 / 0.82 / 0.95 / 1.02 / 1.01 / 0.99 / 0.86 / 1.12 / 0.95 at N(m) = 7 / 13 / 19 / 25 / 31 / 37 / 43 / 61 / 67 / 73 / 79 / 97 / 103,
all ≤ 1.12. The twist rule would give ratios of 7 to 103. So at the record's own plain length the thin-family moment shows no loss up to N(m) = U^{0.335}, far beyond the binding slot length
U^{0.04–0.08}. The proof target is therefore a real truth with a large margin (needed: ratio ≤ N(m)^{0.36}).
At U = 3·10⁶ (X = U^{0.417}, 1.09·10⁷ rows, mmnum/thin_3e6.log): ratio 0.77–1.25 for 23 prime norms N(m) = 7 … 199 (N(m) up to U^{0.355}); random-subset controls 0.96–1.02. Fitted slope of log ratio
against log N(m): 0.066 overall, 0.051 for N(m) ≥ 19 (U = 10⁶ gave 0.10 / 0.084). The fitted loss exponent falls as U grows, consistent with loss N(m)^{o(1)}, far below the 0.36 needed.
- **Gauss expansion of the cubic twist (PROVED bookkeeping, 06:20Z; route closed).** Write (n/m)₃ = τ(χ̄)⁻¹Σ_{b mod m}χ̄(b)e(bn/m). Cauchy over b gives |S(vm²)|² ≤ Σ_b|S_b(v)|², where S_b is the
  plain sum with additive twist e(bn/m) over the full family v ≤ V. At natural size per b, this totals N(m)·V·X, which is the twist rule's loss N(m) again. A saving would need cancellation in
  the b-sum weighted by χ̄(b): a correlation of the N(m) additively twisted plain sums across b. That is the same orthogonality between distinct shifts that marked RD needs (§98). So input 2 and
  input 3 share one missing step: decorrelation of the family across residue shifts modulo the slot modulus.
- **Input 2 from a restricted optimal sieve (PROVED reduction, 06:22Z).** Moving χ_m onto the coefficients turns the thin family into the full family over v ≤ V = U/Z² with coefficients
  c_nχ_m(n), where Z = N(m). An optimal sieve (V + X²)·‖c‖² would make it lossless whenever V ≥ X², i.e. Z ≤ U^{(1−0.834)/2} = U^{0.083}. That covers the whole binding slot range U^{0.04–0.08}.
  The known sieve fails because of its middle term (VX²)^{2/3} (Heath-Brown type, the U^{0.13} loss). That term is attained by theta-type (Gauss-sum) coefficients. So input 2 needs the optimal
  sieve only for coefficients with no theta correlation. Our coefficients are plain times a fixed cubic character mod m. Their theta correlation is Σ_n w(n)χ_m(n)g̃(n), whose X^{5/6} term is
  the residue of the χ_m-twisted Kubota series at s = 5/6. Next step: compute that residue, which may be zero for nontrivial cubic χ_m (Patterson's local factors), and turn "no theta
  correlation" into removal of the middle term (the Dunn–Radziwiłł decomposition, §40, §95).
- **The twisted residue, computed (NUMERICAL, mugauss_t1/twist.py, twist_4e4.log, Z[ω], squarefree c, N(c) ≤ 4·10⁴).** Σ g̃(c)/X^{5/6}: untwisted 0.348–0.373 (stable).
  With χ_m (N(m) = 7, 13) the term is 0.03 → 0.002 as X grows from 625 to 4·10⁴, so the residue vanishes. With χ_m² it is 0.19–0.24 (stable), so the residue survives at about 0.63 of the
  untwisted value. So one orientation of the cubic twist is orthogonal to theta and the other is not. In the sieve both c_nχ_m(n) and its conjugate meet g̃, so the thin family's coefficients
  still correlate with theta through χ̄_m = χ_m². Removing the middle term needs a sieve that separates the two orientations (pairing only c·χ_m with g̃), or a correlation bound
  for χ_m²: the surviving residue is a fixed constant, which may still be small enough. Next: express the χ_m² residue through Patterson's local factor at m and price it against the
  allowed Z^{0.36} loss.
- **The χ_m² residue law (NUMERICAL, twist2_1e5.log, N(c) ≤ 10⁵, N(m) = 7, 13, 19, 31; HEURISTIC derivation).** |Σ g̃(c)χ_m²(c)| / X^{5/6} = 0.2230, 0.2154, 0.2111, 0.1940, against
  0.3559·N(m)^{−1/6}·N(m)/(N(m)+1) = 0.2252, 0.2154, 0.2070, 0.1945. The phase is fixed per m (−26°, +75°, −108°, −23°). The χ_m sums are 0.002–0.005 with random phases.
  Derivation: by twisted multiplicativity g̃(cm) = g̃(c)g̃(m)χ_m^{±1}(c), so the χ_m² sum is ḡ̃(m)·Σ_{N(c′)≤XN(m), m|c′} g̃(c′) ≈ (XN(m))^{5/6}/(N(m)+1) = X^{5/6}N(m)^{−1/6}·N(m)/(N(m)+1).
  The other orientation has no such representation, and its residue is 0.
  **Pricing.** The thin coefficients' theta correlation is N(m)^{−1/6} of the maximal one. Even if the middle term scaled with the squared correlation, the saving would be N(m)^{1/3} ≤ U^{0.027}
  at the binding slot length, against the U^{0.13} middle-term loss. The residue route alone therefore cannot close input 2. The middle term must be shown to be smaller than its theta
  projection, which again is the decorrelation step.
- **Duality remark (PROVED).** In dual form the thin family is Σ_{h≤H′}|Σ_n c_nχ̄_n(hm⁴)|², with H′ = D²Z²/U. Since χ̄_n(hm⁴) = χ̄_n(m⁴)χ̄_n(h), this is the full family with twisted coefficients c_nχ̄_n(m⁴), and its
  operator norm equals the untwisted one. So no general-coefficient argument can see the thinness. The general sextic sieve is orthogonal only for H′ ≤ D^{1/2}, i.e. (3/2)(n₁ + n₂) + 2z ≤ 1, impossible at
  n₁ + n₂ = 0.834. The proof must use the plain structure of c_n twisted by the cubic character (·/m)₃, without reflecting through the thin family's conductor VZ. That is the precise open step of input 2.
- **Why it is hard, and its Dirichlet analogue (PROVED count; analogy labelled).**
  - L(s, ψ_{vm²}) = L(s, ψ_v·(·/m)₃) has conductor ≍ V·Z. So the thin-family fourth moment Σ_{v≤V}|L(1/2, ψ_v(·/m)₃)|⁴ has AFE length VZ for L², against family size V: it is past the conductor by
    exactly Z (the full family is exactly critical).
  - The Dirichlet analogue is the fourth moment over a coset {χψ : χ mod q} with ψ fixed mod m. There orthogonality over χ reduces it to the shifted divisor problem Σ d(n)d(n + qh), h ≲ m, which is solvable
    with power savings by the δ-method and Kloosterman sums.
  - The sextic family has no such orthogonality: its harmonic analysis is Poisson in the row variable, which gives Gauss sums of modulus nn′ and the self-dual structure of §83–§84.
  - So input 2 is "Lindelöf on average past the conductor by Z" for plain products in a thin family. It has a solvable Dirichlet analogue, but needs a sextic replacement for the additive divisor problem.
  - **Correct analogue.** The sextic family is "one character per modulus", parametrised by the row u, like the quadratic family {χ_d}, not like a coset of a fixed-modulus group. For quadratic characters
    even the untwisted fourth moment has only a log-saving asymptotic (Shen–Stucky, arXiv 2402.01497, §83). So the thin-family statement is open in its closest analogue, although the sextic data show no loss.
- **Fitted loss exponent (ratio ≈ N(m)^c).** At X = U^{0.35}: c = 0.10 unweighted, 0.23 weighted. At X = U^{0.417}: c = 0.10 / 0.15. Both are driven by the local factor at N(m) = 7, 13 (ratios 0.66–0.86 rising to about 1).
  Excluding N(m) < 19 the fit is c ≈ 0.08 (X = U^{0.417}). The twist rule is c = 1 and the target needs c < 0.36, so the data put input 2 inside its target in every fit.

## 101. A row-dependent κ from zero density: correct in principle, worthless at the record (MODEL; PROVED caps)
- **Idea.** Lemma 18.1's κ = 2β* − 1 comes from the pointwise prime-sum bound under the uniform zero-free region Re s > β*. On rows u whose L(s, ψ_u) has no zero beyond σ₁ < β*, the slot bound
  improves to κ₁ = 2σ₁ − 1. By the detector's density curve (§96), rows with a zero beyond σ₁ number ≤ U^{R(2σ₁−1)}. That is below the plain-leg count U^{0.59} once σ₁ ≳ 0.75, so they can be excluded from
  the count at no cost.
- **Caps.**
  - (i) The lemma is stated for κ ∈ [7/10, 1] (Cycle61's extension), so κ₁ is capped at 0.7, i.e. 6κ₁ = 4.2.
  - (ii) Using it needs the induction run on the good rows only. That set is not a family closed under the lemma's Poisson transforms; the bad rows can only be removed by positivity at the top
    level, where deletion costs κ₁ > 2q and gains nothing.
  - (iii) Even granting it, the price is negligible: plain weight 4.2 gives θ = 0.874831 (cp/fixh.py, h = 0.796), since the inverse leg binds as soon as the plain leg improves (§92: even 3.75 gives only 0.874661).
- So the plain-side κ-mechanism is exhausted at the record. The binding legs are the inverse slot cost (c_I = 2) and LM, and both need the cancellation statements of §98.

## 102. Two more routes to the trilinear cancellation, both closed (PROVED exponent remarks)
- **Sum over the slot primes first.** Σ_{m~Z} a_m g̃₆(h, m)χ̄_n(m)⁴ is a prime sum of cubic Gauss sums twisted by a character of conductor ≍ H·N ≫ Z. Heath-Brown–Patterson / Dunn–Radziwiłł bounds give
  Z^{5/6} plus an error with a positive power of that conductor, which is useless for Z ≤ U^{0.15}. Replacing primes by all integers via a Heath-Brown identity in m gives smooth Kubota sums with the same conductor problem.
- **Poisson in h.** The modulus mm′nn′ ≍ Z²N² gives a dual of length Z²N²/H = U, the original family (self-duality, §98). Cauchy in (h, n) returns the slot pairs, i.e. T(m, m′).
- So any proof must keep all three of h, n, m active at once, which is the configuration Dunn–Radziwiłł call self-dual. This is the exact frontier of the method.

## 103. How Kloosterman sums of modulus mm′ enter: factor the row variable (HEURISTIC structure; the exponent count is the next step)
- §100 noted that w_u = |M_u(N)|² is not a convolution in u, so BFI's dispersion does not start. But ψ_u(n) is multiplicative in the row: ψ_{u₁u₂}(n) = ψ_{u₁}(n)ψ_{u₂}(n). So the smooth row weight can
  be factored by a hyperbola decomposition W(N(u)/U) = Σ_{u₁u₂=u} α(u₁)β(u₂) with U₁U₂ = U. Then
  Σ_{u ≡ b (q)} w_u = Σ_{n,n′}μ(n)μ(n′) Σ_{u₁} α(u₁)ξ_{nn′}(u₁) Σ_{u₂ ≡ b·ū₁ (q)} β(u₂)ξ_{nn′}(u₂),
  where ξ_{nn′}(u) = ψ_u(n)ψ̄_u(n′) is a character of u mod nn′. The progression condition now carries the inverse ū₁, which is the BFI mechanism.
- Poisson in u₂ modulo q·nn′ factors by CRT into (mod nn′): the Gauss sums that build the Möbius dual B(·) of §98; and (mod q = mm′): a genuine Kloosterman-type sum S(h, b·ū₁; mm′)-shape, after the
  b-sum against χ_{mm′}. Cauchy–Schwarz over u₁ and the dual variable then gives a Type II form in Kloosterman sums of composite modulus c = mm′ with p ≍ p′ ≍ Z. That is exactly the setting of
  the 2025 non-abelian amplification bound (arXiv 2511.08445: saving c^{−1/12} for sums of length √c when c is a product of two same-size primes), and of Milićević–Qin–Wu (arXiv 2511.07550).
  The average over the moduli c = mm′ is the BFI/Deshouillers–Iwaniec average.
- **What remains to check (next session).**
  - (a) The lengths: whether the Type II Kloosterman form produced has length about √c in both variables at the binding geometry (r = 0.745, z ≈ 0.146 for c_I = 7/4).
  - (b) Whether the mod nn′ Gauss-sum factor can be carried through the Cauchy step without destroying the gain (it enters as the dual coefficients B(hm⁴), whose large sieve is available, §98).
  - (c) Whether a c^{−1/12} saving per modulus beats the loss from Cauchy–Schwarz and from opening the row factorisation.
  If (a)–(c) work, this would give c_I < 2 unconditionally, the first unconditional move of the inverse slot cost. That is the most concrete unconditional route found this session.
- **First pass at checks (a)–(b) (06:10Z).** After Poisson in u₂, the sum over (n, n′) of μ(n)μ(n′)ξ_{nn′}(u₁)·g̃(h; nn′) recombines into |B(h·u₁*)|², the Möbius dual evaluated at a row twisted by u₁
  (ψ_{u₁}(n)g̃(h, n) = g̃(h·u₁*, n) by twisted multiplicativity). So the resulting form is Σ_{u₁,h} α(u₁)|B(h·u₁*)|²·S(h, b·ū₁; mm′). The Kloosterman sum has the right composite modulus, but the coefficient
  couples u₁ and h through B(h·u₁*). It is not a clean Type II form A(u₁)·C(h), and separating the variables by Cauchy–Schwarz over (n, n′) reintroduces T(m, m′).
  The bridge therefore needs either a version of the amplification that tolerates a coupling weight |B(h·u₁*)|² (a "twisted Type II" with a structured, positive weight, whose large sieve §98 controls),
  or a row factorisation that keeps u₁ out of the Möbius characters, e.g. u₁ restricted to sixth powers, where ψ_{u₁} ≡ 1. In the latter case the coupling disappears but u₁ ranges over U^{1/6} values only. That
  last variant is the cleanest thing to try next: a Type II Kloosterman form with one variable over sixth powers of length (U₁)^{1/6}.
- **The sixth-power variant fails on lengths (PROVED count).** With u₁ = a⁶ (so ψ_{u₁} ≡ 1 and no coupling) the Kloosterman variable a has length U₁^{1/6}. Reaching √c = Z = U^{0.146} needs U₁ ≥ U^{0.876},
  which leaves U₂ ≤ U^{0.124}. Poisson in u₂ against modulus mm′·nn′ ≍ U^{1.78} then has dual length ≈ U^{1.66}, far past anything controllable. This variant is just the copies amplifier again (rows a⁶u₂).
  So the row-factorisation bridge stands or falls with a coupled ("twisted") Type II Kloosterman bound carrying the weight |B(h·u₁*)|². That is the precise new estimate to aim for, combining
  non-abelian amplification in mm′ with the large sieve for the Möbius dual in h·u₁*.
- **Without the row factorisation, the exact form of the barrier (PROVED count).** Poisson in u for T(m, m′) = Σ_u w_uχ_q(u), q = mm′, gives (U/(qN²))·Σ_{h≤K} G_q(h)·N·|B(h)|², with K = qN²/U and |G_q(h)| = √q,
  where G_q(h) = χ̄_q(h)·G_q(1).
  - With the large sieve (Σ_{h≤K}|B(h)|² ≪ KN), T ≪ √q·N². That meets the needed UN/Z exactly when r + 2z ≤ 1, which is c_I = 2.
  - If instead |B(h)|² were equidistributed in classes mod q, i.e. Σ_{h≤K} χ̄_q(h)|B(h)|² ≪ KN·q^{−1/2+ε}, then T ≪ N² and the condition becomes r + z ≤ 1, cost 1.
  - So the whole inverse slot cost is the square-root equidistribution of the Möbius dual weights |B(h)|² in residue classes to the small moduli mm′ ≤ U^{2z}, for h up to K = qN²/U.
    This is a Bombieri–Vinogradov-type statement for |B(h)|², and the row factorisation of §103 is the attempt to prove it BFI-style.
  - **BV-type route closes the circle (PROVED count).** By the classical large sieve in h over characters mod q ≤ Q = Z², Σ_{q≤Q}Σ_χ|Σ_{h≤K}χ(h)|B(h)|²|² ≤ (K + Q²)·Σ_{h≤K}|B(h)|⁴.
    If Q² ≤ K (true for z ≤ (2r − 1)/2 ≈ 0.245) and the dual fourth moment is natural, Σ_{h≤K}|B(h)|⁴ ≪ KN²U^ε, then the needed equidistribution holds on average over moduli: cost 1 on average.
    The dual fourth moment is past its conductor (|B|² has length N² > K), the dual form of the Möbius fourth moment. So §88 (M4 ⟹ free slots), §100 (dispersion ⟹ shifted M4) and this BV form all meet at
    one input: the Möbius fourth moment at r ≈ 3/4, needed loss < 0.114, best known ≈ 0.40.
- **A caveat on the classical-strength reading of §93 (PROVED structure; corrects the strength of the analogy).**
  - The row factorisation is an identity: summing the residue b against χ_{mm′} collapses the Kloosterman sums back, Σ_b χ(b)S(h, bū₁; q) = G(χ)²χ(u₁)χ̄(h), and χ(u₁)χ̄(h) = χ̄(h·u₁⁵) since χ⁶ = 1.
    So after the b-sum one is back at Σ_k χ̄(k)|B(k)|². Kloosterman structure survives only if a fixed residue class is kept and Cauchy–Schwarz is applied before the b-sum.
  - BFI, Maynard and Lichtman go beyond √x **for a fixed residue class**. For averages over all characters of each modulus, the multiplicative large sieve is sharp, and no beyond-barrier result is possible in general.
  - Our slot structure |Σ_m a_mψ_u(m)|² is a character twist, with one specific sextic character u ↦ ψ_u(m) per prime m. It is a thin set of characters, so beating the barrier is not excluded,
    but it is not the BFI situation.
  - The rows 2″, 2ᵐ, 2ˡ (c_I = 7/4, 5/3, 107/66) therefore read as "what an improvement of the barrier by the classical factors 8/7, 6/5, 132/107 would give". They are not consequences of a sextic transcription of
    the fixed-residue theorems. The certified θ values are unaffected; only the analogy that motivates them is weaker than §93 suggested.
- **M4 is LM for μ∗μ (PROVED reformulation).** Σ_u|M_u(U^r)|⁴ = Σ_u|Σ_{n≤U^{2r}}(μ∗μ)(n)ψ_u(n)|²: the long mean square at s = 2r ≈ 1.49 for the coefficients μ∗μ. Copies give the loss
  (5/6)(2r − 1) = 0.408 (§33). The needed λ < 0.114 is an LM slope ≲ 0.23 at s = 1.49 for μ∗μ, far beyond slope 4/5 for μ. So the hub input of §88/§100/§103 is the strongest form of the LM problem, and the
  whole wall is LM-type mean squares past the conductor plus their self-duality.

**§98 addendum: absorbing the twist into the Gauss frequency (PROVED bookkeeping against the paper's §18.4–18.6, 06:32Z; route closed).**
- Poisson over thin rows k ∈ m²O ∩ (k₀ + abO) gives Gauss sums G(a, h·m̄²) = χ_a(m²)G(a, h) on a dual range Z^{2z} longer. Since χ_a(m²) = χ̄_a(m⁴), this equals G(a, hm⁴): the untwisted
  Gauss-row norm N_C restricted to frequencies on the sublattice m⁴O, over a range Z^{6z}·T₁. Positivity extends this to the full untwisted norm over that range, with the moving support q unchanged.
- The paper's second transform (18.30)–(18.31) bounds N_C by its diagonal j = 0 term, of exponent K + g − ℓ − s₀, against the allowance a₀ + q_e + …. Its excess is (A − M) − d − (K₀ − K) + …,
  so enlarging the range by 6z requires A + 6z ≤ M. That is cost 6, exactly Lemma 18.1's p⁶ amplifier ratio. The frequency shift h ↦ hm⁴ plus sparser rows is the same row extension.
- Keeping χ_a(m²) in τ instead (range +2z, support q_e + z) is the twist rule, with loss Z^z. So neither presentation beats cost 6 / loss N(m). The thin family needs off-diagonal cancellation
  in the second transform's j-sum, which is not just the diagonal. This confirms §92's diagnosis from inside the paper's proof (oai5.txt lines 9773–10200).
- **Why the loop does not close (PROVED, same bookkeeping).** Evaluating the sublattice norm directly, by Poisson in h′ over m⁴O modulo uv, gives a dual j-range Z^{4z} longer and the
  correlation F(u, v; j·m̄⁴) = χ_u(m⁴)χ̄_v(m⁴)F(u, v; j). The twist comes back as a character on the columns. Each Poisson step moves the slot modulus between "sublattice of rows" and "character
  on columns", and positivity in either form charges it in full. This is the self-duality of §83 and §98 at the level of the paper's own second transform.

**§98 addendum: the de Faveri–Dunn–Hoffstein second moment (CITED: arXiv 2607.07911, read 06:37Z; lead closed).**
- They prove unconditionally that the cubic large sieve is optimal: Ξ₃(A, B) ≫ (AB)^{−ε}(A + B + (AB)^{2/3}). Their quartic bound is A + B + A^{3/4}B^{1/2} + A^{1/2}B^{3/4}, and they conjecture
  A + B + A^{1−1/n}B^{2/n} + A^{2/n}B^{1−1/n} for order n. So no general-coefficient sextic sieve can drop the middle term, which confirms §98's duality remark.
- Their upper-bound input (2.9) is Lindelöf on average over the numerator, Σ_{N(a)≤A}|ψ⁽³⁾(a, 1+ε+it)|² ≪ A^{1+ε}. It comes from Heath-Brown's cubic sieve at the balanced point A against
  AFE length A^{1/2}, where the middle term (A·A^{1/2})^{2/3} equals A exactly. The thin family's dual is the sextic Kubota series at numerators hm⁴, i.e. the χ_m-twisted series. Its
  numerator conductor is larger by N(m)⁴ (m⁴ is sixth-power-free), so the same argument loses a power of N(m), and positivity on the sparse numerators hm⁴ loses too. The method gives
  nothing for thin subfamilies.
- Also seen: Liu, arXiv 2610.12234, varies the paper's (b, ℓ) and gets 0.874957. That is weaker than the tracker's 0.874836, so it is not an input.
- **Dunn, arXiv 2403.13151 (CITED, read 06:39Z; lead closed).** This is the cubic Duke–Iwaniec theorem: it gives a power saving for Σ_{ϖ≤X} ρ_f(ϖ) over primes in progressions mod v, for a fixed cubic
  metaplectic cusp form f. Its uniformity in v is N(v)^4 (Thm 1.2, Thm 1.5). Its Type I estimate (Prop. 1.6) still goes through Heath-Brown's cubic sieve. It covers cuspidal coefficients only.
  Our thin family needs theta (Eisenstein-residue) coefficients twisted at level m, with loss below N(m)^{0.36}. The N(v)^4 uniformity is far too weak, and the theta case is
  Dunn–Radziwiłł, which needs GRH for primes. Nothing transfers.

**Family change upward, to order 12 over Q(ζ₁₂) (CITED reading of oai5.txt pp. 4–6, 06:42Z; closed).** The Gauss-sum bias X^{1/2+1/n} falls with n, so a higher-order family would have weaker
theta barriers. But the architecture is not order-generic. Its base sum averages the Fourier coefficients of Patterson's cubic theta function against the sextic characters, and the reflection
(Prop. 5.1) uses their explicit evaluation as cubic Gauss sums. For n ≥ 4 the theta coefficients are unknown (de Faveri–Dunn–Hoffstein 2607.07911 §1.2 calls this "a central open problem").
The only extra input there is the Suzuki-type evaluation at special indices for n = 4, and the quartic port was stuck at Q^{3/2} (§26–§30). Also, the large sieve's middle term for order n is
conjectured to be A^{1−1/n}B^{2/n}, and the proved bound (AB)^{2/3} does not improve with n. So an order-12 port needs both an unknown theta evaluation and an unproved sieve.

**Where Lemma 18.1's q-loss sits, for row-shift twists (PROVED tracing of oai5.txt (18.21)–(18.31), 06:45Z).**
- For τ(n) = χ_n(m²), the twist is a row shift at every level of the induction. After the first Poisson, χ_a(m²)G(a, h) = G(a, hm⁴). After the second, Lemma 13.4's column characters become
  χ_a(j)χ̄_a(m⁴) = χ_a(jm²). So every child family is again {ψ_{jm²}}: rows on the sublattice m²O.
- The diagonal of the second transform (18.31) does not use q. The allowance's q_e = q + R + E is spare room there. q is needed in the width M = m + q only because the terminal Hecke
  reflection sees the conductor of n ↦ χ_n(jm²), which is N(j)·N(m) (the cubic character mod m has conductor m).
- The lemma's q-loss is not sharp here: the numerics (§98 addenda, U = 3·10⁶) show the thin family at natural size. The dichotomy that the proof cannot escape:
  (a) put τ into the coefficients, a_nτ(n). They are no longer smooth, so only the general sieve applies, and its middle term is sharp (de Faveri–Dunn–Hoffstein);
  (b) put τ into the conductor. Reflection, i.e. Poisson in n, must run modulo m, which lengthens the reflected sums by N(m). That is the twist rule.
  (c) The hybrid splits into progressions n ≡ b mod m, where the coefficients are smooth. Cauchy over the N(m) classes costs N(m), and avoiding it needs Σ_k S_b(k)S̄_{b′}(k) = o(diagonal) for b ≠ b′.
  The diagonal of that correlation vanishes (n ≡ b, n′ ≡ b′, nn′⁵ a sixth power forces b = b′ up to sixth powers). Its off-diagonal is the progression-restricted Gauss-row norm, which needs Poisson mod m: back to (b).
- So the missing estimate is concrete: the off-diagonal of Σ_k S_b(k)S̄_{b′}(k) for b ≠ b′ mod m, i.e. decorrelation of progression-restricted plain sums. This is the same
  statement as §98's Gauss-expansion addendum, now located at a specific place in the paper's induction.

## 104. de Faveri's optimal large sieve (CITED: arXiv 2610.04045, 2 Oct 2026; read 06:49Z) and what it does for the record
- **Theorem 1.1 (PROVED there).** For every n ≥ 3, over n-th-power-free moduli: Θ_n(A, B) ≪ (AB)^ε(A + B + A^{1−1/n}B^{2/n} + A^{2/n}B^{1−1/n}). For n = 6 the middle term is
  A^{5/6}B^{1/3} + A^{1/3}B^{5/6}, replacing BGL's (AB)^{2/3}. This is sharp for n = 3, 4 (de Faveri–Dunn–Hoffstein).
- **Effect on the record (PROVED exponent counts; no gain).** The lossless range is B ≤ A^{1/2}, the same as BGL. Beyond it the loss halves, to (B²/A)^{1/6}.
  - Inverse leg at r = 0.745, general sieve on M·Q: count 0.7345 − 0.133z, best 0.70 at s = 1. Lemma 17.1 gives 0.5935.
  - Pure prime witnesses: count 1 + (2s−1)/6 − 2qs, best 0.70 at s = 1.
  - Plain leg: loss U^{0.111 + z/3}. Lemma 18.1 is lossless in its range.
  - Thin family (input 2): loss (D²/V)^{1/6} = U^{0.111+z/3}, against the allowed U^{0.36z}. The figure U^{0.13} in §98 was already this shape.
  So the structured lemmas dominate the optimal general sieve everywhere the detector uses moments. Theorem A's closure stays exact.
- **Corollary 1.5 (zero density for n-th order families; for n = 6 on [2/3, 11/12]).** δ₆(σ) = (2σ−1)((2σ−1)·36 − (6σ−2)·6 + 3)/((2σ−1)(3−2σ)·36 − 24σ + 2). At the binding σ ≈ 0.733:
  1 − δ₆ = 0.752, against the record's R = 0.5935 (A = 2.82 against 2.22). At σ = 7/8 it is 0.466. No gain.
- **Theorem 1.2 (cubic, PROVED there).** Σ_{N(a)≤A, N(b)≤B} μ²(ab)|L(½, ψχ_{ab²})|² ≪_ψ (AB)^{1+ε}: the cubic second moment ordered by conductor, tight. With A = V, b = slot prime m,
  this is the cubic analogue of input 2 averaged over slot primes: Σ_mΣ_v ≪ (VZ)^{1+ε}, natural size. Our family is sextic: χ_n(vm²) = χ_n(v)·(n/m)₃, a sextic row times a fixed
  cubic twist. de Faveri lists the n ≥ 4 analogue as open ("It would be interesting to obtain analogues ... for n ≥ 4"). It uses structural features of n = 3, notably that the cubic
  conductor of χ_{ab²} is ab.
- **Next.** Decide whether the averaged form (sum over m ~ Z) of input 2 suffices for the plain slot cost. If it does, the target becomes "de Faveri Theorem 1.2 for the sextic family, in the
  short-sum form": a precise and newly reachable-looking statement.
- **The averaged form of input 2 suffices (PROVED from §92's own steps, 06:53Z).** §92 step 4 already applies Cauchy over m inside the h-sum. So the off-diagonal is
  ≪ (U/(DZ))·Σ_{m~Z}Σ_{N(h)≲H}|B_{hm⁴}|². §92 then drops the sublattices by positivity (step 5), which gives cost 5. If instead Σ_mΣ_h|B_{hm⁴}|² ≪ U^ε·Z^{1+c}·(D + H),
  the off-diagonal is UZ^c + DZ^{2+c}. At c = 0 that is cost 2 (row 6, θ ≤ 0.8734 under the certificate). No per-prime bound is needed, only the average over slot primes.
- **Why de Faveri's cubic proof works and the sextic one does not: a parity obstruction (PROVED).** Prop. 8.1 for n = 3 opens the square and Poissons in the squared variable b.
  Its character χ_q(b)² = χ̄_q(b) becomes χ_q(e) on the dual side, the same power as the a-variable, so a and e merge into the single row ae with divisor-bounded multiplicity, and
  Theorem 1.1 finishes. In the sextic case the thin twist is χ_n(m)^{2i} (from G(p^iu, h) = R(p,u)^iχ_p(u)^{2i}G(p^i,h)G(u,h)). It is always cubic, with an even exponent. Poisson in m
  maps exponent 2 to −2 ≡ 4 and back, and Poisson in the row variable maps 1 to 5. A merge needs equal exponents mod 6, i.e. an odd exponent on m, which a cubic twist never has.
  This is the self-duality hm⁴ ↔ vm² of §98, now explained: the cubic case escapes it because −2 ≡ 1 (mod 3).
- **The cubic × quadratic split (PROVED identity; no gain).** Write χ = ξ·ε with ξ = χ⁴ cubic and ε = χ³ quadratic. Then χ_n(hm⁴) = ξ_n(hm)·ε_n(h): the cubic part merges (w = hm) and the
  quadratic part keeps the row h = w/m, so rows hm⁴ = w·m³ with w ∈ mO. Dropping m | w costs Z (density). For each m the family over w carries the quadratic twist
  ε_n(m) = χ_n(m)³ of conductor m, and the twist rule charges another Z, on the D-term as well (UZ). That is worse than §92. A quadratic twist has no Gauss bias, so a twist-free treatment
  of quadratic twists of the sextic family (a Heath-Brown-type quadratic argument inside the sextic family) is the one structural opening this leaves.
- **Where reflection enters Lemma 18.1 (CITED reading of oai5.txt (18.9)–(18.16), 06:55Z).** Reflection, the only step that genuinely needs the twist's conductor in the width
  M = m + q, is used in exactly two places: (i) the padded zero-slot core, where each factor longer than M/2 is reflected once (18.10); (ii) the range 5M/6 < A ≤ M, where the comparison
  with lengths L = M/4 and A − L has its longer factor reflected (18.15). The range A ≤ 5M/6 uses the two-transform estimate, whose diagonal (18.31) does not need q.
  For the thin family at the record, A = 0.834 > 5(1 − z)/6 for every z ≥ 0 (the record's plain length 0.417 sits just past the 5/6 edge). So top-level reflection, and with it the
  N(m) conductor, is unavoidable at the binding lengths. At A ≤ 5(1 − 2z)/6 the top level avoids reflection, but the children of the second transform carry the support q + R + E and
  reflect when their own A_c > 5M_c/6. A version of Lemma 18.1 that keeps row-shift supports out of the width until a child actually reflects is the precise modification to try. Its value
  depends on how often children reflect; that is a full re-derivation of §18.6–18.8 (about 900 lines), not done here.
- **Conjugate grouping (P·Q̄)(P̄·Q) (PROVED; closed, 06:58Z).** Since |PQ|² = |PQ̄|², one can try putting the slot on the C side with the opposite conjugation. Within one side, χ_n(m′)χ̄_{m′}(n) = R(n, m′) is a
  fixed phase, so the C-side Gauss factor would be G(n, h)Ḡ_{m′}(h) with no thin twist. With the trivial bound |Q*(h)|² ≤ Z this would give U + DZ², i.e. cost 2. But the complete
  sum mod L = nmn′m′ carries all cross factors ξ_ℓ(L/ℓ). By reciprocity they contain χ_m(n)²·χ̄_{m′}(n′)², which couples each plain column with the slot of the same conjugation. In the
  conjugate grouping that coupling runs across the C/D split (n ∈ C, m ∈ D). It is a genuine sextic symbol, not a finite ray phase, so Cauchy–Schwarz in h cannot separate it. The thin twist
  is intrinsic to the complete sum, and only its placement depends on the grouping. Closed.
- **Recursion model for "support out of the width" (HEURISTIC exponent model of (18.20)–(18.31), generic c = d = R = E = 0).** One two-transform step maps (m, q, A) to
  (m − σ, q + extracted, A): the second transform's rows j have range Z^{2a₀−K−g} = Z^{m−σ}, while the residual columns keep length about A. So A_k/M_k rises along every branch until
  A_k > 5M_k/6, and then the child reflects with the full conductor. Its bound carries Z^z in both N_C and N_D, which returns ½(z + z) = z to the parent's ledger. So the modification
  gives back the twist rule. Closed at heuristic level.
- **Where the paper uses large sieves, against de Faveri (CITED reading of oai5.txt, 06:59Z).** The sextic (BGL) sieve, Lemma 9.1 with (KD)^{2/3}, is used only in Part I's Prop. 9.2,
  whose count R(δ) = min(1, max(1 − δ/2, 4/3 − δ)) gives 11/12. Part II's row counts come from Prop. 19.2, i.e. Lemmas 17.1 and 18.1. Their terminal energy (§14, Lemma 14.2)
  uses only Goldmakher–Louvel's quadratic sieve, which is already optimal, and Heath-Brown's cubic sieve (14.6). For n = 3, de Faveri's bound A + B + A^{2/3}B^{2/3} + A^{2/3}B^{1/3} has the
  same middle term (AB)^{2/3}, which is sharp (de Faveri–Dunn–Hoffstein). So Theorem 1.1 replaces nothing in the proof of the 7/8-type record. Part I's count would improve to
  max(1 − δ/2, …) with the A^{5/6}B^{1/3} shape, but Part I is not binding.
- **Other recent preprints on the 7/8 architecture (CITED, read 07:01Z).** Liu (arXiv 2610.12234, latest HTML "Improved zero-free half-planes") varies b, ℓ and κ, extends the unweighted fourth
  moment to κ < 3/4, and gets B_new = 0.874957019… (root of 1188ξ³ − 2160ξ² + 171ξ + 190). Sang (arXiv 2610.15291) reaches 0.8748867 through an amplified first moment of the
  reflected numerator. Both are above the tracker's 0.874836, and the κ < 3/4 extension is already in the tracker's pipeline through Cycle61 (κ ∈ [7/10, 1]). No new input to combine.
- **Poisson in the long variable, de Faveri's way, on the sextic averaged thin family (PROVED exponent counts, 07:04Z).** In T = Σ_{m~Z}Σ_{v≤V}|Σ_n c_nχ_n(vm²)|², the long variable is
  v (V = U/Z² ≫ Z), with exponent 1. Poisson in v gives dual rows e ≤ E = N²/V with χ̄_n(e), so the family becomes rows em⁴ and T_off ≈ (V/N)·Σ_{e≤E}Σ_m|Σ_n c″_nχ̄_n(em⁴)|².
  Natural size for that dual family (EZ rows, N columns) would give (V/N)(EZ + N) = NZ + V ≤ VZ, i.e. T natural. So the sextic averaged thin family is equivalent to the dual thin family
  {em⁴} being column-orthogonal (EZ < N for z < 0.055). de Faveri's Θ₆ over all rows up to EZ⁴ loses max(Z⁵, U^{0.084}Z²). The Pólya–Vinogradov Gram bound needs E²Z⁵ ≤ N, which fails.
  Poisson in the short variable m (F = N²/Z) is worse: even the natural bound gives T_off ≈ VN ≫ VZ.
- **Cubic versus sextic architecture (PROVED structural remark).** Write χ = ξ·ε (ξ = χ⁴ cubic, ε = χ³ quadratic). Then χ_n(vm²) = ξ_n(v m̄)·ε_n(v): the cubic component is exactly de
  Faveri's Prop. 8.1 shape (rows a·b², with χ_q(b²) = χ̄_q(b) merging after Poisson), and in a purely cubic family input 2 would follow from his argument. The obstruction is the
  quadratic component ε_n(v), carried by the row variable alone. The paper itself splits the sextic family this way in Lemma 14.2 (χ_k(nb)³ quadratic on rows, χ_P(n)^{−2} cubic on slots),
  bounding the parts with separate quadratic and cubic sieves. The sextic family is forced by the architecture: cubic Gauss sums from the theta coefficients, quadratic from reflection.
  So the missing input is concrete: a de Faveri Prop. 8.1-type bound for rows a·b² in the mixed family χ = (cubic)·(quadratic), i.e. "Prop. 8.1 with a quadratic rider on a".
- **Prop. 8.1 with a quadratic rider, first pass (PROVED identities; open).** For rows a·b² in the mixed family, Poisson in b merges the cubic parts, χ_q(ab²) ↦ ξ_q(ae)·ε_q(a), i.e. rows a·e⁴.
  Writing w = ae gives χ_q(w)ε_q(e) = χ_q(a)ξ_q(e), the same type again. Opening the square gives [Σ_a χ_qχ̄_{q′}(a)]·[Σ_e ξ_qξ̄_{q′}(e)], and double Poisson returns the
  conjugate family χ̄(a*e*⁴). So the quadratic rider regenerates under every Poisson step: the mixed family is self-dual as a whole, while its cubic sub-family is not. Any proof has to
  handle ε_q(a) by something other than Poisson, for instance by Heath-Brown's quadratic recursion run jointly with the cubic merge. No such argument exists yet; this is the frontier
  statement for input 2.
- **Correction to the cubic analogy (PROVED exponent count, 07:06Z).** De Faveri's Prop. 8.1 is a general-coefficient bound, AB + M^{2/3}N^{1/3}Q + M^{1/3}N^{2/3}Q^{2/3}, and his Theorem 1.2
  uses it only at column length Q = √(AB) (the AFE length), where it is lossless (term-2 exponent −0.13 to −0.16 below AB). Our thin family has rows A·B = V·Z = U^{1−z} and columns
  Q = U^{0.834}, far beyond √(AB) ≈ U^{0.48}. There even the cubic Prop. 8.1 loses U^{0.19–0.25} (z = 0.02–0.08). So the cubic case is not solved by his methods in our range either. The
  long-column regime is handled only by smoothness-based recursions (Lemma 18.1), and there the obstruction is the slot conductor at reflection (above). The parity remark stands for the
  merge step, but the "cubic input 2 would follow" sentence holds only in the short-column range Q ≤ √(AB).
- **CORRECTED (07:10Z): general coefficients do obstruct input 2 (PROVED lower bound via the Poisson identity).** An earlier version of this bullet averaged |ρ(r)|² over all rows
  and concluded there was no obstruction. That was wrong: the obstruction concentrates on sixth-power copy rows of the dual. Fix one m, so x_n = λ_nχ_n(m)². By §98's identity,
  Σ_v|Σ_n x_nχ_n(v)|² = V‖x‖² + (V/N)(Σ_{e≤W}|Σ_n x_n g̃(n)χ̄_n(e)|² − W‖x‖²), with W = N²/V. Take the dual vector aligned with one row, λ_n g̃(n) = χ_n(m₀⁴), which is legitimate since
  |g̃(n)| = 1 on squarefree n. Then the W^{1/6} copy rows e = f⁶ each give N², and this exceeds the dual diagonal WN. So for m = m₀ the original family exceeds its natural size VN by W^{1/6}
  (per-m loss U^{0.118–0.138} for z = 0.02–0.08), and the average over m ~ Z exceeds VZN by W^{1/6}/Z = U^{0.098–0.058}. So no general-coefficient theorem can give input 2
  in our range; this is de Faveri's A^{5/6}B^{1/3} loss, reduced only by the factor Z from averaging. Any proof must use the plain (smooth) structure of the coefficients, and the thin-family numerics
  (natural size) are a statement about plain coefficients only.
- **Gao's small-width rate is tight within his inequality system (PROVED, exact at ν = q = 839/1250, 07:08Z).** ρ(q) = 0.38831, 4ρ − 1 = 0.55322 > (2 + 22ρ)/21 = 0.50203, so
  γ = 0.11798 comes from the c ≤ d branch. Adding the coupling 4L_c ≤ 5(c + d) to that branch does not help: the worst case is R = c with d ≫ c, where min(5R, 5(c+d)/4) = 5R. A grid over
  R ≤ c, d ∈ [c, 21c] gives worst F/c = −0.117977 = −γ exactly. Raising γ needs a new constraint on d relative to c, or on L_c beyond (4) in the d ≫ c regime, e.g. from the numerator's
  actual length. The large-width end binds through timbalud's τ_s = 13/400. With slacks 5.1·10⁻⁶ and 6.2·10⁻⁶, this whole lever is worth at most the 3.4·10⁻⁵ affine ceiling.
- **All twists mod m give the same fourth moment (NUMERICAL, mmnum/twistall.py, twistall_1e6.log; U = 10⁶, X = U^{0.417}, 3.6·10⁶ rows, 48 columns, 07:14Z).** For every character ψ of
  (O/m)^× (N(m) = 7, 13, 19; 6, 12, 18 characters), T(ψ) = E_v|Σ w w ψ(n₁n₂)χ_{n₁n₂}(v)|² agrees with the untwisted moment on columns coprime to m to ±0.002. That holds for the
  sextic-type ψ (orders dividing 6: these are exactly the thin families, rows vm^j) and equally for the foreign ψ (orders 4, 9, 12, 18). The ratios 0.729 and 0.850 are the column mass
  coprime to m ((1 − 1/7)² = 0.735 and (1 − 1/13)² = 0.852). Structural meaning (PROVED identity): T(ψ) = Diag + Off(ψ), with Diag independent of ψ and Off(ψ) = Σ_t F(t)ψ(t), where
  F(t) = Σ_{N≠N′, N≡tN′ (m)} c_N c̄_{N′} Σ_v χ_{N}χ̄_{N′}(v). Equality of all T(ψ) means F(t) ≈ 0 for t ≠ 1: the off-diagonal kernel carries no information about N/N′ mod m.
  So nothing about the thin twists is special in the data; the special role of the sextic-type ψ (row shifts, conductor at reflection) is purely an artifact of the proof. The target
  statement in its cleanest form is: the plain fourth-moment off-diagonal, restricted to pairs with N ≡ tN′ (mod m) for t ≠ 1, is o(Diag/Z^{0.36}) for small slot moduli m.
- **Asymptotic route, "thin = full − coprime" (PROVED reduction, 07:17Z; circular).** A power-saving asymptotic Σ_u|P(u)|² = U·Diag + O(U^{1−η}) for the full family, combined with the
  same for the rows coprime to m, would give the thin family Diag + O(ZU^{−η}Diag), natural for Z ≤ U^η. The secondary (theta, X^{5/6}-type) terms are below the natural thin size
  U/Z for z < 1/6. But the coprime restriction is invisible to multiplicative row characters mod m (they see only m ∤ u), and additive row twists e(bu/m) are what detect m | u. After
  Poisson in u, the complete sum mod mNN′ gives a coset condition h ≡ −bNN′ (mod m) times G(NN′, hm̄) = χ_N(m)χ̄_{N′}(m)G. Separating the bilinear e(cbNN′/m) uses λ(N)λ(N′) over all
  λ mod m, and that lands exactly on the column-twisted dual norms with twists λ(N)χ_N(m), i.e. the "all twists mod m" family. Numerically every member is natural (twistall above);
  for a proof each member pays the twist rule. So the asymptotic route reduces to the same decorrelation statement (F(t) ≈ 0 for t ≠ 1) and adds nothing new.
- **Hadamard-product form of the averaged thin family (PROVED identity, 07:20Z; circular).** In column (dual) form the averaged thin family is c*(A∘B)c, where A(N,N′) = Σ_{v≤V}χ_Nχ̄_{N′}(v)
  (Gram matrix of the row family) and B(N,N′) = Σ_{m~Z}ξ̄_Nξ_{N′}(m) (Gram matrix of the slot family, rank ≤ Z, eigenvalues about X²). Schur's inequality gives ‖A∘B‖ ≤ Z‖A‖, but ‖A‖ is the
  general-coefficient norm, which carries the copy loss. Evaluating at the actual plain c gives either Σ_m (c∘ξ̄_m)*A(c∘ξ̄_m), the twist rule per m, or Σ_v (c∘χ(v))*B(c∘χ(v)), which is natural iff
  the vectors c∘χ_·(v) are equidistributed against the Z-dimensional span of {ξ_m}, i.e. the per-m thin family again. Same decorrelation statement in a third guise.
- **Why exponent 1 would be free, and why the slots never give it (PROVED, 07:22Z).** A thin family with rows v·m (exponent 1) is natural with no input at all: χ_n(v)χ_n(m) = χ_n(vm), and the pairs
  (v, m) cover every row w ≤ VZ with multiplicity #{m | w, m ~ Z} ≪ log, so Σ_mΣ_v|P(vm)|² ≪ Σ_w|P(w)|²·log, natural by Lemma 18.1. Exponent 2 (rows vm²) covers only a 1/Z fraction of the rows
  up to VZ², so positivity loses Z. The exponent is 2i, where i is the slot prime's valuation in the column (G(p^iu, h) ∝ χ_p(u)^{2i}). It is always even, hence a cubic twist.
  Conjugate slots give i = 5, exponent 4. Cubed slots (i = 3, exponent 0) would carry no twist, but Σ_pψ_u(p)³ is a prime sum of the quadratic character ψ_u³ and does not detect
  zeros of L(s, ψ_u). Q³ and |Q|² put the twist on products of distinct primes, and exponent 0 occurs only on the trivial diagonal. So the parity of the slot exponent is fixed by the detector.
- **Diaconu–Ion–Paşol–Popa, arXiv 2607.27131 (CITED, abstract read 07:24Z).** They give asymptotics for twisted first and second moments of r-th order Hecke L-functions (r ≥ 3) over fields
  containing the 2r-th roots of unity, using multiple Dirichlet series. For r = 6 that needs the 12th roots of unity (Q(ζ₁₂)), so the sextic family over Q(√−3) is not covered; r = 3 is.
  Their "twists" are of the form Σ_χ|L(½, χ)|²χ(a), i.e. multiplicative row twists, which see only m ∤ u, not the thin condition m² | u. They are also full L-value moments, while we need a plain
  fourth moment at product length 0.834. Not an input as stated. A sextic version with twists by a fixed character of conductor m, uniform in N(m), would be one.
  Correction of degree (07:25Z): |S₁S₂|² with two plain sums is a fourth-moment-type quantity. Through Perron each S_i is a t-integral of L(½ + it, χ), so the thin family needs a twisted
  fourth moment Σ_v∫|L(½+it, χ_vξ_m)|⁴ uniform in N(m), not a second moment. The paper's Lemma 18.1 is the untwisted fourth moment, natural for all lengths. Multiple Dirichlet series
  methods (DIPP) reach second moments only, so even a sextic version of 2607.27131 with fixed-character twists would not supply input 2. The needed statement is a twisted fourth
  moment of sextic L-functions with conductor dependence N(m)^{<0.36}.
- **A structured version of de Faveri's recursion (HEURISTIC programme, 07:26Z).** In the proof of Theorem 1.1 the only excess over A + B comes from dual rows a = d·e⁶ with e ≥ 2: copies,
  weighted by |Σ_bλ_bg̃(b)χ̄_b(d)|² at small d. So the loss is driven entirely by the correlation of λ (alternately λ and λg̃ along the recursion) with characters of small rows. For
  λ = c∘ξ_m (plain coefficients twisted by the cubic slot character):
  - the odd levels see twisted smooth theta sums Σ c_bξ_m(b)g̃(b)χ̄_b(d), damped by N(m)^{−1/6} as measured (twist2_1e5.log);
  - the even levels see Σ c_bξ_m(b)χ_b(d), which is a short character sum except at d = m⁴f⁶, where χ_b(m⁴)ξ_m(b) = χ_b(m⁶) gives full correlation. That alignment is the thin family's own
    main term (the analogue of the zero frequency for plain c), not an excess.
  So the programme is: run the Θ₆ recursion for this specific λ, book the m⁴f⁶ rows as main term, and bound every other excess by twisted theta sums uniformly in N(m). If the
  theta damping propagates through the levels, the averaged thin family would come out natural. Each level needs a uniform-in-m error term for twisted Kubota series, the same
  input §94 isolated, now with an explicit place to use it. Not carried out.
  **Programme check (exponent count, 07:27Z): circular as stated.** Pointwise convexity for the twisted Kubota sums (|B_m(d)|² ≲ X²N(d)^{1/2}) makes the odd-level excess VZ·W^{3/2}, i.e. U^{0.168+3z}
  above natural. So the programme needs the mean square Σ_{d≤D}|B_m(d)|² of twisted theta sums with Lindelöf-on-average quality, uniformly in m. That mean square is exactly the paper's
  Gauss-row norm N_C with the twist, which the paper's second transform reduces back to plain coefficients and the twist rule. De Faveri's power-free kernel organises the copies but does
  not supply this mean square. So the structured recursion is a reorganisation of the same two-transform argument and inherits the same N(m) loss.
- **Gao's large-width binding constant τ_s = 13/400 is c_I = 2 again (CITED: submissions/detector-affine-radical-c0125161, VariableSupply.lean, 07:28Z).** inverseWidth S = min(S, 13/400) comes from
  inverse_capacities: r ≤ 467/500 and r + 2s < 1 give s < (1 − 0.934)/2 = 0.033, i.e. Lemma 17.1's cost 2. plainWidth = min(S, 3/200) comes from 2m + 6(3/4 + 2Δ)s ≤ 1 with m ≤ 233/500,
  i.e. Lemma 18.1's 6κ. So both ends of the affine construction bind through the same two slot costs. The affine ceiling is not an independent lever: raising τ_s means c_I < 2 inside the
  HIGH source bins, and raising the plain width means a plain slot cost below 6κ. Both are the open inputs.

## 105. The unconditional Patterson paper (OpenAI, 25 Sep 2026) against the four inputs (CITED reading of lit1017e/patterson.txt §1–§5; PROVED exponent remarks)
- **What it proves.** Σ_{Nπ≤X} G(π) = (6/5)c_* X^{5/6}/log X + o(X^{5/6}/log X) for normalised cubic Gauss sums at primes, plus o(·) cancellation for every fixed angular mode ℓ ≠ 0. The error is little-oh, not a power.
- **New mechanism (§4, Prop. 4.1).** For prime-structured β, Σ_{p,q}|Σ_b β(b)χ_b(vpq²)|² ≪ Y^{7/3−η}, with η > 0 ineffective (limiting-exponent contradiction), at the two configurations where Heath-Brown's cubic sieve and the ordinary character sieve are sharp: (P,Q) ≈ (Y,1) and P = Q ≈ Y^{1/3}. The latter is the Dunn–Radziwiłł Type III wall of §95. Steps: Heath-Brown-type convolution identity; reflect a dominant full sum; moments of several lengths force rigid value exponents v_i = 5a_i/6 on Y^{2/3} rows; Bombieri–Halász–Montgomery with the off-diagonal Gram bound Lemma 4.3 gives the contradiction.
- **Fit against our inputs (PROVED exponent remarks).** The mechanism removes only the extra theta/copies terms of a sieve, and only by an ineffective η.
  - M4 (Möbius fourth moment, §88). L₄ < 0.117 lies below even the perfect large sieve (L₄ = 1/2 at r = 3/4, column term). Out of reach.
  - LM and RD (inputs 1 and 3). These are asymptotics, self-dual under Poisson (§22, §70), and need a saving well beyond the sieve's extra terms. The rigidity argument also needs the μ-specific bound at every length, but powers of factors are no longer Möbius polynomials. Not applicable.
  - Input 2 (★). For z ≤ 0.083 it is a "beat the sieve's extra terms" problem, but the saving needed is U^{0.06–0.10}, i.e. most of the way to natural, not η. No usable bounds exist at the other lengths (the sixth moment is open), so the rigidity step has no input.
  - c_I: the cubic-sieve extra term enters only case (b) of §32, which is not binding.
- **Side remark (PROVED given the 7/8 theorem).** Because ζ_K has no zeros for Re s > θ, the theta model term of Σμ(c)g̃(c)χ(c) is ≪ X^{θ−1/6+ε} = X^{0.708}, below §79's pointwise requirement D^{0.773}. The rest of that requirement, a saving of D^{0.06} below 5/6 for the fluctuating part, would need an explicit and large η, which the method does not supply.
- **Leads it leaves.**
  - Dunn–Radziwiłł's unconditional level-aspect Voronoi formula (their Theorem 2.4 here) is the level-m functional equation for twisted cubic Kubota sums. That is the analytic input the single-weight sublattice bound of OPEN_PROBLEM_DECORRELATION.md needs.
  - "Free quadratic twist" route to plain cost 3 (θ ≈ 0.874313). Rows hm⁴ = (hm)·m³. Drop hm ∈ mO by positivity, losing Z; what remains is the sextic family twisted by the quadratic character ε(m) = χ(m)³. If Lemma 18.1 held for quadratic-twisted smooth coefficients with no loss, cost 3 would follow. This is open.
