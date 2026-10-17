# Attempt at θ = 13/15 (sextic architecture), session 3

Labels: PROVED, MODEL (exponent algebra of the 7/8 architecture), NUMERICAL, HEURISTIC, OPEN.

## 0. What 13/15 requires (recap of REPORT.md §3, §9, §13.4)

- 13/15 is the ceiling of the probe architecture: the generic-row condition C0 < 0 forces ℓ < 1/5, and θ = 11/12 − ℓ/4.
- Reaching it needs a density-hypothesis-strength row count for the sextic twist family {ηχ_u : N(u) ~ U}:
  #{u : L(s, ηχ_u) has a zero with Re s ≥ σ, bounded height} ≪ U^{2(1−σ)+ε} on σ ∈ [0.51, 0.867].
  Conditional Theorem D (REPORT §3): DH with loss λ gives θ = 11/12 − (1−2λ)/(4(5+6λ)); λ = 0 gives 13/15.
- Equivalent sufficient inputs found earlier: MM′ (Möbius × plain mean square at total length up to U^{3/2}),
  or the Möbius fourth moment M4(r) for all r ≤ 1, or plain moments of every order. All are supercritical
  Lindelöf-on-average statements for a thin family (family size ≈ conductor).
- The required density exponent in σ ∈ [0.63, 0.79] is below Guth–Maynard's exponent for ζ there (REPORT §8.1).

## 1. New route tried: bootstrap from the contradiction hypothesis

Inside the contradiction argument every Hecke L-function over F = Q(√−3) is zero-free in Re s > β*.
This is information the row counts do not use. Question: does it improve any row count?

### Lemma B (PROVED, standard tools)

Let ψ be a primitive finite-order Hecke character over F of analytic conductor q, and suppose L(s,ψ) ≠ 0 for Re s > β*.
Put σ₁ = β* + ε and, for ½ ≤ σ ≤ σ₁, λ = (σ − ½)/(σ₁ − ½). Then for every q₀ > 0

    |L(σ+it, ψ)|^{q₀/(1−λ)} ≤ (q(|t|+2))^{O(ε)} · ∫ P(σ, t−v) |L(½+iv, ψ)|^{q₀} dv,

where P(σ,·)/(1−λ) is the (normalized) harmonic measure density of the line Re s = ½ for the strip ½ < Re s < σ₁,
which decays like exp(−π|t−v|/(σ₁−½)).

Proof.
1. Lindelöf beyond the zero-free half-plane: log L(s,ψ) is holomorphic in Re s > β*. Borel–Carathéodory in discs
   centred on Re s = 2, followed by Hadamard's three-circles theorem (Titchmarsh §14.2 with ½ replaced by β*), gives
   log|L(σ₁+it, ψ)| ≤ ε′ log(q(|t|+2)).
2. log|L| is subharmonic and of finite order in the strip, so (Phragmén–Lindelöf) it is bounded by its harmonic
   majorant: log|L(σ+it)| ≤ ∫P(σ,t−v) log|L(½+iv)|dv + ∫P₁(σ,t−v) log|L(σ₁+iv)|dv, with ∫P = 1−λ, ∫P₁ = λ.
3. Bound the second integral by step 1. For the first, Jensen's inequality for the probability density P/(1−λ) gives
   exp(∫P log|L|) ≤ (∫(P/(1−λ))|L|^{q₀})^{(1−λ)/q₀}. Raise to the power q₀/(1−λ). ∎

**Corollary B1 (CONDITIONAL on the family fourth moment at ½, which follows from the paper's plain lemma via the
approximate functional equation).** For the retained rows u ~ U and σ = ½ + λ(σ₁ − ½):

    Σ_u |L(σ+it, ηχ_u)|^{4/(1−λ)} ≪ U^{1+ε}.

So every moment is Lindelöf-on-average strictly to the right of ½ (e.g. the sixth moment at σ = ½ + (σ₁−½)/3).
This is genuinely supercritical information, and it is free inside the contradiction argument.

### Proposition B2 (PROVED, exponent algebra): the bootstrap does not improve any detector count

Write b = β − ½ for the real part of the detected zero and b₁ = σ₁ − ½, so b ≤ b₁. Witnesses are evaluated at
real part β. Moving a witness from the critical line to the line ½ + λb₁ (where Corollary B1 applies) costs a
factor N^{λb₁} in the Mellin shift and gains the higher moment.

- **Plain witness** S_u(N), N = U^m, threshold N^{b}: the count is U^{1 − m·4(b − λb₁)/(1−λ)}.
  d/dλ of the saving is 4m(b − b₁)/(1−λ)² ≤ 0, so λ = 0 (the existing plain fourth-moment count 1 − 2δm) is optimal.
- **Mixed witness** (detector integral term, |L·M_X| at the shifted line against Y^{β−α′}), Hölder with Corollary B1
  and the Möbius second moment: saving 4(b − λb₁)/(3 − λ), derivative 4(b − 3b₁)/(3−λ)² < 0. Again λ = 0.
- **Halász–Montgomery off-diagonal.** Lindelöf beyond β* for the pair characters χ_uχ̄_{u′} (conductor ≍ U²) gives
  Gram entries ≪ N^{σ₁}U^{ε}. A count needs N^{δ} > N^{σ₁}, i.e. 2β − 1 > σ₁, impossible since 2β − 1 < β ≤ σ₁.
- Pointwise convexity between ½ and σ₁ (μ(α) = (σ₁−α)/(4(σ₁−½))) in the detector's integral term gives a count
  U^{1−δ+(2α−1)+2μ(α)}, minimised at α = ½, i.e. U^{3/2−δ}: worse than the existing counts.

**Reason (structural).** Every input derived only from "no zeros right of β*" is scale-compatible with a zero at
real part β ≤ β*: interpolating towards σ₁ ≥ β loses at least as much as the stronger moment gains. A row count
below the critical-line counts needs information at real parts below β, which is exactly zero-density information.

So the contradiction hypothesis cannot be recycled into the missing density input. The rows need an external
zero-density theorem for the sextic family.

## 2. Quartic route as an alternative path to 13/15 (MODEL)

- Quartic bookkeeping (NOTES §3, §6): σ₀ = 7/8 − 3ℓ/8 with completed-row range M = 1 − ℓ − g, and with only the
  large-sieve shape σ₀ = 7/8 + ℓ/8 + 3g/8 ≥ 7/8.
- Reflection plus an optimal dual sieve reaches M ≤ ½ (Gao–Zhao unconditionally: M ≤ 3/8). That is exactly the
  no-gain regime. Any M > ½ lies in the reflection core around the fixed point 2/3, where Theorem ZF shows the
  needed bound already implies a quasi-RH strip for the quartic Hecke family over Q(i).
- So the quartic path to 13/15 is blocked by a quasi-RH-strength input, the sextic path by a DH-strength input.

## 3. Literature check (Oct 2026 search)

No unconditional density-hypothesis-type zero-density estimate for cubic or sextic Hecke families was found.
Known results for these families are moments, one-level density under GRH, and nonvanishing proportions
(Gao–Zhao, David–de Faveri–Dunn–Stucky, Lin–Wong 2025).

## 4. Status

- θ = 13/15: OPEN. The single missing input is DH-strength zero density for the sextic twist family on
  σ ∈ [0.51, 0.867] (Conditional Theorem D then gives 13/15).
- New this session: Lemma B / Corollary B1 (supercritical moments off the line, free inside the contradiction
  argument) and Proposition B2 (they give no gain), which closes the self-bootstrap route.

## 5. Round 2 (after "do not stop")

### 5.1 Block reflection of the Poisson rows (NEGATIVE, structural)
Idea: bound Σ_{u∈B} (row u) for a dyadic block B by reflecting the block-restricted probe, giving cancellation across
rows instead of a density count. Read against paper.tex §"The exact Poisson series" (Sept 30 preprint, lines 3723–3847):
Poisson is in m modulo s·b_*·A with A = cn³ a summation variable, and the frequency is H = ua⁶. Restricting u to B
inserts a kernel in m/(s b_* A). In the low representation (eq. low-separated) the A-sum must stay a completed theta
row Σ_A (coeff)·χ_A(m) for the cubic reflection; a factor depending on m/A destroys that structure. So the
block-restricted probe has no low-side representation, and bounding the non-principal rows by reflection is
equivalent to bounding the principal row (circular, as in REPORT §7(b)).

### 5.2 Quartic route at 13/15 needs (G) at M = 8/15 (MODEL)
σ₀ = 1 − M/4 (L3_REPORT §1) gives M = 8/15. Kubota reflection M ↦ M/(3M−1) sends it to 8/9. Then (G) at 8/9 implies
PO at η = 2 − 8/9 = 10/9, and Theorem Z gives zero-freeness of every quartic Hecke L(s,(k₀/·)₄) over Q(i) for
Re s > max(1 − η/8, ½ + 3η/8) = 0.9167. So the quartic input for 13/15 already contains a new quasi-RH (0.917) for a
family not covered by OpenAI's theorem, and GRH does not give it by any known argument (L3_REPORT §10b).

### 5.3 October 5 route (MODEL)
θ = ½ + 5c/12 from MS(D^c). 13/15 needs c = 0.88 < 1, i.e. square-root cancellation for Möbius sums over fewer
characters than the length (REPORT §7(b)); by the explicit formula this needs zero-freeness near ½ on average,
i.e. stronger than DH.

### 5.4 Smallest input that moves θ at all
REPORT §11: the plain sixth moment Σ_u|S_u(U^{m₁})S_u(U^{m₂})S_u(U^{m₃})|² ≪ U^{1+λ}, m_i ≤ 0.45, λ ≤ 0.12.
- From the plain fourth moment plus the pointwise bound |S_u(N)| ≤ N^{β*−½+ε} (Lemma B, inside the contradiction
  argument): λ = 2(β*−½)m ≈ 0.75m = 0.34 at m = 0.45. Too large by a factor of about 3.
- What is missing is a large-values bound #{u : |S_u(N)| ≥ V} ≪ U^{1.12}V^{−6} for V ∈ [U^{0.06}, U^{0.17}],
  beating the fourth-moment count U V^{−4} by U^{0.12}V^{−2}. Halász–Montgomery fails (off-diagonal pair conductor U²,
  Gram entries ≥ N^{β*}). For quadratic characters the unconditional sixth moment is also not known (Gao 2021: sharp
  bounds for k ≤ 2; Shen–Stucky: fourth moment asymptotic).

### 5.5 Net (round 2)
Every route to 13/15 found so far needs one of: DH for the sextic twist family on [0.51, 0.867] (sextic route),
a new quasi-RH ≈ 0.917 for quartic Hecke L over Q(i) beyond GRH-methods (quartic route), or beyond-DH Möbius
cancellation (Oct 5 route). Each is a new theorem beyond the current state of the art.

## 6. Round 3: critical-length levers (no density hypothesis)

### 6.1 Numerics: plain sixth moment (NUMERICAL, sextic/mmnum/p6x.py)
Rows u sixth-power-free, N(u) ≤ U, excluding u supported on λ (units, powers of 1−ω: principal-type rows, |S| ≈ 2.9).
M6/(6 M2³) at m = 0.5: 0.867 (U=1e4), 1.055 (1e5), 1.126 (1e6); random multiplicative control 0.857, 1.059, 1.133.
Mixed (0.3,0.4,0.5): 1.069, 1.063, 1.332 vs control 1.046, 1.062, 1.338. Max row share ≤ 0.03%.
So the sextic family behaves exactly like the random multiplicative model: no power-size excess up to U = 10⁶.

### 6.2 Value of the two slot costs in the record's framework (MODEL, c61/capc3.py, sens)
Plain lemma: n1+n2+cP·κ·z ≤ M (paper cP = 6). Marked inverse lemma: r + cI·z ≤ m (paper cI = 2).
| change | θ (grid model; baseline 0.87492 ≈ record) |
|---|---|
| cI = 1.98 / 1.95 / 1.9 / 1.8 / 1.7 | 0.874908 / 0.874777 / 0.874714 / 0.874445 / 0.874185 |
| cI = 1.5 / 1 | 0.87328 / 0.87289 |
| cP = 4 / 3 / 4/3 (natural) | 0.87435 / 0.87371 / 0.87179 |
| both natural (cP = 4/3, cI = 1) | 0.86979 |
Both are critical-length statements (total length ≤ conductor), unlike the density-hypothesis inputs.
Both natural capacities together reach ≈ 0.8698; the rest of the way to 13/15 still needs density-hypothesis input.

### 6.3 Where cI = 2 comes from (PROVED, reading of paper.tex 9783–11130)
- Canonical invariant (i) F − M − Q − z₀ ≥ c (F = N+V, N = r+z columns, M = 2N − m dual rows) gives m − r − 2z ≥ c.
  Slots cancel from F − z₀, so (i) says: dual rows ≤ non-slot column + label − puncture.
- Terminal reflected energy, row branch: E_ref = O/2 + N*/2 + v/2 − u + … ≤ (M+N)/2 needs only M ≤ F (z_a cancels).
- Column branch at maximal dual: T_d = 2H + 2A₀ + 2z_a − N* − …; active marked primes enter the reflection
  conductor (Corollary marked-reflection) and lengthen the dual by 2z_a, and u ≤ z_a, so E ≈ 2M + 2z_a − N
  forces M + z_a ≲ N. That is the second copy of z.
- Recursion: F_c − M_c = F − M + D_c and Q_new ≤ Q + D_c + 4η, so a deficit in (i) is never repaired by Poisson
  steps (only branches with A_i + V > 0 gain). Hence cI = 2 is intrinsic to detecting slot primes through the
  cubic-theta reflection.

### 6.4 Where cP = 6κ comes from (PROVED, reading of paper.tex 12477–14985)
- Each recursion step enlarges the Gauss-row ball by g ≥ J₊ + σ to force the width drop M' ≤ M − σ; this raises
  A − M by about σ (child excess). Zero-slot children are rescued by Hecke reflection (any lengths).
- Positive-slot children cannot reflect slots, so the excess is paid by removing slots pointwise (|Q|² ≤ Z^{κz}),
  financed by an allowance ℓ from averaging the row over h ↦ hp⁶ (pool of Z^ℓ primes), which itself costs a row
  extension 6ℓ. Excess 6ℓ per allowance ℓ ⇒ removal must lower the affine form 6× faster than its κ-cost ⇒ cP = 6.
- Sixth powers are the only frequency multipliers fixing every column Gauss sum, so the ratio 6 is forced for
  linear row amplification. Cycle61's affine trick (one top-level amplifier, λ ≈ 4.03 in a narrow window) is the
  only known local improvement.

### 6.5 Numerator line for bad rows (NEGATIVE, exponent algebra)
The high exponent's per-row factor U^{δ/2} comes from the numerator L(w,χ_•(u)) on Re w = 1−a (functional equation).
Moving w to ½ and using the plain fourth moment by Hölder over the bad set B changes the row term by
Z^{h(1−R)/4 − (h−y)δ/2}. At the binding (δ, R, h, y) ≈ (0.467, 0.594, 0.795, 0.4595) this is Z^{+0.0024}: no gain
(a sixth moment would give one; that is the supercritical input again).

### 6.6 Cycle61's affine trick, analysed (PROVED algebra from c61paper.tex §2.2)
Theorem 2.2 = one top-level p⁶ amplifier of width a = K_a·g, then children split by zero parameter:
low bins (δ ≤ q₀): delete slots pointwise at rate 2q₀/9, reserve a(1−4q₀/3) − (2q₀/9)g;
high bins (δ ≥ q₀): pointwise child bound + Lemma 2.1 density ρ(δ)=1−10δ/(3(5−2δ)), reserve
1 − ρ₀ − q₀A − a(5−6ρ₀). Eliminating a gives the affine law n + λz ≤ B with
λ = (q₀+4.5K)/(q₀+K), B = (1−ρ₀+K)/(q₀+K), K = (2q₀/9)(5−6ρ₀)/(1−4q₀/3).
Capacity at the binding m = 0.417 is maximal at q₀ ≈ 0.66 (0.0384; Cycle61's 0.68 is essentially optimal; old 0.0369).
Both reserves are monotone in δ, so two bins are optimal for these two bounds.
Child density improvements: ρ − 0.01 → 0.0391, ρ − 0.03 → 0.0415, DH for children → 0.0692.
So the affine mechanism is near its limit; large plain-capacity gains need a new mechanism.

### 6.7 Large-value tails of Möbius polynomials (NUMERICAL, mmnum/lvstruct.py, lvctrl.py)
U = 1e6, r = 0.75: P(|M|² ≥ t·mean)/e^{−t} = 1.23, 3.12, 11.9, 63.8 at t = 4, 8, 12, 16; the random multiplicative
control gives 1.23, 3.05, 10.8, 46.6. max/mean = 19.2 vs 20.7. The heavy tail is generic multiplicative chaos;
μ shows no extra structure and no power-size large values (consistent with DH and far beyond).
BGL Corollary 1.6 density (A = 8/(7−6σ), 2.86 at σ = 0.7) is weaker than the current counts (2.275).

### 6.8 Family-averaged probe as a density detector (NEGATIVE, structural)
Idea: average the probe identity J_w = f_w + rows_w over targets w = ηχ_u to get density for the targets.
- Twisting the target shifts the theta-row index: B_{m,ηχ_u} ≈ B_{mu,η}, so the low side averages a
  multiplicatively shifted correlation ⟨A, B_{·u}⟩; with random-like rows the diagonal gains a factor Q.
- But the rows of target w are signals of targets w·u'^{∓1} weighted by numerator factors N(u') independent of w:
  J = (I + K)F with K a convolution on the target group. On the dual (ideals n) this is multiplication by 1 + N̂(n),
  and F̂(n) is just μ(n)η(n)n^{−s}×|family|. The averaged identity therefore only re-encodes the multiplicativity
  of μ; the zero information of individual targets is lost. Inverting I + K is equivalent to the original problem.
Conclusion: density for the rows cannot be extracted from the probe by averaging over targets.

### 6.9 Möbius fourth moment via the canonical recursion (PROVED limit)
On squarefree n the coefficient of M_u(D)² is μ(n)·d_D(n), and d_D enters the canonical family only as a full-length
mark of length r (product form via Mellin in the divisor). Invariant (i) with N = 2r, z₀ = r gives 1 − 2r − r ≥ 0,
i.e. r ≤ 1/3: exactly the proven range. M4 beyond 1/3 needs a non-mark treatment of the divisor weight.

### 6.10 Amplitude classes in Theorem 2.2 high bins (MODEL, c61/ampclass.py, ampgain.py)
Capacity at m = 0.417 rises 0.03844 → 0.03960 (m = 0.40: 0.0474 → 0.0501), but θ in the Cycle61 grid model is
unchanged (0.8748997 with and without): the binding m leaves the affine window. Plain capacity helps only if global.
Slot primes kept outside the Möbius–Gauss conversion still enter the n₀-reflection conductor through the
reciprocity factor (n₀/p)₆, so the inverse cost cI = 2 cannot be avoided that way either.

## 7. Work queue (continuation)
Q1. Plain sixth moment with loss: interpolate the plain 4th moment with a lifted 12th moment (sixth-power lift +
    BGL) and check whether any c < 0.6 in loss c(3m−1) is reachable for m ≤ 0.45.
Q2. Inverse long branch: test amplification by paired residue characters (sextic = cubic × quadratic) for Möbius
    polynomials of length U^s, s ∈ (1, 1.3]; target slope α′ ≤ 0.55.
Q3. Mixed moment slightly supercritical with loss c(s−1): what c moves θ, and what does the sixth-power lift give.
Q4. Record squeeze in the Cycle61/radical framework: re-optimise reserves and rationals at the current binding.

## 8. Continuation 1 (12:04Z)
- Q1 (lossy plain sixth moment via lifted 12th moment): sixth-power lift + BGL gives Σ|S³|² ≪ U^{(7/2)m+1/6} (U^{1.74}
  at m = .45); unlifted BGL gives (UL)^{2/3} = U^{1.57}. Both worse than the pointwise route (loss 0.34). DEAD.
- Q2 (long-branch amplification): Vaughan/Heath-Brown decomposition: Type I pieces vanish (smooth length > conductor),
  Type II pieces are supercritical bilinear forms with arbitrary coefficients (§14.2: false in general), so they need
  μ-arithmetic at supercritical length. Linear copy amplification already optimal (§7a). DEAD.
- Q3 (slightly supercritical mixed moment): requires critical MM first, which is itself unproven (theta reflection
  needs r + 2m < 1, empty at r + m = 1). DEAD.
- New: averaged probe as density bootstrap (MODEL). Targets ηχ_u, u ~ U, probe scale Z = U^c (c ≳ 1 needed for the
  target's conductor). Even granting ideal orthogonality of shifted theta rows and ignoring rows:
  D(a) ≤ 1/2 + c(σ₀ − a). At a = 0.73, c = 1: 0.645 > current count 0.594. Rows add a term needing density at the
  larger level U·Z^d. NEGATIVE.

## 9. Continuation 2 (continuous work, 12:10Z–)

### 9.1 Inverse slot cost: exact origin (PROVED, reading of paper.tex 7657–8060)
Reflected energy E_ref = O/2 + max(H, v+ℓ_b) + [max(v, z_a, 2(v+z_a)/3) − v] − S₀ − B₀ − ℓ_b − 2e/3 − ½(T_d−v−3ℓ_b−e)₊.
The bracket is Heath-Brown's cubic large sieve over active slot tuples P (moduli Z^{z_a}) against the dual n (Z^v):
slots are free whenever v ≥ 2z_a. The cost appears in the column branch: there E ≈ O/2 + T_d with
T_d − S₀ − B₀ ≤ 2M − O + Q + V − N* + 2z_a (active slot primes are in the reflection conductor c = c_F F_act R P),
so E ≤ F needs M + Q/2 + z_a ≲ N; with N = r+z, M = 2N − m this is r + 2z ≤ m. Pulling the slot prime out of the
column leaves the cubic twist χ_{n'}(P)⁴ (mod P), which again puts P in the conductor. Intrinsic to theta reflection.

### 9.2 Record framework is locally optimal (MODEL + handoff notes)
Current record 0.8748389 is within 5·10⁻⁸ of the float optimum of its own constraint system; binding is the new-h row
at δ ≈ 0.467, p = 1/2. At d = h: G(h) = (1+δ)/2 + ℓ/2 + ℓδp − 5h/12 + (3/4)hR − (δp/2)C_h (z-contour 17/50 cancels).
The weighted numerator (Hölder with the plain 4th moment, coefficient ¾R + ¼) is optimal given the per-row pointwise
bound |L(½)|² ≤ U^{δ/2} on bad rows: total(k) = R + δ/4 + (1−R−δ)/k is minimised at k = 4 because 1−R−δ ≈ −0.06 < 0.
A numerator sixth moment over all rows with loss λ would help iff λ < (1−R)/2 ≈ 0.205; available λ = β*−½ ≈ 0.375.

## 10. Continuation 3: Option 1 (DH for the sextic family), Guth–Maynard and large-sieve round
All details are in DH_SEXTIC.md §7–§15. Summary:
- Price list at the binding class (δ = .466, p = ½; t* = 1.162, r* = .745, m* = .417, three-way tie 0.5935):
  Möbius count −0.01δ → θ −0.00024; plain −0.01δ → −0.00011; LM slope 0.80 → −0.00025; LM slope 0 → exactly 13/15.
- Benchmark: the family's count R(δ,½) already beats the ζ analogue of Guth–Maynard at every σ in [0.70, 0.80], and
  its near-line slope (≈ 0.82) beats Conrey's mollifier slope (4/7). DH in the window is therefore stronger, relative to
  family size, than any zero-density exponent known for ζ at the same σ.
- Guth–Maynard and Chen–Gupta–Li (q-aspect) do not transfer: Gauss sums of the pair characters factor
  (τ(ψ_{u1}ψ̄_{u2}) = ψ_{u1}(u2)ψ̄_{u2}(u1)τ(ψ_{u1})τ(ψ̄_{u2})), so the cubic trace is exactly self-dual,
  M_B M_B^* = N·I + (N/U)D(M′_B M′_B^* − H·I)D^*, H = U²/N, and no free archimedean/residue variable remains.
- Heath-Brown identity, all-t covering (§12) and Hölder recombinations reproduce R exactly. The single-t count is lossless.
- Large sieve: OpenAI's sextic sieve is BGL K + D + (KD)^{2/3}. Dunn–Radziwiłł / de Faveri–Dunn–Hoffstein: (AB)^{2/3} is
  attained by Gauss-twisted sequences (the reflected-column case). A hypothetical sharp sieve gives at most θ ≈ 0.8696.
- r-th order moment literature (DdFDS 2410.03048, DIPP 2607.27131): twisted second moments with short twists only;
  DIPP needs μ_{2r} ⊂ F, which excludes the sextic family over Q(ω). MM′ needs Möbius twists of length U^{0.745}.
