# Open problem LM: Möbius randomness beyond the conductor for sextic residue symbols

DaysSky, with Claude assistance. 16 October 2026. Details: `DH_NOTE.md`, `DH_SEXTIC.md` §1–§79. Conditional theorem: `DH_THEOREM.md`.

## Statement

Let O = Z[ω]. For a nonzero element u, put ψ_u(n) = ν(n)(u/n)_6^{±1}, where ν is a fixed finite-order character, and fix a
smooth weight W supported in [1/2, 2]. For D ≥ 1 let

  M_u(D) = Σ_n μ(n) ψ_u(n) W(N(n)/D),

with n running over ideals (primary generators). Write Σ* for a sum over **sixth-power-free** u with N(u) ≍ U.

**(LM)** For every fixed s ∈ (1, 3/2] and ε > 0:  Σ*_u |M_u(U^s)|² ≪ U^{1+ε} · U^s.

More generally, **(LM_α′)** is the bound Σ*_u |M_u(U^s)|² ≪ U^{1 + α′(s−1) + ε} · U^s.

**Certified values (17 Oct, CONDITIONAL_THEOREMS.md; exact rationals, relative to the radical record's proof).** LM_{4/5} ⟹ θ ≤ 0.87462; LM_{3/4} ⟹ θ ≤ 0.87425;
LM_{7/10} ⟹ θ ≤ 0.87368. All are below the frontier 0.874836, and copies give α′ = 5/6. LM_{α′} for s ∈ [1, 3/2] is the rows-dominant asymptotic (RD_{c₀}) with
α′ = 5/(6(1 + c₀)) (PLAN_RD.md).

## Why it matters

In the OpenAI quasi-Riemann architecture (zero-free half-plane Re s > θ for ζ and all Dirichlet L-functions):

- **Known:** LM_{5/6}, from sixth-power copies plus the critical-length Möbius moment. This gives the current record
  θ = 0.8748389 (our local record; the upstream frontier is now Gao's 0.874836, PR #62, which contains it).
- **LM_α′ ⟹ DH on σ ≥ (1+α′)/2.** The row count is R = 1 − δ + λ with λ = (α′ − δ)·δP/(2J), the closed form of
  Xu–Sang (2610.15291), which agrees with our model to machine precision. So the density-hypothesis gap is proportional
  to α′ − δ. DH on δ ≥ α′ also follows directly from the paper's no-slot detector at t = 3/2, which is its own endpoint step
  (20.7), with no selected primes (DH_SEXTIC §63). On that slot-free route, LM_α′ with α′ < 0.0544 already beats the
  record (DH_THEOREM.md §3a).
- **LM (α′ → 0) ⟹ θ = 13/15**, which is the architecture's ceiling, provided the paper's bin floor a0 = 51/100 is also
  lowered towards 1/2. With the floor at 0.51, LM gives 0.86979. The implication is stated and margin-checked against every
  estimate of the paper's §20 in `DH_THEOREM.md`. DH is needed only at row scales d ∈ [2/5, h+ζ]. Intermediate slopes give intermediate θ:
  α′ = 0.7 → 0.8737, 0.5 → 0.8706, 0.3 → 0.8684.
- **Any α′ < 5/6 at s ≈ 1.16 improves the model's baseline** (dθ/dα′ ≈ 0.009 there). Passing the actual record, which is
  0.00008 below that baseline, needs α′ ≤ 0.814 on [1, 1.16], or the equivalent U^{0.003} saving below.
- **Smallest form.** LM only just past the conductor (D ∈ [U, U^{1+c}]) already gives θ = 0.8732 (c = 0.05), 0.8714 (0.1) and
  0.8695 (0.16) (DH_SEXTIC §55). LM up to s = 1.16 gives DH for the family on σ ∈ (0.781, 0.875], and up to 1.3 on (0.659, 0.875]. A power-saving *asymptotic* there would also give α′ < 5/6 at every s, via copies plus
  subtraction of scales (§54). Even a partial slope counts: in the same model, slope 0.7 instead of 5/6 on [1, 1.05] gains
  0.00025 over the model's baseline, three times the gap between that baseline and the record (§55, MODEL). The frontier is
  an integrated saving of about 0.003 at the long-witness length. A bound U^{1+(5/6)(s−1)−0.003}·U^s for s ∈ [1.10, 1.25],
  i.e. U^{0.003} better than copies on that window, would already move the record in the model.
- **On the exact record pipeline (DH_SEXTIC §71, 16 Oct 23:50Z)** any saving counts. U^η better than copies on s ∈ [1.10, 1.25]
  gives θ = 0.8747953 / 0.8747080 / 0.8744044 for η = 0.001 / 0.003 / 0.01 (record 0.8748391; dθ/dη ≈ −0.044), and η = 0.03 on
  [1.05, 1.30] gives 0.8735620.

## What does not work (proved barriers)

- **Toolkit closure (Theorem A).** All available mean-value inputs and their Hölder combinations are realised by an
  explicit adversary, so the count is exactly R.
- **Self-duality (Theorems B, C).** The dual coefficients are μ(n)g̃6(n) = ε(n)(n/|n|)·conj(g̃3(n)), i.e. twisted cubic Gauss
  sums (§56). Poisson in u maps MS to UD + (U/D)Σ_h(|B_h|² − D), with
  B_h = Σ μ(n)g̃(n)χ̄_n(h)W. Poisson in h returns MS with coefficient exactly 1. Copies map slope α′ to min(α′, 5/6), so 5/6
  is a fixed point.
- **Decompositions (Theorem D).** Vaughan, Heath-Brown and Linnik identities all leave Möbius-heavy pieces, whose LM is a
  higher moment beyond the conductor.
- **Large sieve is not enough.** By Xu's lower bound, the general-coefficient sextic large sieve exceeds D^{5/6}U^{1/3} for
  D ∈ [U^{16/11}, U²]. That is worse than what copies give for μ, so any proof must use the arithmetic of μ.
- **Formulation.** Over all elements (not sixth-power-free), LM would imply a zero-free region Re s > 7/9 for the
  unit-type rows by itself. The sixth-power-free version is the right, weaker statement.
- **Zero dips.** The adversary needs M_u saturated at U^{0.745} and U^{1.162} but not at U¹. Closing this via
  Turán/Jensen/off-line moments is short by a factor of about 10 in moment order. A cluster of zeros contributes a
  divided difference; a dip needs off-line zeros spread over heights ≳ 0.2 (§40, §47).

## Analogues

- **ζ (t-aspect).** This is the power-saving variance of μ in short intervals, equivalently averaged Chowla with
  H = D/T shifts. Only o(·) savings are known (Matomäki–Radziwiłł; Ramachandra).
- **Explicit formula.** LM reads avg_u Σ_ρ |W̃(ρ)|²/|L′(ρ)|² ≪ U^ε plus DH-type control, i.e. weak Mertens for the family.
  Under GRH it holds.
- **Quadratic characters.** Unconditionally, nothing beyond the conductor improves on copies by squares (slope 1/2 over
  Heath-Brown's sharp quadratic sieve). Under GRH the upper bound is immediate, but the asymptotic is known only below the
  conductor (Gao–Zhao, arXiv 2212.09241). For the sextic family over Q(ω), Gao–Zhao (arXiv 2201.01885) give moments and
  one-level density under GRH, and Blomer–Goldmakher–Louvel give a large sieve for n-th order Hecke characters. Neither
  reaches beyond the conductor for μ.

## What LM asks of the zeros (HEURISTIC, DH_SEXTIC §62)

With |L′(ρ)| ≍ 1, LM_α holds iff the zero counts satisfy R(δ′) ≤ 1 − δ′ − (s−1)(δ′ − α), i.e. strictly sub-DH above α. In a
DH-saturated world the slope would be κ ≈ 0.75. So LM is morally stronger than the DH it yields, and an unconditional proof
must use cancellation or the size of 1/L′(ρ), not just zero counts. Single exceptional rows cannot block LM for s ≤ 1.334,
thanks to the proved strip β* ≤ 0.8748389. This is a real difference from the ζ-analogue, where one zero near Re s = 1
already blocks power-saving variance.

## Evidence (NUMERICAL)

mean|M_u|²/diag stays at 0.88–0.96 at every tested scale:

- U = 10³ (all 3636 rows) and U = 3·10³, up to cols/rows = 76, i.e. 76 times past the effective conductor;
- U up to 10⁶ in the rows-dominant range;
- a random multiplicative control gives the same plateau;
- the squarefree-indicator control collapses (0.14), as smooth sums should beyond the conductor.

The plateau constant 0.94 is exactly the coprimality factor Π_{p∤6}(1 − 1/(N(p)(N(p)+1))) = 0.9399. Measured against each
row's exact coprime diagonal (all rows, no sampling), MS/truediag = 0.95–1.01 at U = 10³ and 0.986–1.006 at U = 3·10³
for s ∈ [0.8, 2.0]. At U = 3·10³, s = 2.0 is 57 times past the effective conductor, and the off-diagonal there is +0.5% of the
diagonal, within the ±1% row-fluctuation scale (§57). At larger U:
- U = 3·10⁴: 1.0013 ± 0.0034 at the effective conductor (all 108 606 rows), 0.9997 ± 0.0064 at 3.3× and 1.0016 ± 0.0091 at
  11× past it;
- U = 10⁵: 0.986 ± 0.006 at the effective conductor and 1.004 ± 0.007 at 3.4× past it;
- U = 10⁴, s = 1.95: 0.964 ± 0.025 at 121× past the effective conductor (2000 sampled rows).
The random multiplicative control behaves the same.

This is what Proposition C′ predicts under GRH. No finite computation can distinguish LM from GRH-on-average.

## The sharpest form of the problem (17 Oct, DH_SEXTIC §69–§79)

**Equivalent dual statement (PROVED reduction).** Put K = D²/U, which lies in (D, D^{1.5}] for s ∈ (1, 3/2], and
c_n = ε(n)(n/|n|)conj(g̃₃(n)) = μ(n)g̃₆(n). LM at length D holds iff

  Σ_{N(h)≤K} |Σ_{n~D} c_n χ̄_n(h) W(N(n)/D)|² = K·Σ|c_n|²·(local factor) + O(U^ε·D²).

This is a sextic large sieve with an exact leading term, ε-free on K, for one explicit sequence in the rows-dominant range.
- **The coefficients are automorphic up to a periodic sign (NUMERICAL, §77).** On 3522 squarefree n, ε(n) depends only on n mod 4, so
  c_n = ε(n mod 4)·(n/|n|)·conj(g̃₃(n)) exactly. Möbius is absorbed in the relation between the sextic and cubic Gauss sums,
  μ(n) = ε(n)(n/|n|)conj(g̃₃(n)g̃₆(n)). The one non-automorphic feature is the row twist: χ̄_n(h) = (h/n)₂·(h/n)₃, whose cubic part
  makes conj(g̃₃(h, n)) and whose quadratic part (h/n)₂ does not fit any Eisenstein series. (By contrast, the automorphic sextic Kubota
  series Σ g̃₆(h, n)N(n)^{−s} is the cubic object with an extra μ.)
- **Measured directly (§76).** With complete squarefree composite columns, Σ_{N(h)≤K}|Σ c_n χ̄_n(h)|²/(rows·‖c‖²) stays at the
  coprimality constant 0.94 from K ≈ D up to K = 32D, like a random sequence. The unweighted g̃₆ climbs to 1.09.
- **No bias to fight.** The sequence carries no Patterson theta bias (angular twist ℓ = 1; NUMERICAL: Möbius-weighted cubic Gauss
  sums stay at O(√X) up to X = 2·10⁵, while Σ g̃ ≈ 0.355 X^{5/6}).
- **Large sieves do not reach it.** Heath-Brown's and de Faveri's large sieves (arXiv 2610.04045) carry (KD)^{2/3} or
  K^{5/6}D^{1/3} terms, attained for biased sequences (Dunn–Radziwiłł, de Faveri–Dunn–Hoffstein), and an ε on the leading term,
  which is fatal here.
- **Pointwise cancellation is insufficient.** Square-root cancellation of each inner sum gives only the trivial bound.
- **No Rankin–Selberg unfolding.** The quadratic factor (h/n)₂ changes with the row.

**Smallest record-moving target.** Any η > 0 in Σ*_u|M_u(U^s)|² ≪ U^{1+(5/6)(s−1)−η+ε}U^s for s ∈ [1.10, 1.25]. On the record pipeline
dθ/dη ≈ −0.044 (η = 0.001 already gives −4.4·10⁻⁵, more than the whole remaining affine game). Certified in exact rationals:
η = 0.001 ⟹ θ ≤ 0.87480 (DH_SEXTIC §73).

**What cannot work.**
- Positivity alone: copies are the positivity shadow of the exact scale identity A = Σ_k F^{(k)}, so two-sided asymptotics are
  needed (§70).
- Zero-side bootstrapping: the explicit formula fed with the record's own counts returns exactly the copies exponent (§75).
- General large sieves: slope 5/6 = 1 − 1/n is the general-sequence optimum (de Faveri's mechanism is the copies mechanism).
- Convexity in the dual (§79): taking absolute values on Re w = 1/2 and using the balanced functional equation beats copies on the direct
  half (U^{1/2}D^{3/2}). The root-number half, however, is the Möbius family over small rows, weighted by N(h)^{−1/2}. It needs pointwise
  |M_h(D)| ≪ D^{0.773} for bounded h, against the proved strip D^{0.8748}.
- First moments (§77): even square-root cancellation in L¹ gives a worse count than copies at the binding class (δ = 0.467 > 0.232).
- Thick analogue (§77): for the ≍ Q² primitive Dirichlet characters of conductor ≍ Q, the same statement is a power-saving
  Barban–Davenport–Halberstam/Montgomery–Hooley bound below the square-root barrier, which is open. The sextic family is special only
  because Poisson turns μ into Gauss sums (there the dual of μ is a Chowla sum).

## Possible lines of attack (none carried out; DH_SEXTIC §54–§57)

1. **Near the conductor.** Prove LM for D ∈ [U, U^{1+c}] only. This is worth θ ≈ 0.8732 at c = 0.05, and as a power-saving
   asymptotic it would also improve the slope globally via copies plus subtraction of scales. The paper's canonical
   recursion (Lemma 17.2) is an ε-lossy upper bound, so it would need main-term extraction.
2. **Cubic dual.** By Hasse–Davenport, B_h(D) = Σ_n ε(n)(n/|n|)conj(g̃3(h, n))(h/n)_2 W(N(n)/D). This is a twisted cubic
   Kubota coefficient, and the dual mean square is a WMDS-type double series of cubic Gauss sums with quadratic twists.
   One needs that series past its convexity range, with power-saving error. The quartic analogue stopped at convexity.
3. **Asymptotic large sieve (Conrey–Iwaniec–Soundararajan type).** For Möbius coefficients this needs Chowla-type
   shifted-convolution input with power saving.
4. **A different long witness (§62).** LM over all rows inherits its slope from the highest occupied bins. Any long witness
   whose mean square is insensitive to zeros above the target bin would avoid the top-down sub-DH requirement. Restricting
   LM to a bin's rows, however, is equivalent to DH itself, so the witness must still allow Poisson or a large sieve over
   all rows.
