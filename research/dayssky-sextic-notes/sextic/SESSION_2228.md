# Session 16 Oct 2026, 22:28Z–17 Oct 06:28Z: frontier levers at the record point
**Abstract (06:09Z).**
- **Frontier.** Unconditional: still 0.874836 (Gao, PR #62).
- **New certified conditional bounds (exact rationals, tight, cross-checked).** At the strength of a sextic analogue of each classical level-of-distribution theorem: BFI 0.87425, Maynard 0.8740, Lichtman 0.87382; with LM slope 4/5 as well: 0.87405 / 0.8738 / 0.87362. Plus the earlier ladder down to 0.8692 and 13/15.
- **Binding inputs, located exactly.** Each is checked numerically at the binding point:
  - inverse slot cost: three-variable cancellation, single inequality in PLAN_RD.md; ρ = 1.00–1.06 at r = 0.745, U = 10⁹;
  - plain slot cost: thin-family fourth moment; loss exponent 0.05–0.07 at U = 3·10⁶ (0.08–0.23 at 10⁶) against the 0.36 needed;
  - Möbius fourth moment: loss 0.40 against the 0.114 needed.
- **Unconditional side results.**
  - Slot cost 5 at κ = 1 (§92).
  - The theta measure has squarefree local densities 1/(N(p)+1) (§94).
  - Thresholds where the problem becomes self-improving: θ₀ < 5/6, 0.755 and 0.575 (§96, §99, §100).
- **Closed routes (§92–§102).**
  - averaging over slot primes;
  - dispersion on |M_u|²;
  - density loops;
  - feeding the proved region back in as pointwise bounds;
  - Hölder;
  - a row-dependent κ;
  - summing over slot primes first;
  - Poisson in h;
  - changing the family.


DaysSky, with Claude assistance. Labels as in DH_NOTE.md. No commits, no pushes.

**Frontier update (§74).** Upstream merged Gao's directional affine at **0.874836** (GTimG/QuasiRiemannTracker PR #62, 16 Oct
18:30Z). It builds on our detector-affine and timbalud's radical count, so our unsubmitted local 0.8748389 (9acccf3) is
superseded. The tables below are measured from our pipeline (record 0.8748391); the slopes carry over to Gao's point.
Details are in DH_SEXTIC.md §65–§75.

## The record point, priced exactly (MODEL, the record's own pipeline)
All numbers come from the record pipeline (detector-affine, timbalud B = B0 + 7·10⁻⁴, all three row constraints and the
geometry checks), not the older grid model. The record value 0.8748391 is reproduced.

| open input | θ | change |
|---|---|---|
| none (record) | 0.8748391 | — |
| long Möbius mean square U^{0.001} below sixth-power copies on s ∈ [1.10, 1.25] | 0.8747953 | −4.4·10⁻⁵ |
| same, U^{0.003} | 0.8747080 | −1.3·10⁻⁴ |
| same, U^{0.01} | 0.8744044 | −4.4·10⁻⁴ |
| U^{0.03} on [1.05, 1.30] | 0.8735620 | −1.3·10⁻³ |
| inverse slot cost c_I = 1.95 / 1.9 / 1.8 | 0.8747323 / 0.8746166 / ≤ 0.8743642 | −1.1·10⁻⁴ / −2.2·10⁻⁴ / ≤ −4.7·10⁻⁴ |
| affine constants of Cycle61 Thm 2.2 (best q0 ≈ 0.67, with timbalud's ledger) | 0.8748384 | −7·10⁻⁷ |

So every positive power saving over copies in that window is a new record, with dθ/dη ≈ −0.044. This is the cheapest open
input.

**Measured from Gao's frontier (§74).** With λ = 19793/5000 and B = 987047/10⁶, my pipeline gives 0.87483584 (certified upstream
0.874836). From there: U^{0.001} / U^{0.003} / U^{0.01} LM savings give 0.8747923 / 0.8747054 / 0.8744030, and c_I = 1.95 gives
0.8747293. The remaining affine game is capped at about 3.4·10⁻⁵ (window cap z ≤ 0.04), and Gao's HIGH comparisons are tight at
both ends of the cell range (slack 5·10⁻⁶).

**Conditional records, certified in exact rationals (§73).** U^{0.001} below copies on [11/10, 5/4] ⟹ θ ≤ 0.87480. Inverse slot cost
39/20 ⟹ θ ≤ 0.87474. All three row inequalities hold and every side condition passes. The η = 0 sanity run reproduces the record
pipeline.

## What was settled (PROVED readings and reductions)
- **The z-contour is not a lever (§68).** In paper (20.4) the coefficient of Re z is h − d, which vanishes at the binding top
  scale. So 17/50 is already optimal.
- **Theorem 2.2 is a one-parameter family in q0 (§68).** B0(q0) = (1 − ρ0 + c)/(q0 + c) and Λ(q0) = (q0 + 9c/2)/(q0 + c), with
  c = K_a(5 − 6ρ0) and K_a = (2q0/9)/(1 − 4q0/3). This reproduces 14200/14361 and 19312/4787 at 17/25. Timbalud's ledger ties its
  saving to κ = q0 − 3ρ0/2, so the whole lever is worth about 7·10⁻⁷ (not pursued).
- **Lemma 14.2 is sharp in the binding regime (§68).** Both sieves sit at their diagonals, so c_I = 2 is lost upstream: the mark
  q | A enters the reflection conductor as a cubic character of conductor q (Cor. 14.1).
- **Sixth moments are beyond multiple Dirichlet series (§68).** Σ_u Π_{i≤r} L(s_i, ψ^{±1}) u^{−w} has diagram K_{1,r}, which is
  finite type for r ≤ 3, affine D̃₄ for r = 4 (the paper's critical fourth moment), and indefinite for r ≥ 5. Both the plain and
  the numerator sixth moments are r = 6.
- **The LM dual exactly (§69).** LM at length N is equivalent to a sharp sextic large sieve with exact leading term K·‖c‖²
  and error N^{1+ε}‖c‖², for the single explicit sequence c_n = ε(n)(n/|n|)conj(g̃₃(n)), in the rows-dominant range
  K = N²/U ∈ (N, N^{1.5}]. Heath-Brown's (KN)^{2/3} term is attained only by biased Gauss-twisted sequences (Dunn–Radziwiłł
  Thm 1.4), and c_n is bias-free. Two further facts:
  - Pointwise square-root cancellation of the dual sums is provably insufficient: it only gives the trivial N².
  - A Rankin–Selberg unfolding on the cubic cover fails, because the quadratic factor (h/n)₂ changes with the row.
- **Copies are the positivity shadow of an exact identity (§70, cf. §54).** A(U′; N) = Σ_k F^{(k)}(U′/N(k)⁶; N). Beating copies by
  any power needs two-sided asymptotics at the intermediate supercritical scales, since a lower bound there is again an
  asymptotic large sieve.

- **The zero-side bootstrap is a fixed point (§75).** R(δ) + t*δ = 1 + (5/6)(t* − 1) exactly, so the explicit formula fed with the
  record's own zero counts returns the copies exponent. Gonek–Hejhal-type control of 1/L′ cannot help without a sub-R count.
- **de Faveri's optimal n-th order sieve (arXiv 2610.04045, §72) does not touch the record.**
  - Part II uses only the quadratic and cubic sieves, and the cubic one is already optimal.
  - The sextic sieve enters only Part I.
  - In the LM range it reproduces the copies slope 5/6 with a worse intercept, which confirms that 5/6 is the general-sequence optimum.

- **Where LM sits (§77).** Its thick-family analogue is the power-saving Barban–Davenport–Halberstam/Montgomery–Hooley variance of μ
  below the square-root barrier, which is open unconditionally. So is the ζ analogue (Matomäki–Radziwiłł). The sextic family differs
  only in that Poisson turns μ into Gauss sums. A perfect L¹ input would also lose to copies at the binding class (δ = 0.467 > 0.232).


## Second half (00:15Z–): the four inputs (DH_SEXTIC §76–§84, FOUR_INPUTS.md, PLAN_INPUT2.md)
- **Every substantial lever reduces to one of four inputs**: LM; the plain slot cost c_P; the inverse slot cost c_I; the sixth moment of the
  HIGH-side central block. Each one's mechanism was traced in the papers' proofs (§78, §80, §81). Two roots sit underneath them all: the copy efficiency 1/6,
  which is forced because the cubic theta series' Poisson frequencies are u·a⁶, and theta reflection squaring a mark's conductor.
- **Input 2, worked to the bottom (§83–§84).**
  - The paper's 6 is 2 + 4: the 4 is the index of the dual sublattice h·m⁴.
  - Cost 2 is equivalent to a fourth moment over the subfamilies u ≡ 0 mod m² at natural size. These are thin by only U^{0.04–0.08}.
  - Beating 6κ needs only a saving N(mm′)^{−θ′} with θ′ > 0.32 in the twisted fourth moments.
  - Those moments have no main term in either representation and are self-dual.
  - Measured: θ_eff ≥ 0.6.
  - Priced on the record pipeline: slot weight 2 gives 0.873382 (−1.45·10⁻³), weight 1 gives 0.870524 (−4.3·10⁻³; re-checked §89).
  - The Dirichlet analogue is proved by Kuznetsov theory. The plan in PLAN_INPUT2.md would use Louvel's and Dunn's metaplectic tools on SL₂(Z[ω]).
- **Input 4 has evidence (§82).** Excluding the principal row, the central block's sixth moment equals a random multiplicative model up to U = 3·10⁶. The
  best unconditional bound is λ = 3/8; the frontier needs λ < 0.20.
- **Negative results.** Convexity in the dual (§79); first moments (§77); Hecke relations at T_{q³} for c_I; the shifts h ↦ hp³ as amplifiers. Theorem 2.2's
  window does not matter; its B₀ < 1 is LM again (§83 addendum 4). The quartic route gives 7/8 even with full (L³). Gao's own certificate has about
  2·10⁻⁷ slack (not substantial, not pursued).
- **Corrections.** The dual coefficients are exactly ε(n mod 4)(n/|n|)conj(g̃3(n)), with no Möbius sign left (§77, which corrects a first draft).

## New numerics
- **The dual sieve measured directly (§76; dualsieve/ds2.py).** With squarefree composite columns and direct Gauss sums, the LM dual
  sequence μg6 sits on the coprimality plateau 0.94 from K ≈ P to K = 32P (P = 1000 and 3000), the same as a random control. The
  unweighted g6 rises to 1.09 (theta bias). So Möbius removes exactly the obstruction that makes the general sextic sieve non-sharp.
- **The dual coefficients are periodic × Hecke × cubic (§77; g2per/g2per.py).** On 3522 squarefree n, both g̃2(n) and
  ε(n) = μ(n)g̃6(n)/((n/|n|)conj(g̃3(n))) depend only on n mod 4. So c_n = ε(n mod 4)(n/|n|)conj(g̃3(n)) exactly, with no Möbius
  factor left anywhere in the coefficients. The only non-automorphic feature of the dual is the row twist (h/n)₂.
- **Möbius-weighted cubic Gauss sums show no Patterson bias (§69; mugauss/, X ≤ 2·10⁵, 70 554 squarefree c).**
  - Σ g̃(c) ≈ 0.355·X^{5/6}.
  - Σ_ϖ g̃(ϖ) ≈ 1.1·X^{5/6}/log X, consistent with Dunn–Radziwiłł's constant ≈ 1.006.
  - Σ μ(c)g̃(c) stays within 0.63·√X on every scale.
  - Hasse's formula gives μ(c) = g̃(c)³·conj(c/|c|) exactly.

## Literature checked this session
Dunn, "Metaplectic cusp forms and the large sieve" (arXiv 2403.13151): bilinear sieve for cusp-form coefficients, not
theta. Ma (arXiv 2606.26440): function-field Gauss-sum bias bounds. Dunn–Radziwiłł (arXiv 2109.07463): Patterson bias and
sharpness of the cubic sieve. de Faveri, "Optimal large sieve for fixed order characters" (arXiv 2610.04045, 2 Oct 2026); see §72.
Gao, "A width-sensitive affine moment bound for Dirichlet L-functions" (upstream PR #62, 16 Oct): the new frontier
0.874836. None of these supplies the missing sharp sieve for bias-free sequences.

## Third segment (02:15Z–): LOW side, slot-twist sieve, joint prices (DH_SEXTIC §87–§89)
- **LOW side (§87).** It has the steepest slope found (dθ/dΔ ≈ −0.77), and its b/12 comes entirely from Prop. 15.2's sixth-power exceptional frequencies.
  Beating Cauchy–Schwarz there bounds the η-Möbius sum itself (REPORT §12.4), so it is no independent lever. The HIGH constants 13/75 and 17/50 cancel at d = h.
- **Large sieve over the slot twists (§88, PROVED).** Σ|M_uQ_u|² ≤ Σ|M_u|² + (UΣ|M_u|⁴)^{1/2} for slots up to U^{1/4}. Consequences:
  - slots are free given M4;
  - unconditionally it loses at the binding point (0.708 against 0.592);
  - for a single plain block, slots up to U^{1/4} are free unconditionally (Lemma 18.1 allows 0.11);
  - it does not reach the record's two-copy legs.
- **Joint prices (§89, MODEL).**
  - Geometry, LOW validity, floor, slot supply and the affine window are all enforced.
  - c_I = 1.8 / 1.5 / 1.2: 0.874353 / 0.873204 / 0.872762, saturating at the supply ℓ/h.
  - LM slope 0.80 / 0.75: 0.874589 / 0.874186.
  - Plain weight 3 / 2 / 1: 0.874313 / 0.873382 / 0.870524.
- **Correction (§89).** The quick price scripts of this session lacked those conditions. Only c_I ≤ 1.2 moved materially. The certificates were unaffected.
- **Literature (FOUR_INPUTS.md).**
  - DIPP 2607.27131: second moments; sextic over Q(ω) not covered.
  - de Faveri–Dunn–Hoffstein 2607.07911: non-orthogonality.
  - The sixth moment with λ < 0.20 is beyond ζ's best (loss 1/4).
- **Conditional bounds certified in exact rationals (§90; c61/certG.py, certG2.py; relative to the radical record's proof; all three rows plus side conditions).**
  - LM slope 4/5 ⟹ θ ≤ 0.87462 (RD with c₀ = 1/24);
  - LM slope 3/4 ⟹ θ ≤ 0.87425;
  - LM slope 7/10 ⟹ θ ≤ 0.87368;
  - c_I = 3/2 ⟹ θ ≤ 0.87325;
  - c_I ≤ 6/5 ⟹ θ ≤ 0.8728 (the supply saturates);
  - plain weight 2 ⟹ θ ≤ 0.8734;
  - c_I ≤ 6/5 with slope 7/10 ⟹ θ ≤ 0.8716;
  - all three together ⟹ θ ≤ 0.8708.
  Sanity: the scripts fail just below the model optima. An independent float re-implementation (cp/xcheck.py) agrees at every point.
  The top row uses the paper's old capacity wherever the Theorem 2.2 window fails. Cycle61's reserves (4.6) and (4.7) need only positivity (its (4.8)–(4.10)).
- **The single statement behind inputs 1 and 3 (FOUR_INPUTS.md).** It is the rows-dominant asymptotic (RD_{c₀}) for the canonical family: c₀ = 1/9 is LM slope 3/4, and with marks it gives c_I.
  It is an asymptotic with exact main term, so Lindelöf-on-average bounds do not supply it. Its shape is Bombieri–Friedlander–Iwaniec for well-factorable moduli pp′.

- **Fourth segment (04:30Z–).**
  - Certified: c_I = 3/2 with LM slope 3/4 (the RD_{1/9} representative) ⟹ θ ≤ 0.8728 (h = 0.7957; 95 979 / 39 447 / 31 boxes; certG_rd19.log).
  - Model: full LM 0.86911 and slope 2/5 0.87042, where the bin floor nearly binds (slack 3·10⁻⁴ at slope 0).
  - **§92 (PROVED reduction): plain slot cost 5 without zero-free input.** Average over the slot primes before positivity, using the bounded overlap of the sublattices m⁴O. This improves Lemma 18.1's κ = 1 case
    from 6 to 5. It does not beat 6κ₀ ≈ 4.5 at the record. The addendum shows that combining it with deletion or with live slots cannot give 5κ: the averaging gain needs the zero-slot lemma's
    supercritical range.
  - **§93: c_I = 7/4 ⟹ θ ≤ 0.87425, certified** (h = 0.806, Theorem 2.2 window valid). This is the BFI-strength analogue (level 8/7 of the large-sieve barrier), equivalently the marked RD with c₀ ≈ 0.048.
  - **c_I = 7/4 with LM slope 4/5 ⟹ θ ≤ 0.87405, certified** (h = 0.795; RD at c₀ ≈ 1/24, marked and unmarked).
  - **Maynard strength (x^{3/5}): c_I = 5/3 ⟹ θ ≤ 0.8740; with LM slope 4/5 ⟹ θ ≤ 0.8738, both certified.**
  - **Lichtman strength (x^{66/107}): c_I = 107/66 ⟹ θ ≤ 0.87382; with LM slope 4/5 ⟹ θ ≤ 0.87362, both certified.** All classical-strength rows have sanity failures just below.
  - **§94:** BFI dispersion needs a smooth variable, which in the sextic setting carries the Patterson bias. Type I pieces have c(D)·X^{5/6}, with c(D) = Σ_{N(d)≤D}μ(d)Π(N(p)+1)^{−1} fitting the data, so Möbius kills theta via 1/ζ_K(1 + w); a direct per-prime test gives theta fractions 1/(N(p)+1) to within noise for N(p) ≤ 25. This kills the main term but gives no rate. The missing input is uniform-in-d power-saving error terms for twisted Kubota series.
  - **§95–§96:**
    - Power-saving Möbius orthogonality to theta needs the Type III sums that Dunn–Radziwiłł call self-dual and open.
    - DR's corrected cubic sieve turns the critical cubic sieve for μ into a zero-density question: fewer than H^{2/3} at σ = 5/6 suffices. Heath-Brown densities are exactly at the threshold, and below θ = 5/6 the loop would self-start.
    - Density estimates cannot give RD: zero counts give LM slope (2σ − 1) + E/(s − 1).
    - The detector's density curve has an excess over DH peaking at 0.084 near σ = 0.73.
  - Row 2′ re-run reproduces identical boxes and slacks. Row 5′ float cross-check −1.3·10⁻⁴; its sanity run at 0.8726 fails.

## Three exact targets for the next session (from §92–§100)
1. **Marked RD at c₀ ≈ 0.048 (c_I = 7/4, θ ≤ 0.87425 certified).**
   - Required: orthogonality Σ_h B_{hm⁴}conj(B_{hm′⁴}) = o(diagonal) for distinct slot primes, just past the self-dual point. Upper bounds are already enough there (§98).
   - Classical model: BFI, Maynard, Lichtman.
   - Obstruction: the Type III self-duality of Dunn–Radziwiłł, and the Patterson pole in Type I pieces, which Möbius kills only through 1/ζ_K (§94–§95).
   - Evidence: balanced Möbius Type II sums are at √X (§97); the RD plateau holds up to K = 32D (§76).
2. **Thin-family plain large sieve with loss < Z^{0.36} (plain weight below 6κ; weight 2 gives 0.8734).**
   - Required: Σ_{v≤U/Z²}|S₁S₂(vm²)|² ≪ (U/Z²)Z^{0.36}. This is a pure upper bound. The general sieve loses U^{0.13}; Hölder with the eighth moment fails (§98, §100).
   - Classical model: twisted fourth moments (BFKMM); the twist range is unchecked.
   - Evidence: the thin-family mean of |S₁S₂|² on u ≡ 0 mod m² is at natural size (ratio ≤ 1.12) for N(m) up to U^{0.335} at the binding plain length (§98 addendum, U = 10⁶). At U = 3·10⁶ the ratio is 0.77–1.25 up to N(m) = 199 = U^{0.355}, and the fitted loss exponent drops to 0.05–0.07.
   - Proof route (§98 addendum, 06:20Z): Gauss-expanding the cubic twist reduces input 2 to decorrelating the N(m) additively twisted plain sums across residues b mod m, the same step as marked RD (target 1). Reduction (06:22Z): input 2 follows from the optimal sieve (V + X²) restricted to coefficients with no theta correlation, lossless for Z ≤ U^{0.083}. Numerically the χ_m-twisted theta residue vanishes while the χ_m² residue survives at about 0.63 of the untwisted value (mugauss_t1/twist_4e4.log). The surviving residue follows N(m)^{−1/6}·N(m)/(N(m)+1) to within 2% (twist2_1e5.log), which explains it by twisted multiplicativity. It saves at most N(m)^{1/3} ≤ U^{0.027} against the U^{0.13} middle-term loss, so this route alone does not close input 2.
3. **Möbius fourth moment at r ≈ 0.75 with loss < 0.114.**
   - Copies and zero counts both give about 0.40 (§100). It frees slots (§88) and is what the dispersion route needs (§100).

## Next steps (for a future session)
1. **RD in its smallest useful form** (PLAN_RD.md). c₀ = 1/24 already beats the frontier (0.87462, certified); c₀ ≈ 0.19 is worth −1.2·10⁻³. The proposed route is BFI-type
   dispersion over the well-factorable slot moduli, then metaplectic Kuznetsov on Kubota's Γ(3). The missing pieces are a spectral gap uniform in level (pp′)² and a
   "Möbius orthogonality to theta" statement.
2. **Twisted fourth moments** (input 2, PLAN_INPUT2.md). Plain weight 2 is certified to give 0.8734.
3. **Lower the bin floor 51/100 toward 1/2** in the certificates (§52–§53 show it is a margin). This is needed only for conditional bounds below ≈ 0.869.
4. **Lean.** If any hypothesis gets proved, certG/certG2 already give the exact-rational row certificates. Those would then need to be ported to the tracker's Lean pipeline (as for 9acccf3).

**Next concrete step (06:25Z).** Split the thin-family coefficients as c_nχ_m(n) = (theta-orthogonal part) + (projection onto g̃ restricted to multiples of m, weight N(m)^{−1/6}). The projection
part's large-sieve contribution can be computed exactly through the theta function at level m (Patterson), and the orthogonal part needs the sieve without its middle term. Test numerically
whether the thin-family loss data (mmnum/thin_3e6.log) are explained by the projection alone. If so, input 2 reduces to the "sieve without middle term for theta-orthogonal coefficients" statement.
