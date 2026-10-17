# θ = 13/15 via the sextic architecture: what it needs, what is established (16 Oct 2026)

DaysSky with Claude assistance. Labels: PROVED (rigorous within the cited framework), MODEL (exponent algebra on the
record's constraint system), NUMERICAL, CITED. Details and scripts: `sextic/DH_SEXTIC.md` §1–§73, `sextic/*.py`;
consolidated note with statements and proofs: `sextic/DH_NOTE.md`.

## 00. Status after the 22:28Z–06:28Z session (details: SESSION_2228.md, DH_SEXTIC §65–§73)
- **Frontier now 0.874836 (Gao, upstream PR #62, merged 16 Oct 18:30Z).** It contains our detector-affine and timbalud's radical count,
  plus a stronger Theorem 2.2 region (λ = 3.9586, B = 0.987047, q = 0.6712, width-coupled HIGH payment). So our local 0.8748389 is
  superseded. No unconditional input that moves 0.874836 by more than about 10⁻⁷ was found (§74).
  - The only unconditional lever found is the free threshold q0 in Cycle61 Thm 2.2. It is worth 7·10⁻⁷ once timbalud's ledger is
    accounted for (§68), and it is not pursued.
- **Exact prices at the record point (MODEL, the record's own pipeline, §71).**
  - A U^η saving over sixth-power copies for the long Möbius mean square on s ∈ [1.10, 1.25] gives θ = 0.8748391 − 0.044η.
    Concretely: η = 0.001 → 0.8747953, η = 0.003 → 0.8747080, η = 0.01 → 0.8744044.
  - Inverse slot cost c_I = 1.95 / 1.9 → 0.8747323 / 0.8746166.
- **Conditional records, certified exactly (§73).** U^{0.001} over copies on [1.10, 1.25] ⟹ θ ≤ 0.87480. c_I = 39/20 ⟹ θ ≤ 0.87474.
  All three row inequalities hold in exact rational arithmetic, and all side conditions pass.
- **The LM dual exactly (§69, PROVED reduction).** LM is a sharp sextic large sieve, with exact leading term, for one explicit
  bias-free sequence ε(n)(n/|n|)conj(g̃₃(n)) in the rows-dominant range. Pointwise square-root cancellation is provably
  insufficient. Rankin–Selberg fails because of the row-dependent quadratic twist.
  - NUMERICAL: Möbius-weighted cubic Gauss sums show no Patterson bias up to X = 2·10⁵.
- **New literature (§72).** de Faveri's optimal n-th order large sieve (arXiv 2610.04045) does not touch the record's inputs:
  - Part II uses only quadratic and cubic sieves, and the cubic one is already optimal.
  - On LM it reproduces the copies slope 5/6 with a worse intercept.
  - It confirms that slope 5/6 is optimal for general sequences, so only Möbius-specific arithmetic can improve it.

## 0. Status after the 8-hour session (14:37–22:37Z)
- **Conditional chain (MODEL + PROVED pieces).** LM (Σ_u|M_u(U^s)|² ≪ U^{1+ε+s} for s ∈ (1, 3/2]) ⟹ DH in every bin
  ⟹ θ = 13/15 (Conditional Theorem D). Without LM, every other known or plausible input stops at or above ≈ 0.869:
  - natural slot costs: 0.8698;
  - natural slot costs with no slot-supply cap: 0.8692.
- **Why LM is the crux (PROVED).** It is a fixed point of every available method:
  - toolkit closure (Theorem A);
  - exact self-duality of Gram matrices and of the long mean square (Theorems B, C);
  - no standard combinatorial identity avoids Möbius-heavy pieces (Theorem D);
  - the convolution identities are satisfied by an explicit adversary.
  LM is the family's weak Mertens problem (explicit formula). Its ζ-analogue, power-saving Möbius variance in short
  intervals, is open.
- **Evidence (NUMERICAL).** LM shows no off-diagonal mass from rows-dominant sizes up to 76× past the effective conductor
  (U = 10³, 3·10³, 10⁴; §45, §57). The plateau is exactly the coprimality constant 0.9399. Exact all-row off-diagonals near
  s = 1 shrink like U^{−1/2} (3% → 1% → 0.6% for U = 10³ → 3·10⁴). This is consistent with LM and with GRH-on-average,
  which no finite test separates. Slot moments show no excess beyond the proven
  ranges (U ≤ 10⁹ at r* = 0.745), including for Gauss-correlated slots.
- **Correction (16 Oct, DH_SEXTIC §50–§53):** 13/15 also needs the paper's bin floor a0 = 51/100 lowered towards 1/2
  (with the floor at 0.51, DH gives 0.86979), and the principal Euler region extended slightly below 0.87. Both were
  checked against the paper's source and are margin changes. DH is needed at the top row scales. Below them the
  existing counts suffice, though not Liu's cruder count (§51).
- **Conditional theorem fully margin-checked (§61, DH_THEOREM.md).** Every margin of the paper's §20 was re-evaluated in
  exact arithmetic at the 13/15 geometry; all hold. DH is needed only at row scales d ∈ [2/5, h+ζ], since below 2/5 the
  paper's own count works. Two region edges move (D1 to 1/2, D2 to 0.86), and the small rows take ε0 = 1/3 in D1(ε0). Each is checked against the
  exact local Euler exponents. The same code reproduces the paper's own margin numbers exactly.
- **Smallest high-value inputs (§55–§56).** LM only on s ∈ [1, 1.05] already gives 0.8732, and on [1, 1.16] gives 0.8695.
  By Hasse–Davenport (checked numerically), the dual coefficients of LM are twisted cubic Gauss sums,
  μ(n)g̃6(n) = ε(n)(n/|n|)conj(g̃3(n)). So the open input is an asymptotic sextic large sieve for weight-1 cubic Gauss sums
  just past the self-dual point. Even slope 0.7 instead of 5/6 on [1, 1.05] moves the model to 0.87467 (§55). DH windows
  alone are weaker: DH above σ = 0.75 is worth nothing (§59).
- **Literature (§48).** Xu–Sang (0.8748867) and Liu (0.874957), October 2026 preprints, add no DH/LM input. Our 0.8748389 is
  below both.
- **Record unchanged:** 0.8748389 (Lean-verified, commit 9acccf3). No commits or pushes this session.

## 1. The target
Conditional Theorem D (PROVED relative to the 7/8 architecture). A density-hypothesis count U^{1+λ−δ} in every bin gives
θ = 11/12 − (1−2λ)/(4(5+6λ)). λ = 0 gives θ = 13/15, which is also the architecture's floor (Part II generic-row
ceiling ℓ < 1/5). At θ ≈ 0.870 the needed count exponent is within 0.07 of DH on all of σ ∈ [0.65, 0.85]. So 13/15
means DH for the sextic family {ψ_u : u ~ U} on the whole window σ ∈ (1/2, 0.875], and the paper's bin floor
a0 = 51/100 lowered towards 1/2. With the paper's floor, DH on [0.51, 0.875] gives θ = 0.86979 (DH_SEXTIC §50). The precise statement
with a proof outline is Conditional Theorem D′ in DH_NOTE §1.

## 2. Where the record stands (MODEL, sextic/bindpt.py)
Record 0.8748389, binding class δ = 0.466 (σ = 0.733), p = ½. Count R = 0.5935 against DH 0.534, a three-way tie:
- Möbius witness at length U^{0.745} (marked inverse moment, slot cost 2): 1 − δr − q(1−r)
- plain witness at U^{0.417} (plain 4th moment, slot cost 6κ / Cycle61 Thm 2.2)
- long Möbius witness at U^{1.162} (sixth-power copies, slope 5/6)

Price list (θ in the Cycle61 grid model, baseline 0.87492):

| input | θ |
|---|---|
| long slope α′ = 0.80 / 0.70 / 0.50 / 0 | 0.87467 / 0.87372 / 0.87055 / 0.86667 (= 13/15) |
| Möbius count −0.01δ / −0.05δ | 0.87468 / 0.87333 |
| plain count −0.01δ / −0.05δ | 0.87481 / 0.87439 |
| inverse slot cost c_I = 1.9 / 1.8 / 1 | 0.874714 / 0.874445 / 0.87289 |
| full sharp sextic sieve (a0 = 1, c_I = c_P = 1) | 0.8696 |

## 3. Established this round
1. **Benchmark (PROVED comparison).** At every σ ∈ [0.70, 0.80] the family's count with slots already beats the ζ
   analogue of Guth–Maynard (δ = .466: 0.5935 vs 0.6162). Its near-line slope (0.82) beats Conrey's mollifier slope
   (4/7). DH for this family in the window is a stronger statement, relative to family size, than every known ζ
   zero-density exponent at the same σ (Bourgain's DH range for ζ starts at σ = 25/32).
2. **Toolkit closure (PROVED by explicit adversary, sextic/closure.py).** With all available inputs and the detector
   for every t ∈ [1, 3/2], the count equals R(δ,p) exactly. Inputs: marked inverse, plain lemma, copies for any product,
   M4 via copies, marked Möbius×Möbius, theta-reflection MM′, pure slots, pointwise sixth moment, and all Hölder
   combinations. The witness row has Möbius saturated at {r*, t*} and plain at t* − r*. It is checked at
   δ = .3, .4, .466, .55, .65.
3. **The long slope 5/6 is an exact fixed point (PROVED).** Poisson in u with factorised Gauss sums gives an exact
   involution MS(U;D) ↔ dual Möbius–Gauss sums at K = D²/U, with coefficient exactly 1. Copies are the only contraction,
   which gives (1+5s)/6. Any α′ < 5/6 needs a third ingredient (spectral/metaplectic or GRH-type).
4. **Guth–Maynard and Chen–Gupta–Li do not transfer (PROVED).** τ(ψ_{u1}ψ̄_{u2}) = ψ_{u1}(u2)ψ̄_{u2}(u1)τ(ψ_{u1})τ(ψ̄_{u2})
   makes the cubic trace exactly self-dual, M_BM_B* = N·I + (N/U)D(M′_BM′_B* − H·I)D*. No free residue or archimedean
   variable remains to produce their linear constraint.
5. **Single-t counting is lossless (PROVED).** Using every t at once is equivalent to a covering-chain condition, and its
   value equals min_t max(short, long).
6. **Large sieve side (CITED + PROVED reading).** Dunn–Radziwiłł and de Faveri–Dunn–Hoffstein (Jul 2026): the n-th order
   sieve is not orthogonal; conjecturally Ξ_n = A + B + A^{1−1/n}B^{2/n} + A^{2/n}B^{1−1/n}. For sextic rows/columns this is
   sharp only to length U^{1/2} (worst case) and U^{5/6} (positive sequences). Pure prime-slot counts are closed even
   conditionally (they need U^{0.87}). A hypothetical full sharp sieve still stops at 0.8696.
7. **Other routes evaluated (MODEL), none below R:**
   - Jensen/Cauchy dichotomy at the zero gives 0.628 at δ = .466.
   - Soundararajan-type central bounds reproduce |L(½)|² ≤ U^{δ/2}.
   - Explicit-formula/Turán routes fail on uncontrolled residues 1/L′(ρ_j).
   - Heath-Brown-identity decompositions only add worse witnesses.
   - Restricted-Gram numerics cannot reach the DH thresholds at computable U.
8. **Literature to Oct 2026 (CITED).** Liu (arXiv 2610.12234, 0.874957) adds nothing beyond the record's κ-dependent
   plain cost. Moment results (DdFDS 2410.03048; Diaconu–Ion–Paşol–Popa 2607.27131, which needs μ_12, not available
   over Q(ω)) have short twists only. No density-hypothesis result exists for any thin family (quadratic, cubic, sextic).

## 4. Minimal new inputs, any one of which moves θ
| input | needed | known | nature |
|---|---|---|---|
| LM: Σ_u\|M_u(U^s)\|² ≪ U^{1+α′(s−1)}, s ∈ (1, 3/2] | α′ < 5/6 (α′ = 0 gives 13/15) | 5/6 | Möbius variance beyond the conductor (single-modulus analogue is also not known for Dirichlet characters) |
| inverse slot cost | c_I < 2 | 2 | slot primes inside the theta-reflection conductor; needs a Type II bound for sextic Gauss sums with one very short variable |
| Möbius large values at U^{0.745} | count < 1 − δr − q(1−r) | second moment | restricted Gram s1(M_B)² ≤ U^{0.94} for \|B\| = U^{0.59} |
| plain 6th moment at U^{0.417} | loss < 0.17 | 0.375 (pointwise) | sixth moment of sextic L-values |
| numerator 6th moment (all rows) | loss < 0.20 | 0.375 | same |

## 5. Next concrete attacks
(a) Type II sums Σ_{p~P, k~K} α_p β_k g̃(pk)χ(·) for sextic Gauss sums with P = U^{0.04}: square-root cancellation in p
    uniformly in the k-coupling would lower c_I.
(b) Metaplectic Kuznetsov input for Möbius-twisted Gauss-sum series Σ_n μ(n)g̃(h,n)W(n/D): any power saving in the
    smoothed fluctuation Σ_h ŵ(h/K)(|B_h|² − D) at K = D^{1+η} lowers α′. η > (s−1)/(6(2s−1)) ≈ 0.02 already beats copies
    at s = 1.162.

## 6. Added in the 8-hour session (DH_SEXTIC.md §24–§37)
1. **Decompositions are exhausted (PROVED, exponent algebra).**
   - Vaughan/Heath-Brown identities plus the in-bin pointwise caps (plain factors > 1/2 cannot be saturated) give only
     descending witness chains. These are exactly the detector's covering chains, so R is unchanged (§24, §27).
   - Möbius-heavy pieces (μ_{≤y})^{∗k} occur in every standard identity for μ, so LM cannot be proved by such
     decompositions (§30a).
   - The convolution identities μ∗1 = δ and μ·log = −μ∗Λ are satisfied by an explicit extension of the adversary.
2. **Where r < 1 enters OpenAI's Möbius second moment (PROVED reading).** The canonical invariant F − M − Q − z0 > 0
   becomes m − r − 2z > 0, and the puncture cancels. At r > m the child is rows-dominant and only the diagonal bound
   survives (§26). Extending the raw moment to D = H^{1+c0} would give slope 5/(6(1+c0)) and is the same third
   ingredient as §22.
3. **LM numerics (NUMERICAL).** mean|M_u(U^s)|²/diag = 0.93–0.95, flat for s ∈ [0.5, 1.5] at U = 3·10³ and 10⁴,
   matching the random multiplicative control. Effective columns/rows only reaches ≈ 1–2, so this is modest support (§44).
4. **Slot costs look like artefacts (NUMERICAL, §31).** The marked moment Σ|M·P|² and the plain fourth moment Σ|S²P|²
   show no excess (ratio 1.000 ± 0.003 at U = 10⁶) far beyond the proven ranges, for every slot-coefficient type including
   Gauss-correlated ones. The source of c_I = 2 is structural to theta reflection: active slot primes lengthen the dual by 2z_a
   while the saving is capped at z_a. A sharper bilinear sieve alone removes only a secondary small-dual penalty (§32).
5. **Literature (CITED, §29).**
   - Gao–Zhao 2212.09241: the quadratic analogue of LM has an asymptotic only for Y ≪ X^{1−ε}, and only under GRH.
   - Dunn 2403.13151: a Type II large sieve for cubic cusp-form coefficients; not applicable to theta/Gauss sums.
   - Sextic theta coefficients are not explicitly known, which blocks a spectral route to LM.
6. **Sextic Möbius–Gauss sums have no visible bias (NUMERICAL, §30c).** Σμ(n)g̃(n) changes sign and stays at size X^{1/2}
   up to X = 3·10⁴. Twisted multiplicativity was verified by brute force.
7. **Cost-1 slots ⟸ Möbius fourth moment (PROVED reduction, §33).** For one slot of length z < 1/4,
   Σ|MP|² ≤ (diag) + Z(U·Σ|M|⁴)^{1/2}. At r* this beats A_I only if the Möbius fourth-moment loss is λ < 0.114
   (copies give 0.408).
8. **LM = the family's weak Mertens problem (PROVED reformulation, §34).**
   M_u(D) = Σ_ρ W̃(ρ)D^ρ/L′(ρ, ψ_u) + O(D^{−ε}U^{−1/2}). LM needs DH-type zero control plus a family-averaged
   Gonek–Hejhal negative moment. Its ζ-analogue (power-saving variance of μ in short intervals) is open: only
   (log)^{−A} (Ramachandra) and o(H²) (Matomäki–Radziwiłł) are known (§36).
9. **One constant behind α′ = 5/6 and c_P = 6κ (§37).** Sixth-power copies have efficiency 1/6. A replication mechanism
   with efficiency e would give α′ = 1 − e and c_P = 1/e. c_I = 2 has a different source (the BGL cubic term, §32).
10. **θ under combined conjectures (MODEL, §35).** LM at α′ = 0 alone gives 13/15. Natural slot costs alone give ≈ 0.8698.

## 7. Updated ranking of the frontier inputs (θ in the Cycle61 grid model; record-equivalent baseline 0.87492)
| input | kind | θ | status |
|---|---|---|---|
| c_P κ = 1 and c_I = 1 (natural slot capacities) | critical length | 0.86979 | numerically supported (§31, rows-dominant range); proof needs a non-reflection treatment of slots (§32) |
| c_I = 1 | critical length | 0.87289 | same |
| c_P κ = 2 (marked plain lemma) | critical length | 0.87339 | numerically supported |
| LM slope α′ = 0 | beyond conductor | 0.86667 = 13/15 | numerically supported, modest range (§28, §44); open even for quadratic characters |
| MLV / MM′ / plain sixth moment | beyond conductor | 13/15 (jointly) | open |
