# S52. Mean-value Bochner: (L³) ≪ Q^{1.5417+ε} without Chinta–Gunnells (DaysSky with Claude assistance)

Labels as before. Scripts: l3f/hull3.py … hull6.py, fixa.py, fixa6.py, prof*.py, combo6.py, second.py.

## 0. Idea in one paragraph
* **The input S50 needs.** After the Cauchy route, (L³) needs bounds for first moments X(h₀), i.e. for the double series
  Z_{h₀}(s,w) = Σ_c T(c)N(c)^{−s}ψ(h₀c,w) near (Re s, Re w) = (½, ¾) (the theta residue) and (½, ½) (the dual term).
* **The known route.** Chinta–Gunnells (S51) would continue Z to C² by Bochner's tube theorem, starting from the cone of absolute
  convergence and its 12 images under W(G₂). That needs the full CG machinery, and it is conditional.
* **New observation.** Bochner does not need absolute convergence. Any open set where a series representation converges locally
  uniformly can be used, and mean-value theorems give much larger sets:
  * LoA from the quartic large sieve, in the d-variable;
  * S49's (L³)-type mean square, in the c-variable;
  * pointwise convexity on short blocks.
* **What makes it work.** Each chamber's series converges on such mean-value regions. Bochner's tube theorem plus an envelope-of-holomorphy
  sup argument gives continuation with polynomial growth on the convex hull, and the convex envelope of the h₀-exponents.
* **Result.** Bounds for X(h₀) that beat Cauchy–Schwarz (P4-type input), unconditionally, and a bootstrap through S50.

## 1. Normalisation facts (PROVED from DDHL (4.24)–(4.28))
* DDHL's clean functional equation (4.28) is for ψ̂ = G_∞(s,ℓ)·ζ_{Q(i),λ}(4s−3,ℓ)·ψ, not for the Kubota series ψ itself.
  ζ_{Q(i),λ}(s,ℓ) is the Hecke L-function of (m/|m|)^{4ℓ}.
* In our normalisation (centre ½), the full Kubota series therefore carries 1/L(4s−1, α^{4ℓ}) on the left of the centre. For x ∈ (−¼, 0)
  this has poles at L-zeros.
* Consequences:
  * (i) Bochner is applied to F = N(s,w)·(w−½)(w−¾)·E₂(s,w)·Z_{h₀}. Here N is the product over the 6 positive roots β of
    L(1+4β(x,y), α^{4ℓ_β}), with the root twists of S51 §17, and E₂ is the 2-adic factor clearing the A_{ji}-poles on y = 0, ¼.
  * (ii) Z itself may only be evaluated where every β(x,y) ≥ 0. Without fix (a) that means x ≥ 0 and y ≥ 0. With fix (a) (convolution
    T = T*∗λ, S51 §20) the constraint x ≥ 0 can be dropped, keeping 3x+y, 2x+y, x+y, y ≥ 0. Evaluating E₂ at y = δ > 0 instead of on
    y = 0 costs Z^{O(δ)}.

## 2. The 12 chambers (PROVED bookkeeping, l3f/hull3.py)
* Chamber w ∈ W(G₂) has coordinates v′ = M_w v and twist t_w = M_w(1,0) (S51 §17).
* It has an h₀-factor N(h₀)^{φ_w(v)}: each d-reflection at chamber u contributes −(M_u v)_y, from N(h₀c)^{½−w}. c-reflections contribute
  nothing, since the c-node index is d³ and free of h₀.
* In every chamber the c-node is a twisted χ̄-Kubota family at index (d-variable)³, with |ℓ_c| ∈ {1,2} (never 0). The d-node is a
  Kubota family at index h₀·(c-variable), with ℓ_d ∈ {0, ±3}.
* Non-coprime corrections: for squarefree c, the only non-coprime d-terms are d ∋ π² with π ∥ c, of normalised size √N(π). On the c-side,
  the terms c ∋ π^{3j+1}, π^{4i} with π^j ∥ d. Each is a convergent local correction in every region used below (their exponents are
  −1−x−2y, −½−4x, etc., checked to be negative on the regions). STANDARD bookkeeping, not written out.

## 3. Mean-value regions (each PROVED from the stated input; constants and corrections STANDARD)
In chamber coordinates (x′, y′):
* **R_d (LoA, d-variable).** Σ_{c∼C}|ψ(h₀c,½+y′+it)|² ≪ C^{1+ε}N(h₀)^{1/3}(1+|t|)^A for y′ ≥ 0. This follows from the quartic large
  sieve (C + D + (CD)^{2/3}) with D = √(N(h₀)C), plus the theta polar term. For y′ < 0 use the FE of ψ̂ (no 1/ζ).
  Region: x′ > ½ + max(0,−y′), with h₀-exponent b = 1/6 + max(0,−y′).
* **R_d′ (blockwise, l3f/hull6.py).** For blocks c ∼ C, d ∼ D = C^δ, take the better of:
  * LoA: exponent ½ − x′ − δy′;
  * convexity of the c-series ψ_α(d³,·) (|Σ_{c∼C}T(c)χ̄_d(c)| ≪ C^{½}N(d)^{¾}): exponent δ(5/4 − y′) − x′.
  This gives the region x′ > ½ − 0.4y′ for 0 ≤ y′ ≤ 5/4, and x′ > 0 for y′ ≥ 5/4.
* **R_c ((L³)-type, c-variable).** Let Λ(x) be the mean-square exponent of the c-node family at Re s = ½+x over d ∼ D:
  * Λ(0) = λ, the current (L³) exponent;
  * Λ(x) = max(1, 13/8 − 3x) for x > 0, from the S49 large sieve with the AFE;
  * Λ(x) = max(1−6x, 13/8−3x, λ−6x) for x < 0, from the FE of ψ̂.
  Region: y′ > Λ(x′)/2 (for x′ ≤ ½), and y′ > ½ for x′ > ½. b = 0.
* Pull each region back by M_w and add φ_w to b.

## 4. Bochner with growth (PROVED, standard several-complex-variables argument)
Let U be a bounded connected open set where F is holomorphic with |F| ≤ C_U(1+|Im z|)^A N(h₀)^{b(Re z)}.
* **Continuation.** By Bochner's tube theorem, every holomorphic function on T(U) extends to T(conv U).
* **Sup bound.** If |f| ≤ M on T(U), then for |λ| > M the function 1/(f−λ) extends, so f ≠ λ on T(conv U). Hence |f| ≤ M there.
* **Polynomial growth.** Apply the sup bound to f = F·Π_j(z_j+M₀)^{−A}·N(h₀)^{−ℓ(z)}, with ℓ affine and ℓ ≥ b on U. This gives
  |F| ≪ (1+|Im z|)^{2A}N(h₀)^{ℓ(Re z)} on T(conv U). Minimising over ℓ yields the convex envelope of b (LP duality).
* **Vector-valued version.** The same argument works for ℓ²(h₀∼K)-valued F, applying it to each linear functional. Region bounds may
  therefore be taken in RMS over dyadic h₀-ranges, which is all S50 needs. It does not change the numbers below.

## 5. Evaluation and bookkeeping (NUMERICAL LP over sampled regions; samples lie inside the regions, so the hull is under-estimated and the bound is valid)
Write N(h₀) = Z^η and Q = Z^{2/3}. Then:
* Main term: X_main ≪ Q^{3/4}·Z^{½+x}N(h₀)^{b} at (x, ¼) in the hull.
* Dual term: X_dual ≪ Z^{½+x}Q^{½+y}N(h₀)^{b} at (x, y) in the hull with 0 ≤ y < ¼.
* No shift past the pole: same formula with y > ¼.
* S50: s = min over k₀ of max(2 + (⅓−k₀)/4, (1 + max_{η≤k₀}[⅓ + η + (⅓−η)/4 + 2e_X(η)])/2).

| regions used | (L³) exponent, first pass from λ = 13/8 | after bootstrap |
|---|---|---|
| R_d, R_c (identity chamber only) | 1.6124 | 1.6102 |
| + σ_d-image of R_c | 1.5885 | 1.5833 |
| + σ_c-image of R_d | 1.5729 | 1.5700 |
| all 12 chambers | 1.5469 | **1.5417** |
| all 12 chambers + fix (a) | 1.5417 | **1.5365** |

* **Profile (all chambers, fix (a)).**
  * X(h₀ = 1) ≪ Z^{0.917}, against Z^{1.29} needed. Small twists are no longer the problem.
  * The h₀-slope is ≈ 1.29, so twists N(h₀) ≥ Z^{0.236} fall back to the large sieve with the j⁴-multiplicity (H/K₀)^{1/4}.
    That is the only loss left before Q^{3/2}.
* **Why large twists resist.** The optimal hull combinations (l3f/combo6.py) reach y = ¼ only by mixing the b = 0 point (0, λ/2) of
  R_c with d-reflected points, whose FE factor N(h₀c)^{½−w} costs N(h₀)^{|y|}.
* **Second moment.** The analogous hull for the (2,6,6) second-moment series (l3f/second.py) gives nothing beyond the bootstrap. The
  m-reflection turns T(c)g̃(χ̄_c) into μ(c)·(Hecke), i.e. 1/L, so the obstruction is genuine.

## 6. Status
* **Theorem (modulo the STANDARD items marked in §2–§3).** Σ_{N(m)≤Q}|ψ_α(m³,½+it)|² ≪ Q^{1.5417+ε}(1+|t|)^{O(1)}. With fix (a) the
  exponent is 1.5365; that additionally needs S50 run uniformly over c-lengths Y = Z/N(e)⁴.
* **Comparison.** This replaces 13/8 = 1.625 (S49) and needs no Chinta–Gunnells or GRH input. Inputs: DDHL (4.28), the quartic large
  sieve, S49, S50, Bochner.
* **Remaining gap to Q^{3/2}.** Twists Z^{0.24} ≤ N(h₀) ≤ Z^{1/3} need X(h₀) ≪ Z^{1.29−0.375η}. The MDS gives Z^{0.92+1.29η}.
* **Below Q^{3/2}.** The Cauchy route has a hard floor at Q^{3/2}; going below needs a non-Cauchy argument (OPEN).

## 7. All row boxes, not just the main one (01:01Z; l3c/boxes_rest.py, boxes_decay.py)
* **The issue.** S49's box grid has worst box = the main box (13/8), but other boxes are close: 1.6208 for near-main boxes with tiny
  m₂, m₃, m₄ parts, and 1.6000 for, e.g., box (x₁,x₂,x₃,y) = (0.6,0,0,0.1) in Q-units. So §6's 1.5417 is the main-box exponent only.
* **Reduction of any box to the main box — PROVED modulo the twisted input.** Write m = m₁m₂²m₃³m₄⁴ and ν₀ = m₁³m₂²m₃ as in S49's
  Corollary (Möbius over d | rad x). Shift the c-contour of each inner smooth sum to Re s = ½; the 1/L normaliser is fine on Re = ½. Then
  S_box ≤ Z·Σ_{d}N(d)^{−1}Σ_{m∈box}|ψ_α(ν₀d²-index, ½+it)|², integrated against |W̃|. So S_box ≤ Z·X₂X₃X₄·max_κ Σ_{m₁∼X₁}|ψ_α(m₁³κ,½+it)|²,
  with κ = m₂²m₃d² (orientation as in S49).
* **Twisted (L³).** Suppose Σ_{m₁∼X}|ψ_α(m₁³κ,½+it)|² ≪ X^{λ}N(κ)^{θ}. The box is then at most the main-box bound whenever
  2λ ≥ 1+2θ, 3λ ≥ 1+θ and 4λ ≥ 1, using x₁+2x₂+3x₃+4y ≤ 1. That holds for θ ≤ 1 and λ ≥ 1.
  * With S49 alone (AFE length (X³N(κ))^{1/2}) one gets λ = 13/8 and θ = 1/2. This already disposes of every box whose S49 value is below
    1.5417 + slack, including the 1.6000 box: its recursive value is 0.1 + 0.6·1.5417 = 1.025.
  * The near-main boxes (x₁ → 1) need the twisted version of §2–§5, with the c-node index d³κ, and κ-exponent θ ≤ 1. Inside the MDS, κ
    enters only through c-reflection factors N(κ)^{−x′}, so θ is the convex envelope of max(0, −x′) over the hull points used.
    Expected θ ≤ 1, but NOT yet computed.
* **Corrected status line.** (L³) ≪ Q^{1.5417+ε} holds for the main box (modulo §2–§3 STANDARD items). For the full sum it is CONDITIONAL
  on the twisted-κ version of S52 with θ ≤ 1, which is the next item to compute. Without it the full-sum exponent is the max of 1.5417 and
  the near-main S49 boxes, i.e. still 13/8 in the limit x₁ → 1.

## 8. Correction to §7: reflect the boxes instead of reducing them (01:04Z)
* **§7's crude route fails for square parts.** Shift to the centre and use twisted (L³) with κ = m₂²d². The MDS κ-exponent comes out at
  ≈ 0.75 for X (l3f/kappa.py), so ≈ 1.26 for the (L³)_κ sum once the AFE-length factor N(κ)^{1/2} is included. With N(κ) up to X₂⁴, the
  m₂-boxes would then cost 1.5417 + 2.97x₂. Only S49's own box bound (≈ 1.625 − 0.47x₂) caps that, at ≈ 1.614 overall. Not acceptable.
* **Better: the c-FE (with fix (a), i.e. the T* family, so the FE is clean).**
  * A box m = m₁m₂²m₃³m₄⁴ has effective index ν₀ = m₁³m₂²m₃ and conductor N(ν₀) ≤ Z²·X₂^{−4}X₃^{−8}X₄^{−12}. So the length-Z c-sum is
    longer than self-dual.
  * Reflecting gives |R|² ≈ (Z/Z′)|R′|², where R′ is a T*-type partial sum of length Z′ = N(ν₀)/Z.
  * Pure m₂-box, in Q-units with S49 on the dual: (L³) contribution ≤ 13/8 − 3x₂. With S52 on the dual: ≈ 1.5417 − 3x₂ + (continuity error).
  * m₃- and m₄-boxes gain even more: −(≈8)x₃, −(≈11)y.
* **What remains for the full sum (not yet done).**
  * (i) S52's bookkeeping at general (Q, Y) with Y slightly below Q^{3/2}. Continuity in Y/Q^{3/2} is expected and makes the near-main
    boxes ≤ 1.5417 + o(1).
  * (ii) Fix (a) uniformly (T ↔ T* transfer, S51 §20).
  * (iii) The box sum, with the multiplicities X₂X₃X₄.
  None of these looks like a new obstruction, but none is written out. The full-sum statement (L³) ≪ Q^{1.5417+ε} is therefore
  PLAUSIBLE, not PROVED. The main-box statement stands as in §6.

## 9. General (Q,Y) and the full sum (l3f/general_q.py; 01:07Z) — NUMERICAL bookkeeping
* **S52 at c-length Y with Q = Y^q, fix (a).**

  | q | Y | S ≪ Y^s | in Q-units |
  |---|---|---|---|
  | 2/3 | Q^{1.5} | 2.0250 | Q^{3.0375} |
  | 0.68 | Q^{1.471} | 2.0315 | Q^{2.9875} |
  | 0.70 | Q^{1.429} | 2.0404 | Q^{2.9149} |
  | 0.75 | Q^{1.333} | 2.0688 | Q^{2.7583} |
  | 0.80 | Q^{1.25} | 2.0904 | Q^{2.6130} |

  The bound is continuous and decreasing in Q-units as Y drops below Q^{3/2}.
* **Reflected boxes.** For an m₂-box, after the c-FE reflection (§8):
  * the (L³)-exponent is x₂ + 4x₂ + (3/2 − 4x₂)·s(q′) − 3/2, with q′ = (1−2x₂)/(3/2−4x₂);
  * its derivative at x₂ = 0 is ≈ 5 − 4·2.025 + 1.5·0.49·0.444 ≈ −2.8 < 0.
  So every box with a square part is below the main box; cube and fourth-power parts decrease faster still.
* **Fix (a) transfer.** Near q = 2/3, S ∝ Y^{≈1.7}, so √S ∝ Y^{0.85} and Σ_e N(e)·(N(e)^{−4})^{0.85} converges. Hence
  R_T = Σ_e μ(e)α⁴(e)N(e)·1_{(e,m)=1}·R*_{Z/N(e)⁴} costs only a constant.
* **Full-sum status.** (L³) ≪ Q^{1.5365+ε}, with the same STANDARD items as §2–§3 plus the routine bookkeeping in (i)–(iii) of §8.
  The only remaining loss before Q^{3/2} is still the large twists Z^{0.24} ≤ N(h₀) ≤ Z^{1/3}.

## 10. Correction: the j⁴-coprimality in S50 Step 2 is not free (l3f/jcorr.py) — exponents revised
* **What S50 assumed.** S50 wrote h = h₀j⁴ and used max_j|X(h₀,j)|² ≈ |X(h₀)|², where X(h₀,j) restricts to (m,j) = 1 and
  N(j) ≤ (H/N(h₀))^{1/4}. That is harmless for Cauchy-strength bounds (h₀-exponent 3/8). It is not harmless here, because the S52 bounds
  depend strongly on the twist (slope ≈ 1.0–1.3).
* **Exact structure (PROVED algebra).** Möbius over e | j. The piece with e | m is, by twisted multiplicativity and reciprocity
  ((e/m′)₂ = ψ_{m′}(e²)), a first moment with:
  * twist h₀e² on the d-node;
  * c-twist κ = e (χ̄_c(e));
  * m′-range Q/N(e).
  The κ-cost (≈ 0.75 per unit, l3f/kappa.py) cancels the N(e)^{−3/4} from the shorter m-range on the main term.
* **Effective first-moment exponent.** e_X^{eff}(η) = max_{0≤ε≤(⅓−η)/4} e_X(η + 2ε). Each dyadic N(e) ∼ Z^ε contributes J·|X_e|²:
  the count Z^ε times the multiplicity J/Z^ε.
* **Revised numbers (with bootstrap).**

  | variant | before §10 | after §10 |
  |---|---|---|
  | all chambers, no fix (a) | Q^{1.5417} | **Q^{1.5625}** (k₀ = 0.167) |
  | all chambers, fix (a) | Q^{1.5365} | **Q^{1.5573}** (k₀ = 0.181) |

* **What still holds.** The §7–§9 box and continuity arguments are qualitative and unaffected. The headline unconditional exponent is now
  1.5573 (fix (a)) or 1.5625 (without), still below S49's 13/8 = 1.625. The remaining loss is again the large effective twists, now
  h₀e² with e up to (H/h₀)^{1/4}.

## 11. Block-optimal regions, RMS twists, and the S(q) recursion (l3f/blocks.py, blocks_rms.py, hull8.py, recur2.py; 01:40Z)
* **Blocks.** Decompose the double sum Σ_{c,d} T*(c)χ̄_d(c)g̃(h₀;χ_d)N(c)^{−½−x}N(d)^{−½−y} into smooth dyadic blocks c ∼ C, d ∼ D.
  T* is the fix (a) family, so the c-series at fixed d is ψ̂-normalised and entire with a clean FE. Write t_C = log C, t_D = log D,
  t_K = log N(h₀) (RMS) and t_E = log N(e²) (pointwise). Each block is bounded by the better of:
  * **B1 (large sieve over c).** Bounds Σ_{c∼C}|Σ_{d∼D}·χ̄_d(h₀c)|. With the RMS over h₀ ∼ K it is
    B1 = −x t_C − y t_D + ½[max(t_K+t_C, t_D, ⅔(t_K+t_C+t_D)) − t_K].
    It uses Ξ(KC, D) for the moduli h₀c; e² enters only as a coefficient twist.
  * **B2 (mean square over d of the c-blocks).**
    B2 = −y t_D − (½+x)t_C + ½σ(t_C,t_D) + ½[max(t_K, t_D, ⅔(t_K+t_D)) − t_K − t_D].
    The last bracket is the large sieve in h₀ (vector-valued Bochner), with the true quartic Ξ(K,D) = K + D + (KD)^{2/3}.
    log S(D,C) = σ is the minimum of:
    * convexity, t_C + 2.5t_D;
    * the S(q) table at q = t_D/t_C;
    * C ≫ N(d)³: negligible, since ψ̂ is entire. Here S(Q,Y) for Y > Q^{3/2} is obtained by c-FE reflection,
      S ≈ (Y/Y′)S(Q,Y′) with Y′ = Q³/Y.
  * **D-range.** t_D ≤ (t_K+t_E+t_C)/2. This is the AFE of ψ(h₀e²c,·); for y ≥ 0 the dual half has smaller blocks.
* **Regions.** Region_d is {∀δ ∈ [0,½]: min(B1,B2)(1,δ) < 0} (directions with fixed c-length growth). Region_c is
  {∀ρ ∈ [0,1.5]: min(B1,B2)(ρ,1) < 0}. With t_K = t_E = 0 and the S49 table:
  * Region_d is x > ⅓ − y/2 for y ∈ [0,1] (old R_d′: x > ½ − 0.4y);
  * Region_c at x = 0 is y > 0.81.
  * The binding block is D = √C, where the reflected c-sums form a T*-family problem with rows = length (q = 1). S49 gives D^{7/3} there,
    which is why the recursion helps.
* **Twist exponents.** b(t_K,t_E) is the max over all blocks in the D-range of min(B1,B2). It is computed jointly, since it is not
  additive in (t_K,t_E). At the Region_d point (0.208, ¼), the RMS exponent is b = 0: the large sieve in h₀ absorbs the long d-ranges.
* **Intermediate result (hull8.py; RMS/pointwise additivity assumed, Ξ(K,D) taken optimal).** (L³) ≪ Q^{1.5234}. Both assumptions are
  removed in recur2.py; numbers there supersede this line.
* **Recursion.** S(q) for q ∈ [0.66, 1.2] is recomputed from the S50 bookkeeping, at general q (H = 2q−1 twists, j-pieces), with
  the hull built from the current table, and fed back into σ.

## 12. Recursion result (recur2.py; corrected large sieves, joint twist exponents; 01:28Z)
* **Iterations.** Each row uses the S(q) table from the previous one; S52's own bound feeds back into the block regions through σ.

  | iteration | input λ | output (L³) exponent | s(2/3) | s(0.8) | s(1) | s(1.2) |
  |---|---|---|---|---|---|---|
  | 0 (S49 table) | | 1.5253 | 2.0169 | 2.1396 | 2.271 | 2.3882 |
  | 1 | 1.5253 | 1.5156 | 2.0104 | 2.1077 | 2.25 | 2.3704 |
  | 2 | 1.5156 | **1.5156 (fixed point)** | 2.0104 | 2.1041 | 2.25 | 2.3621 |

* **Profile at the fixed point.** X_eff(η) at q = 2/3 runs from 1.128 (η = 0) to 1.196 (η = ⅓). The requirement for Q^{3/2} is
  1.29 − 0.375η, which fails only for η ≳ 0.27.
* **Why it stops.**
  * The residue-line point of Region_d is x_r ≈ s(1)/4 − 3/8. With s(1) = 2.25 that gives ≈ 0.19, and Q^{3/2} needs ≤ 1/6.
  * s(1) is itself stuck at 2.25, the pure j⁴-multiplicity value 2 + H/4 at q = 1. There the j-pieces (pointwise twists e² up to Y^{1/2})
    make the MDS useless.
* **Joint RMS over (h₀, e) does not help (checked algebraically).** For fixed e the large sieve over n = h₀c gives Ξ(KC,D), and summing
  over e cancels the family size exactly. So e² enters only through the d-range, as already modelled.
* **Status.** Σ_{N(m)≤Q}|ψ_α(m³,½+it)|² ≪ Q^{1.5156+ε}, unconditional, modulo Lemma 4 of S52_PROOF.md (non-coprime and cusp bookkeeping).
  A finer-grid verification (recur2_fine.py) is running, because the block max b is computed on a grid (§13).
* **Large-prime j-pieces, bounded trivially (checked, no gain).** |X_p|² ≤ (Q/N(p))S_p with S_p ≤ S contributes S ≤ Z²(1 + J/P). This
  only helps if P ≥ J·Z^{−(s−2)} ≈ J at s ≈ 2.015, so it changes nothing. The j⁴ structure (quartic characters cannot see fourth powers,
  the dFDH extremal mechanism) is the hard limit of the Cauchy-plus-MDS approach above Q^{3/2}.

## 13. Why the method stops near Q^{1.51} (01:33Z) — PROVED as exponent algebra
* **The residue-line point.** In Region_d it is set by the block D = √C. There the reflected c-sums form the rows = length problem
  S(Y,Y), giving x_r ≈ s(1)/4 − 3/8.
* **The circularity at q = 1.** The S50 main term is Q^{3/4}Z^{½+x_r} = Y^{1.25+x_r}. It fits the q = 1 budget (X ≤ 1.375 − 0.375η) only if
  x_r ≤ 1/8, i.e. only if s(1) ≤ 2 already. So the recursion cannot lower s(1) below its generic value 2 + H/4 = 2.25, and x_r stays ≈ 0.19.
* **Consequence at q = 2/3.** X(η = ⅓) ≈ 1.19 > 7/6, so the large twists keep s(2/3) > 2. Q^{3/2} by this route therefore needs
  S(Y,Y) ≪ Y^{2+ε}: sub-problem 1 of L3_REPORT §12, the optimal T-large sieve at M = 1. That is proved on GRH (§10b) and open
  unconditionally.
* **Summary.** S52 turns any unconditional progress on sub-problem 1 directly into progress on (L³). With sub-problem 1, x_r ≈ 1/8 and
  the route plausibly reaches Q^{3/2}; the intermediate-η j-pieces would still need checking (HEURISTIC). Going below Q^{3/2} still needs
  a non-Cauchy idea.
* **Fine-grid verification (recur2_fine.py; 01:33Z).**
  * Grids: block c-grid 141 points, d-grid 81 per column, 17 η-points, 5 ε-points at q ∈ {2/3, 0.7, 0.8, 0.9, 1}.
  * The fixed point is unchanged: s(2/3) = 2.0104, s(1) = 2.2496, (L³) ≪ Q^{1.5156}. Grid coarseness was therefore not
    under-estimating the block maxima b.
* **e-tail cut (recur4.py; 01:34Z).** Bound the j-pieces with e > Y^{ε₀} trivially, |X_e|² ≤ (Q/N(e))S_e. The tail is then affordable
  once 2ε₀ ≥ η + (H−η)/4 − (s−2), and s is solved self-consistently by bisection.
  * First pass: s(1) drops from 2.2496 to 2.2188, s(0.8) to 2.1004, s(1.2) to 2.3063. s(2/3) is unchanged at 2.0104 on that pass,
    because the lower s(1) enters the regions only on the next iteration.
  * The continued run is recur5.py.

## 14. Correction to §12: the generic twist range at q > 2/3 (01:36Z)
* **The error.** In recur2/recur4, S50's generic large-sieve part (twists beyond the MDS range) was costed as 2 + (H−k₀)/4. That is
  correct only at q = 2/3. In general it is J·Ξ₄(K,Q) with Ξ₄ ≥ (KQ)^{2/3}, giving s ≥ 2 − q + (H−κ)/4 + log_Y Ξ₄(Y^κ, Y^q).
  * At q = 1 this forces s(1) ≥ 7/3, the S49 value. So the values s(1) = 2.25 and 2.2188 were too optimistic.
  * The tables at q ≳ 0.7 used in the block regions were also too optimistic.
  * The q = 2/3 bookkeeping itself is unaffected, since (KQ)^{2/3} ≤ Q there.
* **Consequence.** The §12 fixed point 1.5156 is superseded. The true fixed point lies between 1.5156 and 1.5253; 1.5253 is iteration 0,
  which used the S49 table for q ≥ 1 and is unaffected.
* **Fix (recur6.py).**
  * The generic cost is computed per twist range with Ξ₄.
  * The e-tail is bounded by the large sieve over h₀, Ξ₄(K, Q/E)·S_e, instead of K copies of the trivial bound.
* **e-tail via the large sieve over h₀ (PROVED algebra).** For fixed e,
  Σ_{h₀∼K}|X_e(h₀)|² ≤ Ξ₄(K, Q/E)·S_e, where e² is absorbed into the m′-coefficients. Hence the pieces with e ∼ E cost
  s ≥ 2 − q + (H−η)/4 − ε + log Ξ₄(η, q−ε). At q = 2/3 that is 2 + (H−η)/4 − 2ε, so only e ≤ Y^{((H−η)/4 − (s−2))/2} need the
  twist-costly MDS bound.
* **A bug in the first recur6 pass.** It printed s(2/3) = 2.0011 (Q^{1.5016}), but costed the generic range (k₀,H] only at grid points
  strictly beyond k₀, missing the sup at κ → k₀⁺. The boundary is fixed and the η-grid refined; that run is in progress. A hand estimate
  puts the corrected value near s(2/3) ≈ 2.012, i.e. Q^{≈1.518}.

## 15. Corrected recursion (recur6.py: generic range via Ξ₄, e-tail via the large sieve over h₀, boundary fixed; 01:41Z)
* **Iteration 0 (S49 table for q ≠ 2/3).** s(2/3) = 2.007, i.e. (L³) ≪ Q^{1.5104}. Also s(0.8) = 2.1334, s(1) = 7/3, s(1.2) = 2.5333.
* **Direct check (check_x.py).** At q = 2/3:

  | η | 0 | 0.222 | 0.278 | 0.306 | 0.333 |
  |---|---|---|---|---|---|
  | X(η,0) | 0.917 | 1.104 | 1.156 | 1.182 | 1.208 |

  The MDS budget at k₀ = 0.3056 is 1.184, so the result is consistent. With a continuous k₀ the optimum is s − 2 ≈ 0.0067, i.e. Q^{1.510}.
* **Why Q^{3/2} is out of reach here (exponent algebra).**
  * It needs X(⅓,0) ≤ 7/6, i.e. x_r ≤ 1/6, i.e. s(1) ≤ 13/6.
  * At q = 1, twists K ∼ Y force the balanced quartic large sieve. That gives s(1) ≥ 1 + log Ξ₄(Y,Y) = 7/3 today.
  * Even an ideal general-coefficient large sieve gives only s(1) ≥ 9/4, by the dFDH lower bound Ξ₄ ≫ Y^{5/4}.
  * So Q^{3/2} via this route needs T-specific input at rows = length with twists up to Y.

## 16. Third correction: the theta polar term in the d-AFE (recur7.py)
* **The gap.** In the 4 chambers whose d-node is untwisted (ℓ_d = 0, including the identity), the AFE of ψ(h₀c,w) contains the theta
  polar term ρ(h₀c)X^{¾−w}, with X = N(h₀c)^θ. Summed over c, using Cauchy–Schwarz and the Rankin–Selberg mean square of ρ, it
  converges only for x > 3/8 + θ(¼ − y). On the theta line itself that means x > 3/8 for every θ.
* **What it invalidates.** The block Region_d (§11) ignored this term, so its points with x < 3/8 near y = ¼ in polar chambers were
  unjustified. Every number from §11 onward (1.5234, 1.5253, 1.5156, 1.5104, and recur6's 1.5016) is invalid until recomputed.
* **What it leaves alone.** The older regions satisfy the constraint (R_d: x > ½; R_d′: x > ½ − 0.4y ≥ ½ − y/2 for y ≥ 0), so the
  §10 values (1.5625 without fix (a), 1.5573 with it) remain valid.
* **Fix.** Region_d is now a union over the AFE balance θ ∈ [0.1, 1]:
  * main blocks for δ ≤ θ, and FE-dual blocks for δ′ ≤ 1−θ, evaluated at (x+y, −y);
  * in polar chambers, the polar constraint, with its own RMS twist exponents:
    t_C(3/8 − x + θ(¼−y)) + t_K(θ(¼−y) − 1/8) + t_E(θ(¼−y) + 3/8).
  * With the S49 table: polar chambers have x_min = 0.475, 0.375, 0.205, −0.131 at y = 0, ¼, ½, 1; non-polar chambers keep
    0.333, 0.208, 0.083, −0.167.

## 17. Audit fixes, the corrected recursion, and the limit Q^{3/2} (recur9–11, cone_check.py, cocycle.py; 02:19Z)
* **Three more pipeline errors, found by re-reading recur6/recur7 against the audit list.**
  1. *Unproved start.* recur6 and recur7 began the bootstrap at λ = 1.5253, a value that was never proved. A bootstrap may only start from
     a proved bound. recur9 starts from S49's 13/8.
  2. *Table extrapolation.* S(q) was read with np.interp, which clamps outside [0.40, 1.6]. For q > 1.6 (directly, or through the
     reflection from q ∈ (1/3, 0.421)) it returned S(1.6) = 2.933 instead of s49(q); at q = 3 the right value is 4.33. Fixed: s49 is used
     outside the table.
  3. *Polar term on reflected R_d points.* In polar chambers, an FE-reflected R_d point whose image has y′ < ¼ now carries the polar
     exponent of that image.
  * None of the three binds. With them fixed, the sequence is the same as recur6's: 1.5104, 1.5053, 1.5027, 1.5016.
* **§10 re-run from the proved start (jcorr_valid.py).** With fix (a): 13/8 → 1.5625 → 1.5573. Without: 1.5677 → 1.5625. Both §10
  values stand.
* **Chamber bookkeeping check (cocycle.py, PROVED by exhaustive check).**
  * The h₀-factor rule (s_d at chamber w gives N(h₀)^{−y′}; s_c gives nothing) is path-independent: 12 chambers, every reduced word of
    length ≤ 12 checked, 0 inconsistencies, and (s_c s_d)⁶ returns φ = 0.
  * First-principles reason: h₀ enters every chamber only through χ̄_d(h₀) in the d-node. The Kubota d-FE keeps the numerator, and the
    c-FE does not touch h₀.
* **recur9 (proved start; table read correctly; samples at the Region_d boundary + 10⁻³).**

  | iteration | 0 | 1 | 2 | 3 | 4 |
  |---|---|---|---|---|---|
  | (L³) exponent | 1.5104 | 1.5053 | 1.5027 | 1.5016 | 1.5016 |
  | s(2/3) | 2.0070 | 2.0035 | 2.0018 | 2.0010 | 2.0010 |

  The stop at 2.0010 comes from the 10⁻³ sampling offset (next bullet). recur10 continues from recur9's table with offset 10⁻⁴ and
  bisection 2⁻²⁰. It gives 1.5008, 1.5004, 1.5002, 1.5002, with s(2/3) = 2.0001.
* **The mechanism (exponent algebra, PROVED given the regions; diag9.py, diag9b.py, limit9.py).**
  * *Binding twist.* At q = 2/3 the only binding constraint is the MDS at the top twist η = H = ⅓. All smaller η have slack ≥ 0.02 at
    the final state, and there the e-tail is empty.
  * *Binding hull point.* It mixes one Region_d point (x′, 0), with x′ just above Region_d's boundary at y′ = 0, from two non-polar
    chambers:
    * chamber s_c (twist (−1,3), φ = 0), pulled back to (−x′, 3x′) with b = b_loc;
    * chamber s_d s_c s_d s_c (twist (−2,3), φ = (−3,−3)), pulled back to (x′, −3x′) with b = b_loc + 6x′η.
    * Hitting y = ¼ forces weights 5/8 and 3/8. That gives x_r = −1/12 (exactly on the (3,1) normaliser line 3x+y = 0) and
      b = b_loc + (3x′ − ¼)η.
    * Hence X(η) = 11/12 + (3x′ − ¼)η + b_loc(η), and X(⅓) = 5/6 + x′ + b_loc.
  * *Boundary of Region_d at y′ = 0.* It is x′ = s(1)/4 − ¼ = ⅓ with S49's s(1) = 7/3. The binding direction is D = C^{1/2}, through
    the reflection S(½) = ½ + ½S(1).
  * *b_loc at x′ = ⅓ + δ′.* It is set by the block (t_C, t_D) = (3η, 2η), i.e. D = C^{2/3}, which is the (L³) shape itself. There
    B1 = ½ − x′ is not binding and B2 = 3η(½(s(2/3) − 2) − δ′).
  * *Contraction.* With k₀ = H (as in the code's grid), s_{n+1} − 2 = max(δ′, (s_n − 2)/2): the observed halving, with fixed point
    2 + δ′. With k₀ chosen continuously, the MDS covers η ≤ ⅓ − 4δ (MDS(η) = 1.625 + 1.125η + O(δ)), the generic large sieve covers the
    rest, and the factor improves to about 1/11.
  * *Why §13 and §15 were wrong to call this impossible.* They only tested the residue-line point of Region_d itself (x_r ≈ s(1)/4 − 3/8).
    The cross-chamber combination reaches y = ¼ at x_r = −1/12 and pays b = ¾η instead. At η = ⅓ that is exactly the needed
    x_r + b ≤ 1/6.
* **The limiting table and the cone inequality (cone_check.py: fresh code, independent of recur*).**
  * The limit is S(q) = q + 4/3 for q ≥ 2/3, which equals s49 for q ≥ 5/6 and is the generic large-sieve floor gen(H) for
    q ∈ [2/3, 1]. By reflection it is 2q + 2/3 on [⅓, ⅔] and 1 + q below ⅓.
  * The recursion drives the table onto this line. In recur9's final state, s(0.7, 0.72, 0.75, 0.8) = 2.0334, 2.0533, 2.0834, 2.1334.
  * For q ∈ [2/3, 5/6], the generic large sieve gives q + 4/3 for η ≥ 5/3 − 2q. The two-chamber combination gives
    X = 0.75q + 5/12 + ¾η on the rest, which fits the MDS budget exactly when q ≥ 2/3.
  * Along the boundary ray of blocks at x′ = ⅓, B2 = (1+t)/4·(S(q′) − q′ − 4/3) + (terms that vanish on the ray): every block is exactly
    marginal on the limiting line.
  * At x′ = ⅓ + δ′ with the limiting table, the worst asymptotic direction is −δ′, so the point is strictly interior. The maximum of
    min(B1,B2) over the whole cone {t_D ≤ ½(t_K + t_C)} is 0, so b_loc ≡ 0 for every η. Checked for δ′ = 10⁻², 10⁻³, 10⁻⁴, on a slice
    of 1201 × 201 points plus a tail to t_C = 40.
  * With the table raised by Δ on [2/3, 5/6), b_loc ≈ 2.6Δ·η, from blocks reflecting to q′ ≈ 0.81. The q = 2/3 bookkeeping then returns
    s − 2 ≈ 0.16Δ, a strict contraction.
* **Interpolation (rigor).** Only step-up interpolation is justified (S(Q,Y) is increasing in Q); linear interpolation of the table is not.
  Because the binding blocks sit exactly on the line, this matters. recur11 reruns from the proved start with a step-0.01 grid on
  [2/3, 0.9] and step-up reads; see the result line below. A step h costs Δ ≈ h in the table, hence ≈ 0.16h in s(2/3), so any ε is
  reached with a fine enough grid.
* **∞-types.** α(c) = c/|c|. The chamber families are α^ℓ-twisted, with ℓ_c ∈ {±1, ±2}, never 0. Everything used is uniform in |ℓ| up
  to polynomial factors: the large sieve, Poisson against the smooth angular weight, DDHL (4.28), and the normalisers L(·, α^{4ℓ|f|}),
  which are nontrivial. So one table serves all ℓ, and each stage only calls finitely many ∞-types.
* **S51's (5,2) obstruction does not enter.** It is the convergence of E′ = naive·ΠL/CG, a factorisation through Chinta–Gunnells. S52
  continues the naive series itself. Non-coprime local terms enter only inside the region representations (Lemma 4(b)), and at the
  binding local point (y′ = 0) they give Π_{π|ν}(1 + O(1)) ≪ N(ν)^ε.
* **Status.**
  * The S52 framework gives Σ_{N(m)≤Q}|ψ_α(m³,½+it)|² ≪ Q^{3/2+ε} for every ε > 0, which is the Cauchy floor.
  * PROVED as exponent algebra: the mechanism and the contraction above.
  * NUMERICAL: the LP hulls at the non-binding twists, and the block maxima on grids.
  * Still open, as before: Lemma 4 constants, boxes and the T↔T* transfer, and a rigorous-numerics pass (exact PL arithmetic for the cone
    and hull checks, plus a continuous-q bookkeeping certificate).
  * (L³) itself (Q^{1+ε}) remains OPEN. Below Q^{3/2} the Cauchy step in S50 has to be replaced.
* **recur11 / recur12 (step-up reads, proved start).**
  * recur11, grid step 0.01: 1.5104 → 1.5055 (fixed).
  * recur12, step 0.004 plus extra twists near H, continuing from recur11: 1.5055 → 1.5020 (iteration 0; the run was then cut by a
    worker restart).

## 18. The limit as a theorem, without LP or grids (S52_LIMIT_PROOF.md, cone_exact.py; 02:34Z)
* **The proof is analytic.** It uses one explicit hull point: the Region_d boundary point (⅓ + δ′, 0) pulled back from chambers s_c and
  "dcdc". Its inputs are:
  * an exact cone lemma, b_loc(t_K,t_E) ≤ 3Δ(t_K + t_E) + ½t_E. This is a rational vertex enumeration at δ′ = Δ = 0 (7323 vertices,
    max 0), plus monotonicity in δ′, plus ∂B2/∂Δ ≤ 3(t_K + t_E) from the sector geometry;
  * closed-form S50 bookkeeping at every q ∈ [⅔, 5/6].
* **The recursion is one number.** For the table offset Δ (S ≤ q + 4/3 + Δ on [⅔, 5/6)): Δ ↦ (δ′+Δ)/(5.5 + 12(δ′+Δ)), starting from
  S49's Δ₀ = 1/12.
* **Result.** For every ε > 0, (L³) ≪ Q^{3/2+ε}. This is PROVED given the framework (S52_PROOF Lemmas 1–3 and 5, the S50 bookkeeping,
  S49) and the open items Lemma 4 and the box / T↔T* transfer.
* **What is no longer needed for the limit.** The LP hulls, the polar chambers, the R_c bootstrap and table interpolation. The numerical
  recursions remain as cross-checks.
