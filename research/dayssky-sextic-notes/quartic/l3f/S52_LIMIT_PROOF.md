# S52 limit theorem: Σ_{N(m)≤Q}|ψ_α(m³,½+it)|² ≪ Q^{3/2+ε} (DaysSky with Claude assistance, 16 Oct 2026)

> **CORRECTION (16 Oct, 03:38Z): the results here are CONDITIONAL on task #26.** That is the Chinta–Gunnells G₂ WMDS over Q(i) with the BD
> sign and the α-twist. The hull must be run on the CG series, because the naive series cannot reach the chamber-dcdc points. See
> S52_CORRECTION_NAIVE_VS_CG.md.


Labels:
* **PROVED** means a complete argument, given the S52 framework.
* **EXACT** means verified by exact rational computation, with the script named.
* **FRAMEWORK** means the S52 inputs listed in §0, used as stated there.

This replaces the numerical recursion (recur*.py) for the limit statement. It uses no LP hulls, no polar chambers, no R_c bootstrap,
and no interpolation of tables.

## 0. Framework inputs (S52_PROOF.md, S52_AUDIT.md)
* **(F1) Bochner with growth.** S52_PROOF Lemma 1 and Corollary 1a (vector-valued, RMS over h₀).
* **(F2) Block bounds B1 and B2.** S52_PROOF Lemma 2, with σ = log S*(D,C). The inputs for σ are:
  * reflection: σ = t_C − t_C′ + t_C′·S(t_D/t_C′) with t_C′ = 3t_D − t_C, for t_D/t_C ∈ (⅓, ⅔);
  * negligibility for t_C > 3t_D;
  * any proved bound S(q) for the T* family, uniformly over the ∞-types α^ℓ with ℓ ≠ 0 (S52 §17).
* **(F3) Region_d.** A chamber-coordinate point (x′, y′) lies in Region_d(θ) if every asymptotic direction has min(B1,B2) < 0. Its RMS
  twist exponent b_loc(t_K, t_E) is the maximum of min(B1,B2) over the blocks with t_D ≤ θ(t_K + t_E + t_C) (S52_PROOF Lemma 3).
* **(F4) Chambers and h₀-factors.** The 12 W(G₂) chambers, with Z(v) = N(h₀e²)^{φ_w·v}·Z_w(M_w v). φ_w is path-independent
  (cocycle.py, EXACT).
* **(F5) Lemma 4 bookkeeping** (non-coprime terms and cusps), and the box and T↔T* transfer of S52 §8–§9. These are the open items of
  S52_AUDIT.
* **(F6) Evaluation (Lemma 5).** For the first moment X with twist h₀e² (RMS over h₀ ∼ Y^η, e fixed with N(e²) = Y^{2ε}):
  X ≤ max(main, dual), where:
  * main = ¾q + ½ + x_r + b at a hull point (x_r, ¼);
  * dual = ½ + ½q + x + qy + b at a hull point with 0 ≤ y < ¼;
  * every root functional 3x+y, 2x+y, x+y, y is ≥ 0 at both points (fix (a)).
* **(F7) S50 bookkeeping at Q = Y^q, H = 2q−1.** S ≪ Y^s holds if there is k₀ ∈ [0, H] such that:
  * (MDS) for η ≤ k₀ and ε ≤ e₀(η): ½(1 + (1−q) + η + (H−η)/4 + 2X(η,ε)) ≤ s;
  * (tail) for ε ≥ e₀(η): 2 − q + (H−η)/4 − ε + log Ξ₄(η, q−ε) ≤ s;
  * (generic) for η ≥ k₀: 2 − q + (H−η)/4 + log Ξ₄(η, q) ≤ s;
  * here log Ξ₄(a, b) = max(a, b, ⅔(a+b)).
* **(F8) S49.** S(q) ≤ s49(q) = 1 + max(1, q + ⅓, q/2 + ¾). In particular s49(q) = q + 4/3 for q ≥ 5/6, and S(1) ≤ 7/3.

Write L(q) = q + 4/3. A table offset Δ means a proved bound S(q) ≤ L(q) + Δ for every q ∈ [⅔, 5/6) (a continuum, with no grid),
together with S ≤ s49 = L on [5/6, ∞).

## 1. Lemma A (region membership) — PROVED
For x′ > ⅓, the point (x′, 0) lies in Region_d(½) of every non-polar chamber, for every table offset Δ ≥ 0.

*Proof.*
1. At y′ = 0 the FE-dual blocks coincide with the main blocks, since (x′+y′, −y′) = (x′, 0) and θ = 1−θ = ½.
2. Take an asymptotic direction (t_C, t_D) = (1, δ), δ ∈ [0, ½], with t_K = 0.
3. For δ < ⅓ the block is negligible.
4. For δ ∈ [⅓, ½], σ comes from the reflection with q′ = δ/(3δ−1) ≥ 1, where S(q′) ≤ q′ + 4/3 by S49. This gives σ ≤ ⅔ + 2δ, with
   no Δ term.
5. Hence B2 ≤ −(½ + x′) + ⅓ + δ ≤ ⅓ − x′ < 0. ∎

## 2. Lemma B (cone lemma) — PROVED (EXACT plus two observations; cone_exact.py)
For every δ′ ≥ 0, every table offset Δ ≥ 0, and all t_K, t_E ≥ 0, at (x′, y′) = (⅓ + δ′, 0) with θ = ½:
  b_loc(t_K, t_E) ≤ 3Δ(t_K + t_E) + ½t_E.

*Proof.*
1. **Homogeneity.** Everything is homogeneous of degree 1 in (t_K, t_E, t_C, t_D), so it suffices to take t_K + t_E = 1, i.e.
   t_K = 1−ρ, t_E = ρ. The blocks are then 0 ≤ t_D ≤ ½(1 + t_C).
2. **Conservative σ.** Take σ = t_C·S_Δ(t_D/t_C), with S_Δ = L + Δ·1[q′ ∈ [⅔, 5/6)] read through the reflection for q ∈ (⅓, ⅔).
   Dropping the convexity cap only enlarges σ. This σ is linear on the four sectors q ∈ [⅓, 5/9], [5/9, ⅔], [⅔, 5/6], [5/6, ∞]; it is
   2t_D + ⅔t_C plus, respectively, 0, Δ(3t_D − t_C), Δt_C + (⅔t_C − t_D), and (⅔t_C − t_D). The sectors are taken closed, which is
   conservative at the jumps.
3. **(i) Exact case δ′ = 0, Δ = 0.** max(min(B1,B2) − ρ/2) = 0.
   * The function is piecewise linear, so its maximum is attained at a vertex of one of the 72 cells (4 sectors × 3 × 3 max-pieces ×
     2 sides of B1 = B2).
   * All 7323 feasible vertices are enumerated in rational arithmetic, with the cap t_C ≤ 10⁵ (the same result with 10³).
   * Every recession direction has slope ≤ ⅓ − x′ ≤ 0 (Lemma A), so the cap loses nothing.
4. **(ii) Monotone in δ′.** B1 and B2 depend on x′ only through −x′t_C and −(½ + x′)t_C.
5. **(iii) Δ-sensitivity.** Only B2 depends on Δ, through ½σ.
   * On the sector [5/9, ⅔], the block condition t_D ≤ ½(1 + t_C) together with t_D ≥ 5t_C/9 forces t_C ≤ 9. So
     ∂σ/∂Δ = 3t_D − t_C ≤ 3/2 + t_C/2 ≤ 6.
   * On [⅔, 5/6] it forces t_C ≤ 3, so ∂σ/∂Δ = t_C ≤ 3.
   * Hence ∂B2/∂Δ ≤ 3 = 3(t_K + t_E), and min(B1,B2) rises by at most 3Δ(t_K + t_E).
6. **Spot checks (EXACT).** With δ′ ∈ {1/100, 1/1000} and Δ ∈ {0, 1/1000, 1/50, 1/12}, the maximum of min(B1,B2) − 3Δ − ρ/2 is
   max(−3Δ, −9δ′) < 0. It is attained at the origin, or at (t_C, t_D) = (9, 5), the block that reflects to q′ = 5/6. ∎

## 3. Lemma C (the two-chamber hull point) — PROVED
* **The two points.** Let x′ = ⅓ + δ′. Then:
  * chamber s_c (twist (−1,3), φ = 0) contributes v₁ = (−x′, 3x′), since s_c v₁ = (x′, 0);
  * chamber s_c s_d s_c s_d (word "dcdc", twist (−2,3), φ = (−3,−3)) contributes v₂ = (x′, −3x′), since M v₂ = (x′, 0) and φ·v₂ = 6x′.
* **Exponents.** By (F1)–(F4), the h₀-exponents are:
  * b(v₁) = b_loc;
  * b(v₂) = b_loc + 6x′(t_K + t_E).
* **The convex combination.** Take w·v₁ + (1−w)·v₂ with 2w − 1 = 1/(12x′). It has y = ¼, x_r = −1/12, and
  b = (3x′ − ¼)(t_K + t_E) + b_loc.
* **Evaluation constraints.** On the whole segment 3x + y = 0, 2x + y ≥ 0, x + y ≥ 0 and y ≥ 0 hold for y ∈ [0, ¼]. So (F6) applies.
* **Dual term.** It runs along the same segment (w ∈ [½, 5/8)), and dual ≤ main, because the dual expression is affine in w and equals
  the main one at w = 5/8.
* **The bound.** With Lemma B, writing a = ¾ + 3δ′ + 3Δ:
  X(η, ε) ≤ ¾q + 5/12 + a(η + 2ε) + ε.
* **Side costs.** The c-twist κ = e costs at most (5/8)x′ε < ¾ε. The shortening of the m-range to Q/N(e) saves ¾ε on the main term and
  ≥ ½ε on the dual term. Ignoring both is therefore conservative.

## 4. Lemma D (one step of the bookkeeping) — PROVED (algebra)
If the table offset is Δ, then for every q ∈ [⅔, 5/6], S(q) ≤ L(q) + Δ′ with
  Δ′ = κ/(5.5 + 12κ),  κ = δ′ + Δ.

*Proof.*
1. **Setup.** Fix q ∈ [⅔, 5/6], write u = q − ⅔ and s = L(q) + σ.
2. **Generic range.**
   * For η ≤ q/2, log Ξ₄(η, q) = q, so generic(η) = 2 + (H − η)/4. This is ≤ s iff η ≥ k₀ := 5/3 − 2q − 4σ.
   * For η ∈ [q/2, H], generic(η) is affine increasing, with generic(H) = L(q). So the generic range [k₀, H] costs ≤ s.
3. **Tail.** For η ≤ k₀ and ε ≤ (H−η)/4, log Ξ₄(η, q−ε) = q − ε, so tail = 2 + (H−η)/4 − 2ε. It is ≤ s for
   ε ≥ e₀(η) = (k₀ − η)/8.
   * The identity log Ξ₄(η, q−ε) = q − ε needs η ≤ (q−ε)/2. For ε ≤ (H−η)/4 that follows from η ≤ (2q+1)/7.
   * And k₀ ≤ (2q+1)/7 is equivalent to 16q + 28σ ≥ 32/3, which holds for every q ≥ ⅔.
4. **MDS.** By Lemma C, MDS(η, ε) = 31/24 + q/2 + (3/8 + a)η + (2a+1)ε.
   * Along ε = e₀(η) its slope in η is (3/8 + a) − (2a+1)/8 = ¼ + ¾a > 0, so the worst case is η = k₀, ε = 0.
   * The condition MDS(k₀, 0) ≤ s reads (3/8 + a)(5/3 − 2q − 4σ) ≤ q/2 + 1/24 + σ, i.e.
     σ(5/2 + 4a) ≥ (a − ¾)/3 − u(5/4 + 2a).
5. **Conclusion.** The right side is largest at u = 0, where σ = (a − ¾)/(3(5/2 + 4a)) = κ/(5.5 + 12κ). For larger u it is smaller, and
   it is 0 once u ≥ κ/(2.75 + 6κ). ∎

*Cross-check (limit_map.py, grids of 401 η × 41 ε, before the restart):* Δ = 1/12 → 0.012917 (formula: 0.012834) → 0.002346
(formula: 0.002302). The grid values are slightly larger, as expected for a discretised search.

## 5. Theorem — PROVED, given the framework (F1)–(F8)
For every ε > 0, Σ_{N(m)≤Q}|ψ_α(m³,½+it)|² ≪_ε Q^{3/2+ε}(1+|t|)^{O(1)}.

*Proof.*
1. **Start.** S49 gives the table offset Δ₀ = max over [⅔, 5/6] of (s49 − L) = 1/12.
2. **Iteration.** Fix δ′ > 0. Apply Lemma D repeatedly: Δ_{n+1} = (δ′ + Δ_n)/(5.5 + 12(δ′ + Δ_n)). Each step uses only the previous
   step's proved table.
3. **Convergence.** The map κ ↦ κ/(5.5 + 12κ) has slope ≤ 1/5.5 = 2/11. So Δ_n ≤ Δ* + (2/11)ⁿ/12, where the fixed point satisfies
   Δ* ≤ δ′/4.5. After finitely many steps, S(⅔) ≤ 2 + δ′/2.
4. **Translate.** By S50, (L³) ≪ Q^{λ+ε} with λ = 1.5(S(⅔) − 1) ≤ 3/2 + δ′.
5. **Constants.** The implied constants depend on δ′ (through Region_d's margin δ′ and the number of steps) and on the finitely many
   ∞-types used. ∎

## 6. What this does and does not give
* **Q^{3/2} is the floor of S50's Cauchy step.** At η = H, the generic large sieve forces s(⅔) ≥ 2. It equals the GRH value
  (L3_REPORT §10b).
* **(L³) itself, Q^{1+ε}, is not reached.** Going below 3/2 needs a replacement for the Cauchy step in S50.
* **θ is unchanged.**
