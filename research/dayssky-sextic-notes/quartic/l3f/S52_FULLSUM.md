# S52, full sum: Σ_{N(m)≤Q}|ψ_α(m³,½+it)|² ≪ Q^{3/2+ε} over all rows (DaysSky with Claude assistance, 16 Oct 2026)

> **CORRECTION (16 Oct, 03:38Z): the results here are CONDITIONAL on task #26.** That is the Chinta–Gunnells G₂ WMDS over Q(i) with the BD
> sign and the α-twist. The hull must be run on the CG series, because the naive series cannot reach the chamber-dcdc points. See
> S52_CORRECTION_NAIVE_VS_CG.md.


S52_LIMIT_PROOF.md treats squarefree rows, the main box. This note adds every other row m = m₁m₂²m₃³m₄⁴. It uses the same framework
(F1)–(F8) and the same open items (Lemma 4 constants and the T↔T* transfer). Labels as there.

## 1. Reduction of a box to twisted main-box problems (PROVED; S49 "Corollary for rows")
* **The row's series.** Let m = m₁m₂²m₃³m₄⁴, with m₁, m₂, m₃ squarefree and pairwise coprime, and m_i ∼ X_i. Then m³ = ν₀x⁴ with
  ν₀ = m₁³m₂²m₃.
* **AFE lengths.** ψ_α(m³,½) is a sum over dyadic lengths Y ≤ L₀ := N(ν₀)^{1/2} of R_Y(m)/√Y, plus the same with dual coefficients. The
  dual coefficients are a T-type family with the same index, by the self-dual c-FE.
* **Removing the m₄ coprimality.** Use Möbius over d | rad x, then T(dc′) = ±χ_{c′}(d²)T(d)T(c′). This gives
  |R_Y(m)| ≤ Σ_d |R^{ϑ}_{Y/N(d)}(m₁)|, where:
  * R^{ϑ}_{Y′}(m₁) := Σ_{c′} T(c′)ϑ(c′)χ̄_{c′}(m₁)W(N(c′)/Y′);
  * ϑ(c′) = χ̄_{c′}(n₁) is a character of c′, with n₁ = m₂²m₃³d² up to fourth powers.
* **The twisted family.** For fixed (m₂, m₃, m₄, d), the box slice is the main-box problem for the twisted family T·ϑ, with squarefree
  rows m₁ ∼ X₁. Its index-conductor is F := N(m₂)²N(m₃)N(d)², i.e. the c-series at row m₁ has conductor N(m₁)³F.

## 2. The twisted limit lemma
Fix ϑ with index-conductor F, and write t_F = log F (log_Y units where noted). Let Δ be the current offset of the untwisted squarefree
table (S52_LIMIT_PROOF §0), and κ the current twist coefficient of the twisted table (defined in (D^ϑ)).

* **(A^ϑ) Region membership — PROVED.** Asymptotic directions do not see F (it is fixed while t_C → ∞). So Lemma A holds verbatim.
* **(C^ϑ) Hull point — PROVED, and EXACT for the cocycle (fcocycle.py).**
  * The twist ϑ lives on the c-node and persists in every chamber with the same F. The d-FE does not touch it, and the c-FE maps a
    Kubota series of index d³n₁ to one of the same index.
  * So s_c at chamber w contributes F^{−(M_w v)_x}. This F-cocycle is path-independent: 12 chambers, every reduced word of length ≤ 12,
    0 inconsistencies.
  * At the two-chamber point of S52_LIMIT_PROOF Lemma C, the F-exponents are ⅓ (chamber s_c) and 1 (chamber dcdc). With weights 5/8 and
    3/8, the cost is 7/12 per unit t_F.
* **(B^ϑ) Twisted cone lemma — PROVED (reduction to Lemma B).**
  * **Claim.** b_loc^ϑ(t_K,t_E;t_F) ≤ 3Δ(t_K+t_E) + ½t_E + (5/6 + κ/2 + 2Δ)·t_F.
  * **Proof.** Put t_D′ = t_D + t_F/3, noting that D³F = (D·F^{1/3})³. The twisted σ^ϑ is one of four cases:
    * direct (q ≥ ⅔): t_C·L_Δ(q) + κt_F;
    * twisted reflection (t_C′ = 3t_D + t_F − t_C ∈ [0, t_C), with t_D/t_C′ ≥ ⅔): t_C − t_C′ + t_C′L_Δ(t_D/t_C′) + κt_F;
    * the Cauchy floor 2t_C + κt_F in the remaining gap;
    * t_D (single term) at t_C = 0, and negligible for t_C > 3t_D + t_F.
  * In every case σ^ϑ(t_C,t_D) ≤ σ_{L_Δ}(t_C, t_D′) + (⅔ + κ)t_F:
    * direct, because L is increasing;
    * reflection, because t_C′ = 3t_D′ − t_C is the untwisted reflection length at t_D′, and t_C′ < t_C forces t_D′/t_C < ⅔;
    * gap: either t_D′/t_C ≥ ⅔, so σ_L ≥ 2t_C, or t_C < 1.5t_D + t_F, so 2t_C − (2t_D′ + ⅔t_C) < ⅔t_F;
    * at t_C = 0, trivially.
  * The bracket g(t_D) = max(t_K, t_D, ⅔(t_K+t_D)) − t_D is non-increasing with slope ≥ −1. Hence
    B2^ϑ(t_D) ≤ B2(t_D′) + ½(1+κ)t_F, and B1(t_D) ≤ B1(t_D′).
  * The block (t_C, t_D′) lies in the untwisted cone with t_E′ = t_E + ⅔t_F. Lemma B there gives
    3Δ(t_K + t_E′) + ½t_E′ = 3Δ(t_K+t_E) + ½t_E + (2Δ + ⅓)t_F. Adding the pieces gives the claim. ∎
  * **Check (twist_cone3.py, grid).** With the twisted table taken literally, the grid maximum is c₃ = 0.666 ≤ 5/6 + κ/2.
* **(D^ϑ) Twisted bookkeeping — PROVED (algebra as in Lemma D).**
  * By (B^ϑ) and (C^ϑ), X(η,ε) gains c·t, with c := 7/12 + 5/6 + κ/2 + 2Δ and t = t_F/log Y. S50's Cauchy step drops ϑ by positivity
    (|ϑ| ≤ 1), so the rest of the bookkeeping is the untwisted one at (rows, length).
  * Lemma D's algebra then gives:
    * for q = ⅔ + u ≥ ⅔: S^ϑ ≤ Y^{L(q) + Δ′ + e}, with e = max(0, (c·t − 2.75u)/5.5);
    * for q = ⅔ − v < ⅔: S^ϑ ≤ Y^{2 + Δ′ + e′}, with e′ = max(0, (c·t − 2.75v)/5.5) (the Cauchy floor, with the same k₀ and tail
      algebra).
  * In absolute terms the twist costs at most (c/5.5)·t_F. That is the twisted table's κ for the next stage.
* **Closing the induction.** κ_{n+1} = c_n/5.5 and c_n = 17/12 + κ_n/2 + 2Δ_n. The fixed point is c* = (11/10)(17/12 + 2Δ*) ≤ 1.56 + 2.2Δ*.
  * Start: S49 (all rows, positivity: the twisted rows m₁n₁ are rows of norm ≤ N(m₁)N(n₁)) gives κ₀ ≤ 2/3.
  * Δ_n follows the untwisted map of S52_LIMIT_PROOF, so c_n → c* ≤ 1.57.

## 3. Summation over boxes — PROVED
* **One slice.** For a slice, every dyadic length Y ≤ L₀ (with Y′ = Y/N(d) ≤ Y) satisfies S^ϑ/Y ≤ L₀·F^{c/5.5}·Q^{O(Δ+δ′)}:
  * if Y′ ≤ Q₁^{3/2} (q ≥ ⅔), then S^ϑ/Y ≤ Q₁Y′^{1/3}F^{c/5.5} ≤ Q₁^{3/2}F^{c/5.5};
  * otherwise S^ϑ/Y ≤ Y′F^{c/5.5} ≤ L₀F^{c/5.5};
  * and Q₁^{3/2} ≤ L₀.
* **Fullest boxes.** With X₁ ≤ Q/(X₂²X₃³X₄⁴), L₀ ≤ (X₁³X₂²X₃)^{1/2} and F ≤ X₂²X₃X₄², a box contributes ≪ Q^{3/2+O(δ′)} times
  X₂^{−1+0.364c} · X₃^{−3+0.18c} · X₄^{−5+0.364c}.
  This already includes the counts X₂X₃X₄ of (m₂, m₃, m₄) and the τ(m₄) choices of d.
* **Convergence.** For c ≤ 1.57 all three exponents are negative. So the dyadic sum over (X₂, X₃, X₄) is O(1), and the boxes add
  O(Q^{3/2+ε}) in total. Boxes with smaller X₁ only shrink L₀.

## 4. Theorem — PROVED, given the framework (F1)–(F8)
For every ε > 0: Σ_{N(m)≤Q}|ψ_α(m³,½+it)|² ≪_ε Q^{3/2+ε}(1+|t|)^{O(1)}, over all rows m (all four box types).

* **Squarefree rows.** S52_LIMIT_PROOF (main box, the binding case).
* **All other rows.** §1–§3 above, with exponents to spare. S49's direct box bounds are not needed, though they also cover every box with
  2x₂ + 4x₃ + 5.5y ≥ 1/8.
* **Still open** (unchanged): the constants in Lemma 4 (non-coprime and cusp bookkeeping, now also for twisted series, where the twist only
  multiplies the local factors), and the T↔T* transfer (S52 §9: a convergent e-sum, costing a constant).

## 5. The T↔T* transfer — PROVED
* **Setting.** Fix (a) runs the recursion on the CG-normalised family T*, where T = T* ∗ λ, λ is supported on fourth powers, and
  |λ(e⁴)| ≍ N(e) after normalisation (S51 §20). All tables in S52_LIMIT_PROOF are T*-tables.
* **S49 into T*.** T* = T ∗ λ′, with λ′ the coefficients of L(4s−3, α⁴), also supported on fourth powers with |λ′(e⁴)| ≍ N(e). Hence
  ‖R*_Y‖ ≤ Σ_e N(e)·‖R_{Y/N(e)⁴}‖, using χ̄_{e⁴}(m) = 1_{(e,m)=1} and positivity in the rows. The S49 bound decays at least like
  N(e)^{−8/3} in this norm, so the e-sum converges and the S49 inputs (F8) hold for T* up to a constant.
* **Result back to T.** ‖R_{T,Z}‖ ≤ Σ_e N(e)·‖R*_{Z/N(e)⁴}‖ ≤ Σ_e N(e)·(Q·(Z/N(e)⁴)^{4/3+ε})^{1/2} = Q^{1/2}Z^{2/3+ε}·Σ_e N(e)^{−5/3}.
  This uses the T*-table at q ≥ ⅔: the length Z/N(e)⁴ is below critical. So S_T ≪ QZ^{4/3+ε} = Z^{2+ε} at Q = Z^{2/3}. The diagonal
  stays ≍ Z, since Σ_e N(e)²/N(e)⁴ < ∞.

## 6. What remains
* **Lemma 4** (S52_PROOF): non-coprime local types and cusp components. The argument is written (a)–(f). What the proofs above actually
  use is narrower:
  * the local factors at y′ = 0 in the chambers s_c and dcdc, which are Π_{π|ν}(1 + O(1)) ≪ N(ν)^ε;
  * positivity for sublattice restrictions in B1;
  * finitely many cusp components.
  A line-by-line check of (f) (the transport of the local types through the termwise FEs) is the one item not yet done independently.
* Everything else in the chain (L³) ≪ Q^{3/2+ε} is proved above, or is a published input (BGL, GL, Huxley, DDHL (4.28), Bochner).
