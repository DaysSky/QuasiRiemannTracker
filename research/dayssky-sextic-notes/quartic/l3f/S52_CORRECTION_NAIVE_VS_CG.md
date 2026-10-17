# Correction (03:38Z, 16 Oct 2026): S52 needs the clean (Chinta–Gunnells) series in every non-identity chamber

DaysSky with Claude assistance. Labels as before.

## The gap — PROVED (exponent check, exact CG p-part exponents from l3e/coef_5.txt)
* **The naive series.** S50's first moment X(h₀) is a Mellin transform of the naive double series: squarefree, coprime (c, d), built
  from the actual (L³) coefficients. By S51 §20 it equals CG × Π_β 1/L_β × E′, where E′_p = (1 + Ax₁ + Bx₂)/N_p and N_p is the CG p-part.
* **The q-exponents of N_p's 13 non-constant monomials** are
  (0,1) ½, (1,0) ½, (1,2) 2, (4,1) 7/2, (4,4) 6, (5,2) 6, (6,2) 13/2, (5,4) 15/2, (6,5) 9, (9,4) 21/2, (9,6) 12, (10,5) 12,
  (10,6) 25/2.
* **The binding hull points of S52_LIMIT_PROOF Lemma C**, with E′ absolute convergence requiring every H-exponent < −1 and x, y > −½:
  * v₁ = (−1/3, 1): worst H-exponent −1.167. Inside.
  * the evaluation point (−1/12, ¼): worst −1.083. Inside.
  * v₂ = (1/3, −1): worst +1.0, and |Bx₂| = q^{1/2} > 1. Far outside.
* **No chamber-dcdc point is inside.** Every point of chamber dcdc's Region_d pulls back to identity y = −3x′ − 2y′ < −½, since
  x′ > ⅓ − y′/2. That is beyond the line |Bx₂| = 1, where even an Estermann-type continuation of E′ stops.
* **Consequence.** The naive series cannot be continued to the chamber-dcdc regions (expected: a natural boundary, since E′'s Euler
  factors are not cyclotomic). So S52_PROOF Lemma 4(f), "other chambers verbatim", is false for the naive series. S52's Bochner hull
  is valid only for a series whose single-variable slices are exact Kubota series in every chamber, namely the CG series.

## What survives, and on what condition
* **Corrected argument.**
  1. Run S52 (regions, Bochner with growth, the cone lemma, the two-chamber hull point) on the CG series Z_CG(·; h₀). Its chamber-w
     series is Z_CG at M_w v with the same p-parts. At the binding local point (x′, y′) = (⅓, 0), every CG monomial has exponent ≤ −1.33
     (worst (1,2)), so the non-coprime corrections converge absolutely there.
  2. Transfer back with naive = Z_CG × E′ × (normalisers)^{−1}. Use it only on the straight contour from the region of absolute
     convergence to the evaluation point (−1/12, ¼), and at the dual points (x, −3x) with −1/12 < x < 0. All of these lie in E′'s
     absolute-convergence region (checked: on 3x + y = 0 it converges exactly for x ∈ [−⅓, 0)). At p | h₀ the local factors of E′ are
     1 + O(p^{−δ}) at these points, so the twist costs N(h₀)^ε.
* **Status of the results.** S52_LIMIT_PROOF (Q^{3/2+ε}, squarefree rows) and S52_FULLSUM (all rows) hold **conditionally on task
  #26**: a G₂ WMDS for n = 4 over Q(i), with the Brylinski–Deligne sign and the α-twist, whose single-variable slices are exact Kubota
  series (CG's local-to-global property). The same applies to every earlier S52 bound that used non-identity chambers, including 1.5573.
* **The unconditional state of the art is unchanged:** S49's Q^{13/8+ε}.
* **What #26 still needs** (S51 §17–§19 did the twist rule, the Gamma normaliser, and the exact check of the local BD-signed action):
  * the slice property: for each d-exponent k₂, the CG p-part's x₁-polynomial equals the Kubota local factor of the c-slice, and
    symmetrically for the d-slices;
  * the global twisted-multiplicativity bookkeeping (class functions Ψ) over Q(i).

## Where naive and CG differ, and what the literature gives for #26 (03:42Z)
* **The c-slices agree — PROVED (exact p-part).**
  * D has no x₂¹ terms. So for d ∋ p exactly once, the CG coefficient of x₂ is H = N_{x₂¹}/D_{x₂⁰} = q^{1/2}(1 − q³x₁⁴)/(1 − q³x₁⁴) = q^{1/2},
    with no c ∋ p terms.
  * Hence CG's c-slice at squarefree d is Σ_{(c,d)=1} T*(c)χ̄_d(c)N(c)^{−s}: the naive series under fix (a).
* **The d-slices differ — PROVED (exact p-part).**
  * At c = p, H_{x₁¹} = (A + q²·x₂²)/(1 − q³x₂⁴): the full Kubota d-series at index h₀c, including d ∋ p² (local type j = 1, k = 2).
  * The naive coprime d-series lacks these terms. The ratio (1 + β_p x₂²)^{−1} has poles at |β_p x₂²| = 1, i.e. y = −¼ in identity
    coordinates.
  * Summed over c these poles are dense, so the naive series has a natural boundary at y = −¼, and the chamber-dcdc point (y = −1) lies
    beyond it.
* **The literature on #26.** Y. Chen, "Eisenstein series on metaplectic covers and multiple Dirichlet series" (arXiv 2411.13143v2,
  Jun 2025) proves the Brubaker–Bump–Friedberg Eisenstein conjecture: Whittaker coefficients of Borel Eisenstein series on n-fold covers
  are the CG WMDS. Its hypotheses are k ⊇ μ_{2n} and an adjoint metaplectic dual datum (G₂ qualifies). Its local theory uses 2n | (q−1),
  so that (π,π)_n = 1.
  * Q(i) with n = 4 lacks μ₈. At q ≡ 5 mod 8, (π,π)₄ = −1: exactly the primes where S51 §19 found raw CG fails and the
    Brylinski–Deligne-signed action passes.
  * So #26 = Chen's theorem extended to k ⊇ μ_n with the BD sign, plus the archimedean twist α and the h₀ twists. Supporting evidence:
    S51 §19 (the exact local relations hold with the BD sign) and DDHL (4.28) (the rank-1 Kubota functional equation over Q(i), n = 4,
    with ∞-type ℓ).

## #26 local theory settled; S51 §19 corrected (04:35Z; full write-up in l3g/S53_QI_LOCAL.md)
* **The correct local action over Q(i).** It is CG's Definition 3.1 with raw Gauss sums and Q_{β,i} multiplied by
  s^{μ(β)(B_i(β)+1)}, where B_i(β) = b_i + Σ_{j>i}c(j,i)b_j and s = (ϖ,ϖ)_S = (−1/ϖ)₄.
  * It comes from the self-symbol (η,ϖ^{k′})_S that CG's Lemma 5.9 drops when μ_{2n} ⊂ F.
  * PROVED by exact computer algebra in Q(√q, γ(1), x):
    * it is a W-action for all ℓ;
    * its p-parts are polynomial;
    * they satisfy the true one-variable slice functional equation in 4191/4191 classes.
* **S51 §19's BD-signed action** is a W-action but fails the slice functional equation (343/554). The two p-parts differ only by monomial
  signs, so every exponent statement in this note stands.
* **What #26 still needs.** Only the Ψ-bookkeeping of CG Prop 5.7 over Q(i), from DDHL (4.28). Self-symbols there have the form
  (y,−1)_S and preserve M(Ω,Φ). This item is labelled STANDARD, not written line by line.
* **05:30Z: the condition is discharged** (l3g/S53_QI_LOCAL.md §4). The Ψ-bookkeeping holds over Q(i) by a consistency argument from
  DDHL (4.28) and Remark 4.2 (published). So S52's Q^{3/2+ε} no longer depends on #26 beyond published inputs and STANDARD Bochner.
  The remaining items are those of S52_PROOF Lemma 4 (non-coprime and cusp constants).
