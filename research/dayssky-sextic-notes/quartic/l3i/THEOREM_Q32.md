# Theorem: the quartic family second moment at Q^{3/2+ε}, with its complete dependency list

DaysSky with Claude assistance, 16 Oct 2026.

Labels:
* **PROVED** means the argument is written in the cited file.
* **EXACT** means checked by exact computer arithmetic (script named).
* **PUBLISHED** means a cited theorem from the literature.
* **STANDARD** means routine and not written out.

## Statement
Let ψ_α(ν, s) = Σ_{c primary} g̃(ν; χ_c) α(c) N(c)^{−s} be the α-twisted quartic Kubota series over Z[i], with α(c) = c/|c|. For every ε > 0:

  Σ_{N(m)≤Q} |ψ_α(m³, ½+it)|² ≪_ε Q^{3/2+ε}(1+|t|)^{O(1)}.

* **Previous record:** Q^{13/8+ε} (S49, unbalanced quartic large sieve).
* **Lindelöf on average:** would be Q^{1+ε}, i.e. (L³), which stays OPEN.

## The proof, as a chain

| step | content | where | status |
|---|---|---|---|
| 0 | AFE: (L³) is the row mean square S(Q, Z) of R(m) = Σ_c T(c)χ̄_c(m)W(N(c)/Z) at Z = Q^{3/2}; all rows m = m₁m₂²m₃³m₄⁴ | S49, S52_FULLSUM §1 | PROVED |
| 1 | S49: S(q) ≤ s49(q) = 1 + max(1, q + ⅓, q/2 + ¾), via the unbalanced quartic large sieve (dFDH Conj. 1.3 for max ≥ min^{6/5}) | l3c/S49_PROOF.md | PROVED (from Gao–Zhao, BGL: PUBLISHED) |
| 2 | S50: Cauchy step plus Poisson in c. S² ≤ Z(ZS + Off), with Off a sum of first moments X(h₀e², j) over twists h₀ ≤ H = Q²/Z | l3d/S50_REDUCTION.md Steps 1–2; S52_AUDIT A1–A4 | PROVED |
| 3 | X(h₀) is the Mellin evaluation of a G₂ double Dirichlet series (n = 4, α-twist), with twist rule and Gamma normaliser from DDHL (4.28) | l3d/S51_G2.md §1–§2, §10, §17 | PROVED (DDHL: PUBLISHED) |
| 4 | That G₂ series exists over Q(i) with all 12 chamber functional equations: the Chinta–Gunnells construction, corrected for μ₈ ⊄ Q(i) | l3g/S53_QI_LOCAL.md §1–§4 | PROVED (local: EXACT, FLINT; global: from DDHL (4.28) + Remark 4.2 and Hecke equidistribution) |
| 5 | Comparison with the actual series: naive = Z_CG × E′ × normalisers⁻¹. E′ converges absolutely at the evaluation point (−1/12, ¼) and on (x, −3x), −1/12 < x < 0 | l3f/S52_CORRECTION_NAIVE_VS_CG.md | PROVED (exponents EXACT) |
| 6 | Mean-value Bochner: block bounds B1/B2, regions, h₀-cocycle (12 chambers) | l3f/S52_PROOF.md Lemmas 1–3, 5; cocycle.py, fcocycle.py | PROVED (cocycles EXACT) |
| 7 | Lemma 4 (non-coprime local factors) for the CG series at the points used | l3g/S53_QI_LOCAL.md §6 | exponents read from exact arithmetic (half-integers, so rounding is exact) |
| 8 | Limit mechanism: region membership, cone lemma, two-chamber hull point, bookkeeping map Δ′ = κ/(5.5 + 12κ), iteration from S49's Δ₀ = 1/12 to S(⅔) ≤ 2 + δ′/2 | l3f/S52_LIMIT_PROOF.md Lemmas A–D, Thm | PROVED (cone lemma EXACT, cone_exact.py) |
| 9 | All other rows (boxes) and the T ↔ T* transfer | l3f/S52_FULLSUM.md | PROVED |
| 10 | Conversion: (L³) ≪ Q^{1.5(S(⅔) − 1)} = Q^{3/2+O(δ′)} | S52_LIMIT_PROOF Thm step 4 | PROVED |

## Published inputs used
* David–Dunn–Hamieh–Lin (arXiv 2306.11875): (4.28) and Remark 4.2, from Diaconu, Invent. Math. 157 (2004), and Suzuki (1983). This is the
  Kubota functional equation over Q(i), n = 4, with ∞-type, as a 24-cusp vector.
* Gao–Zhao: the quartic large sieve with (MN)^{2/3}.
* Blomer–Goldmakher–Louvel: the quartic large sieve.
* Hecke: equidistribution of Grössencharacters on primes in ray classes.
* Bochner's tube theorem (Hörmander Thm 2.5.10), and the three-lines lemma.

## Items labelled STANDARD
* Polynomial growth (finite order) of the normalised G₂ series. This is Phragmén–Lindelöf with the functional equations, as in CG Thm 6.1
  and [BBF06].
* S50 Step 3's class and unit bookkeeping (signs on classes mod 16).

## What is not claimed
* (L³) at Q^{1+ε}. The bound Q^{3/2} is the floor of S50's Cauchy step (S² ≤ Z²S + …). It is also the value given by GRH and by the
  conjectured optimal quartic large sieve (L3_REPORT §10b, S38). Going below it needs an input that breaks the self-duality of
  BELOW_3_2.md §3. l3i/SESSION_5H.md identifies that input: an asymptotic large sieve, with power saving Z^{−1/3}, for μ over the thin family
  of quartic residue characters.
* Any change to θ (L3_REPORT §2).

## Before external use
An independent end-to-end check is still needed, in particular of:
* the reading of DDHL Remark 4.2 in S53 (R1), re-read at 08:10Z and confirmed factor by factor (S53 §4);
* the S50 Step 2 Poisson weights at non-squarefree and non-coprime moduli.
