# S52 audit: every formula in the pipeline (recur7.py), with derivation and status

All logs are base Y, where Y is the c-length of the second moment being bounded. Q = Y^q, H = 2q−1, a twist h₀ has N(h₀) = Y^η, and the
coprimality piece e has N(e²) = Y^{tE}.

## A. S50 reduction at general q
* **A1. Cauchy.** S² ≤ Y·Σ_c W|Σ_m γ_mχ̄_c(m)|², with γ = R̄ on the row set.
* **A2. Poisson in c (S49 Lemma B).** The c-diagonal gives Y·S. The off-diagonal is (Y/Q)Σ_{N(h)≤H}w(h)|X(h)|², with Mellin-separated
  weights and X(h) = Σ_m γ_m g̃(χ̄_m)ψ_m(h).
* **A3. Writing h = h₀j⁴.** X(h₀j⁴) = Σ_{(m,j)=1}(…). Möbius over e | j gives pieces X_e(h₀) with twist h₀e², c-twist κ = e, and
  m-range Q/N(e). Their j-weights are J/N(e), with J = (H/N(h₀))^{1/4}.
* **A4. Bookkeeping.** S ≤ Y^s requires, for every twist level η, one of:
  * (MDS) (1 + (1−q) + η + (H−η)/4 + 2X_eff(η))/2 ≤ s, where X_eff is the MDS exponent for e ≤ Y^{ε₀};
  * (generic) 2 − q + (H−η)/4 + log Ξ₄(Y^η, Y^q) ≤ s, from the large sieve over h₀ with multiplicity J.
  * For the e-tail (e > Y^{ε₀}): 2 − q + (H−η)/4 − ε + log Ξ₄(Y^η, Y^{q−ε}) ≤ s, from the large sieve over h₀ for each e.
  * Here Ξ₄(A,B) = A + B + (AB)^{2/3} (Gao–Zhao).

## B. First-moment exponent X(η, ε)
* **B1. Main term.** X = Q^{3/4}Y^{½+x_r}N^{b}: the theta residue at (x_r, ¼).
* **B2. Dual term.** Y^{½+x}Q^{½+y}N^{b} at a hull point with 0 ≤ y < ¼.
* **B3. No shift.** The same as B2, at a hull point with y > ¼.
* **B4. Evaluation constraints.** All normaliser arguments must be ≥ 0. With fix (a): 3x+y, 2x+y, x+y, y ≥ 0.
* **B5. Hull and h₀-exponent.** The hull is the LP over pulled-back samples from all 12 chambers. b = local block exponent plus
  (c_w·v)(t_K + t_E), where c_w·v is the h₀-factor from the d-reflections.

## C. Local regions in chamber coordinates
* **C1. Block bounds B1, B2.** As in S52_PROOF Lemma 2: RMS over h₀, Ξ₄ in h₀, D-range t_D ≤ θ(t_K + t_E + t_C).
* **C2. σ = log S*(D,C).** The minimum of:
  * convexity, t_C + 2.5t_D;
  * the table, t_C·S(t_D/t_C);
  * negligible for t_C > 3t_D.
  The table at q < 2/3 is obtained by reflection: S = (1 − (3q−1)) + (3q−1)·S(q/(3q−1)).
* **C3. Region_d(θ).** Main directions δ ≤ θ, and FE-dual directions δ′ ≤ 1−θ evaluated at (x+y, −y). The dual-block twist exponent
  carries an extra −y(t_K + t_E).
* **C4. Polar chambers (ℓ_d = 0; 4 of 12).** The constraint x > 3/8 + θ(¼ − y), with polar exponent
  t_C(3/8 − x + θ(¼−y)) + t_K(θ(¼−y) − 1/8) + t_E(θ(¼−y) + 3/8) (Rankin–Selberg plus Cauchy).
* **C5. Region_c.** Directions ρ ∈ [0, 1.5], b = 0. Only used for x ≥ 0; its boundary has y > ½, so it meets no polar issue.
* **C6. Older regions, kept as extra samples.** R_d at y ≥ 0 (x > ½); R_d via the FE at y < 0 (x > ½ − y, b = −y(t_K + t_E) + t_E/6);
  R_c (y > Λ(x)/2).

## D. Recursion
S(q) for q ∈ [0.666, 1.2] is replaced by min(S49, A4 at q), and λ by 1.5(S(2/3) − 1). Each iteration uses only previously established
bounds.

## E. Fixes made in the 02:19Z audit (S52 §17)
* **E1. Bootstrap start.** The recursion must start from a proved bound, λ = 13/8 (S49). recur6 and recur7 started at 1.5253; recur9 onward
  start at 13/8.
* **E2. Table reads.** Outside the grid, S(q) is s49(q); np.interp's clamping is not allowed. Inside the grid only step-up reads are
  justified, by monotonicity in Q (recur11). recur9 and recur10 used linear reads.
* **E3. Polar term.** Polar term on FE-reflected R_d points whose image has y′ < ¼.
* **E4. φ_w path-independence.** Verified (cocycle.py).
* **E5. Limit mechanism.** Re-derived by hand and by fresh code (cone_check.py).

## Known open items
* **Lemma 4** (non-coprime and cusp bookkeeping): argument given, constants untracked.
* **Boxes and T↔T* transfer: resolved.** S52_FULLSUM.md §1–§5 handles twisted main-box problems, the F-cocycle (EXACT) and the
  transfer in both directions.
* **Independent re-implementation.** The pipeline has been patched three times (§10, §14, §16), plus E1–E3. The binding mechanism has
  now been re-derived independently (§17 and cone_check.py). The non-binding twists still rely on the recur* LP hulls (NUMERICAL).
* **Rigorous numerics: resolved for the limit statement.** S52_LIMIT_PROOF.md proves Q^{3/2+ε} with an exact cone lemma
  (cone_exact.py) and closed-form bookkeeping on the continuum q ∈ [⅔, 5/6]. The LP-based intermediate numbers (recur*) remain
  NUMERICAL.
