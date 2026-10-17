# S50. The large-sieve limit Q^{13/8} for (L³) is exactly the T–theta correlation

Setting and notation as in S49 (l3c/S49_PROOF.md). Q = Z^{2/3}, and R(m) = Σ_c T(c)χ̄_c(m)W(N(c)/Z) with T(c) = α(c)g̃(χ̄_c).
S = Σ_{m∈𝔅}|R(m)|² for the main box 𝔅 (squarefree m ∼ Q). This is the box that binds at Q^{13/8} in S49 Theorem 2: its exponent there is
1/2 + 3η/4, and every other box is strictly smaller.

## Step 1. Cauchy route (PROVED)
S = Σ_c T(c)W·[Σ_m R̄(m)χ̄_c(m)] ≤ ‖TW‖·‖A*R‖, so S² ≤ Z·Σ_c W|Σ_m γ_m χ̄_c(m)|² with γ_m = R̄(m).

## Step 2. Poisson in c with γ = R̄ (the identity behind S49 Lemma B; main box is the coprime-radical case written out there)
Σ_c W|Σ_m γ_mχ̄_c(m)|² = Z·Ŵ(0)Σ|γ_m|² + Off, with
  Off ≍ (Z/Q)·Σ_{0<N(h)≤H} |Σ_m γ_m g̃(χ̄_m) ψ_m(h)|²·(Mellin-separated weights),   H = Q^{1+ε}·Q/Z = Z^{1/3+ε}.
Here ψ_m(h) = (h/m)₄. The weight g̃(χ̄_m) is τ(ψ̄_m)/√N(m), and the cross factors are signs by reciprocity.
Write h = h₀j⁴ with h₀ fourth-power-free. Then ψ_m(h) = ψ_m(h₀)·1_{(m,j)=1}, and
  Off ≍ (Z/Q)·Σ_{h₀} Σ_{N(j)≤(H/N h₀)^{1/4}} |X(h₀,j)|²,   X(h₀,j) := Σ_{m∈𝔅,(m,j)=1} R̄(m) g̃(χ̄_m) ψ_m(h₀).
* For h₀ > K₀: the S49 Lemma A bound, applied with the j-multiplicity (H/K₀)^{1/4}, gives Off_gen ≪ Z(H/K₀)^{1/4}S + Z·S.
* Hence S ≪ Z^{2+ε}(1 + (H/K₀)^{1/4}) + (Z·Off_T)^{1/2}, where Off_T = (Z/Q)Σ_{N h₀≤K₀}(H/N h₀)^{1/4}max_j|X(h₀,j)|².

## Step 3. Double reflection of X (PROVED from the Kubota functional equations; constants and class/unit bookkeeping STANDARD)
|X| = |Σ_m R(m)·g̃(h₀;χ_m)| (conjugate, using conj g̃(χ̄_m) = χ_m(−1)g̃(χ_m)). Swap the sums:
  X = Σ_c T(c)W·Σ_m g̃(h₀c; χ_m)ν(m)V(N(m)/Q)   (reciprocity sign ν on classes mod 16).
* **m-reflection.** The inner sum is a smooth partial sum of the untwisted (ℓ = 0) Kubota series ψ(h₀c, w). Its functional equation
  (DDHL (4.28), index preserved) gives:
  * the theta-pole term Q^{3/4}V̂(3/4)·ρ(h₀c), where ρ(ν) = Res_{w=3/4} ψ(ν,w) is the quartic theta coefficient in the normalisation
    of dFDH §9: the "B^{3/4}τ₄" term;
  * plus (Q/√(N(h₀)Z))·Σ_{N(m')≲N(h₀)Z/Q} g̃(h₀c;χ_{m'})·(smooth weight).
* **c-reflection of the dual piece.** After swapping back, the c-sums are α-twisted Kubota sums R_w(m′) at index m′³, which are entire.
  Reflecting each gives length N(m′)³/Z, which is negligible for N(m′)³ < Z^{1−ε}.
  A large sieve (S49 Theorem 1) on the reflected rows bounds the rest:
  |X_dual| ≪ Z^{1+ε}N(h₀)^{1/2} for N(h₀) ≤ Z^{1/15}, and ≪ Z^{17/18+ε}N(h₀)^{4/3} beyond.
* **Result.** X(h₀) = Q^{3/4}V̂(3/4)·𝒞(h₀) + X_dual, where

    𝒞(h₀) := Σ_c T(c) W(N(c)/Z) ρ(h₀c)   (the T–theta correlation).

## Step 4. Exponent bookkeeping (l3d/reduction.py)
Suppose |𝒞(h₀)| ≪ Z^{π+ε}N(h₀)^{β}. Optimising K₀:

| input on 𝒞 | π | (L³) sum bound |
|---|---|---|
| Cauchy–Schwarz plus the Rankin–Selberg mean square Σ_{c∼Z}|ρ(h₀c)|² ≪ (N(h₀)Z)^{3/4} (dFDH §8) | 7/8 | **Q^{13/8}**, exactly S49's limit |
| power saving | 0.80 | Q^{1.597} |
| power saving | 0.75 | Q^{1.572} |
| square-root (P4) | ≤ 0.60 | Q^{1.543}; the cap comes from the crude X_dual for large h₀ |

**Conclusion (PROVED modulo the STANDARD items above).** The large-sieve limit is not a structural barrier for the T-family. It is exactly the
event that T correlates with the quartic theta coefficients at full Cauchy–Schwarz strength. This matches dFDH's extremal example, whose
column coefficients are Gauss sums with theta bias. Any power saving

  (P4)   Σ_{c∼Z} T(c) ρ(h₀c) ≪ Z^{7/8−δ}, uniformly for N(h₀) ≤ Z^{κ},

improves (L³) beyond Q^{13/8}, with no other new input.

## Step 5. Numerics for (P4) (l3/lin_stats.py, l3/lin_twist.py; T-family data phase_1e5/1e6.bin)
* With S ≈ QZ, the linear statistic |X(h₀)|/√S is O(1) in every window: 0.6–1.3 for h₀ = 1, and 0.2–1.8 for ten twists h₀ with N(h₀) ≤ 26.
  The cumulative partial sums in m show no growth.
* If 𝒞 were Cauchy-sized (Z^{7/8}), the ratio would be ≈ Q^{1/4}Z^{3/8}/… ≈ 10³ at Z = 10⁶.
* The data give |𝒞(1)| ≲ 3Z^{1/2}Q^{−1/4} ≈ Z^{0.41}, consistent with square-root size Z^{3/8}.
* Status: (P4) is NUMERICAL, with strong support.

## Step 6. What (P4) is, structurally (HEURISTIC / OPEN)
* Σ_c T(c)ρ(h₀c)N(c)^{−s} is the residue along w = 3/4 of the double series Ξ(s,w) = Σ_d g̃(h₀;χ_d)ψ_α(d³,s)N(d)^{−w}.
  Its partial sums at length Q are exactly X.
* In BBF orientation, Ξ pairs Gauss sums of opposite orientation on the two variables: χ on d, χ̄ (power 3) on c, coupling power 3.
  A Weyl-group multiple Dirichlet series needs a W-invariant quadratic form, which means equal orientation (Q(α₁) = Q(α₂)).
  So Ξ is not a standard A2 WMDS. This corrects the HEURISTIC identification in S44: the orientation there should be re-derived.
* The analogous same-orientation object Σ_c g(χ_c)τ(c)N(c)^{−s} is a Whittaker coefficient of the GL3-cover Eisenstein series induced from
  theta. That object has continuation by Langlands–Moeglin–Waldspurger theory for covers. Ours is its mixed-orientation twin, for which
  no continuation is known.
* (P4) is therefore a new, precisely stated input. It is a correlation between the α-twisted Gauss-sum family and the quartic theta
  coefficients, whose values at primes are unknown (dFDH: "a central open problem"). It is weaker than (L³): it only needs a power saving
  over Cauchy–Schwarz, and numerics show square-root size.

## Step 7. What is known about ρ = τ₄, and why it does not close (P4) (lit/dfdh.txt §6.2)
* Suzuki [Suz83] (dFDH (6.14)–(6.16)), for m ≡ 1 mod λ³ squarefree and coprime to δ:
  * τ₄(m⁴δ) = τ₄(δ);
  * τ₄(m³δ) = 0 for m ≠ 1;
  * τ₄(m²δ) = ±g₄(δ,m)N(m)^{−3/4}τ₄(δ; mu), an explicit Gauss sum;
  * τ₄(mδ) = 0 unless g₂(mδ/m₁, m₁) = N(m₁)^{1/2} for all m₁ | m (a quadratic-residue support condition, density 2^{−ω(m)}).
* **Exponent 3 or 2 would close (P4).** If the theta index were h₀c³, 𝒞 would vanish identically. If it were h₀c², 𝒞 would reduce to
  Σ_c α(c)χ_c(−h₀)μ²(c)N(c)^{−1/4}, a Hecke-character sum with square-root cancellation.
* **But the exponent is forced to be 1.** The dual Gauss sum of the column character c ↦ χ̄_c(m) is g(χ̄_m), and
  g(χ̄_m)χ_m(h₀c) = g(h₀c; χ̄_m) has index exponent 1. In the conjugate formulation it is again exponent 1, at h̄₀c̄.
  So only Suzuki's support condition applies, and it gives no cancellation.
* **Heuristic picture.** Patterson–Eckhardt conjecture: N(π)^{1/4}τ₄(π)² ∝ g̃₄(π). Since T(π)² = −α(π), both T(π) and τ₄(π) are
  "square roots" of Gauss-sum quantities. 𝒞(h₀) asks whether their relative signs are equidistributed: a quartic analogue of Heath-Brown–Patterson.
  * In the cubic case (Patterson: τ₃ ∝ ḡ₃) the same-orientation pairing has a main term (pole) and the opposite one cancels.
  * Ours is the opposite-orientation pairing (Step 6), so cancellation is the expected behaviour. The numerics agree (Step 5).
* **Averaging over h₀ (PROVED, exponent algebra).** Use the large sieve in h₀ ∼ K: Σ_{h₀∼K}|X(h₀;Q₂)|² ≤ Ξ(K,Q₂)·S(Q₂,Z) at Q₂ = Z^{1.2},
  where X_dual vanishes. This gives Σ_{h₀∼K}|𝒞(h₀)|² ≪ Z^{1.8+ε}, which beats Cauchy–Schwarz on average once K ≥ Z^{0.03}.
  But the bounded twists h₀ carry the j⁴ multiplicity H^{1/4} and by themselves reproduce Q^{13/8}. So (P4) is needed pointwise
  for N(h₀) = O(1), above all h₀ = 1.
* Positivity of the smooth majorant W means the j⁴ main term ∝ ∫W(x)|x|^{a−2}dx > 0 cannot be made to vanish by choosing W.

**Status of the (L³) programme after S50.**
* Unconditional: Q^{13/8} (S49).
* Q^{1.543}: conditional on (P4), the T–theta cancellation at bounded twists.
* Q^{3/2}: needs (P4) plus square-root cancellation in the dual m′-sums.
* Q^{1+ε}: needs Lindelöf on average for the T-family (S46), unchanged.
(P4) is a single, explicit, numerically supported statement: square-root cancellation in Σ_c α(c)g̃(χ̄_c)τ₄(h₀c).
It concerns the quartic theta function itself, not the T-family's second moment.

## Step 8. Refinement of the orientation analysis (23:20Z)
* **Identity (NUMERICAL check of a known identity).** T(c) = α(c)g̃(χ̄_c) = μ(c)·conj g̃(χ̄_c) = μ(c)χ_c(−1)g̃(χ_c). Brute force on all 66 squarefree
  primary c with N(c) ≤ 300 gives error 5e-15. So the T-family is the ℓ = 1 (α-twisted) object on the inverse-cocycle cover (χ̄), and
  equivalently the Möbius twist of the untwisted χ-cover Gauss sums. This is why ψ_α is entire.
* **The two formulations of X agree.** Tracking orientations (conj τ₄(δ) = τ₄^{χ̄}(−δ) and τ₄(δ̄) = τ₄^{χ̄}(δ), both from
  conj(x/m)₄ = (x̄/m̄)₄), the theta term is
  𝒞(h₀) ∝ Σ_c T(c)τ₄(−h₀c) = conj Σ_c T̄(c)τ₄(h̄₀c̄), with T̄(c) = ᾱ(c)χ_c(−1)g̃(χ_c).
  So it is either:
  * Gauss sums on the χ̄-cover against theta on the χ-cover at index c; or
  * Gauss sums on the χ-cover against theta at the conjugate index c̄.
* **Not a GL3-cover Eisenstein coefficient.** The BBF residue Σ_c λ(c)g̃(χ_c)τ₄(c) (Hecke character λ, same cover, same index) is a Whittaker
  coefficient of the GL3-cover Eisenstein series induced from θ ⊗ λ, and has continuation by Langlands–Moeglin–Waldspurger.
  Our pairing differs by index conjugation c ↦ c̄, which is not a Hecke twist. Brylinski–Deligne covers of GL3 need a W-invariant Q,
  so no single cover realises it.
* **Möbius form.** 𝒞(h₀) = Σ_c μ(c)·a(c) with a(c) = χ_c(−1)g̃(χ_c)τ₄(−h₀c). Heath-Brown–Patterson-style Type I/II decompositions stall
  at Type II, because τ₄(h₀de) has no factorisation in (d,e). Suzuki's relations factor only through fourth, third and second powers.
* **Status.** (P4) remains OPEN, precisely located: equidistribution of T(c) against theta coefficients across the two covers.
  It is numerically supported (Step 5).
