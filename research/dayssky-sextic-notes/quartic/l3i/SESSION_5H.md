# 5-hour (L³) session, 16 Oct 2026, 07:27Z → 12:27Z (DaysSky with Claude assistance)

Goal set by the user: a proof of (L³), Σ_{N(m)≤Q}|ψ_α(m³,½+it)|² ≪ Q^{1+ε}. Labels: PROVED / NUMERICAL / HEURISTIC / OPEN.

## 07:28 Numerics: is there structure in the off-diagonal? (l3/scale*.txt, multi*.txt) — NUMERICAL
* S/diag − 1 over 30+ windows (rows 62 … 12961) has the scale of a mean of independent Exp(1) rows, ≈ rows^{−1/2}.
* So S − diag ≈ Q^{1/2}Z, and no secondary main term is visible. There is no hidden identity to exploit; a proof needs genuine
  equidistribution.

## 07:31 Literature — the analogue is open
* For GL(3) twisted by Dirichlet characters, a 2026 paper (Math. Ann., "Strong hybrid subconvexity for twisted selfdual GL3
  L-functions") states that the family second moment is known only at (qT)^{3/2+ε} unconditionally (Gallagher's hybrid large sieve).
  It is (qT)^{1+ε} only under generalised Lindelöf.
* (L³) is exactly this shape: family Q, conductor Q³, large-sieve value Q^{3/2}.

## 07:33 Rankin–Selberg on the G₂ cover — idea checked, does not avoid the floor (HEURISTIC analysis)
* **The cube structure is a cocharacter.** Index (m³,1) = t_m·(1,1) with t_m = (m⁶,m⁹), since α₁(t_m) = m³ and α₂(t_m) = 1. So sums
  over m are Mellin transforms along a cocharacter, which is the domain of Rankin–Selberg integrals.
* **But unfolding does not isolate the d = 1 slice.** ∫^{reg}|E_{G₂}|²E_{P}(w) over P (Levi ∋ α₂, 3-step unipotent radical) gives
  Σ_{a₁,a₂}|W_{(a₁,a₂)}(E)|²·(weights). These are second moments of G₂ first moments over twisting parameters, i.e. S50's Off.
  Off enters only after S50's Cauchy step, S² ≤ Z(ZS + Off), whose floor S ≤ Z² is Q^{3/2}.
* **The Cauchy-free form.** S − diag = (Q/Z)·Σ_h X(h)P(h), with X(h) S50's first moment and P(h) = Σ_c μ(c)ξ_h(c)W, from
  T(c)g̃(h;χ̄_c) = μ(c)χ_c(h). Bounding ⟨X,P⟩ needs correlation information for 1/L(·,ξ_h) across the quartic Hecke family
  (ratios-conjecture level). GRH alone gives only the large-sieve shape (L3_REPORT §10b).

## 07:36 What (L³) is equivalent to: an asymptotic large sieve for μ over a thin family (PROVED reformulation, HEURISTIC sizes)
* **Dual identity (BELOW_3_2 §1).**
  S − c₀QZ = (Q/Z)·Σ_{N(k)≤H′}(|P(k)|² − C_kZ), with P(k) = Σ_c μ(c)ξ_k(c)W(N(c)/Z) and H′ = Z²/Q = Z^{4/3}.
* **So (L³) needs** avg_k |P(k)|² = C_kZ + O(Q) = C_kZ·(1 + O(Z^{−1/3})): a power-saving asymptotic for the variance of μ over the thin
  family {ξ_k}, which has one quartic character per modulus.
* **Contrast with the full family.** For all Dirichlet characters mod q ≤ Q′ with Q′ ≥ N, the same variance statement is trivial:
  orthogonality kills every off-diagonal pair with q > N, and the rest is ≤ N². The thin family has no such orthogonality. Its
  off-diagonal kernel is the incomplete sum K(n) = Σ_{N(k)≤H′}(k/n)₄ at n = cc′³, of length H′ = Z^{4/3} against modulus ~Z².
* **Random model.** K(n) ≈ √H′, so the off-diagonal per k is ≈ Z·H′^{−1/2} = Z^{1/3}, far below the allowed Q = Z^{2/3}. This matches the
  numerics (07:30).
* **Why GRH is not enough.** GRH gives P(k) ≪ Z^{1/2+ε} pointwise, hence Σ_k|P(k)|² ≪ H′Z^{1+ε}. That yields S − c₀QZ ≪ Z^{2+ε} and
  Q^{3/2}, but not the needed precision |P|² = C_kZ + O(Q) on average, which is a ratios-conjecture-type statement.
* **Quadratic analogue, for contrast.** The quadratic analogue Σ_{m≤N^{2/3}}|Σ_{n≤N}μ(n)(m/n)|² ≪ N^{5/3+ε} follows from GRH, because
  there the original side is itself a Möbius–character sum (quadratic Gauss sums are trivial). In the quartic case the original side
  carries quartic Gauss-sum phases (T(c) = μ(c)χ_c(−1)g̃(χ_c)), so GRH applies only on the dual side, where it falls short by Q^{1/2}.
* **Conclusion.** (L³) ⟺ an asymptotic quartic large sieve with power saving Z^{−1/3} for the Möbius function. No asymptotic large
  sieve is known for any thin character family (Conrey–Iwaniec–Soundararajan's needs the full family plus a t-average).

## 07:50 Self-contained statement of the proved bound: l3i/THEOREM_Q32.md
* The chain is steps 0–10 with a status for each, the published inputs, and the STANDARD items.
* What remains before external use is an independent end-to-end check.

## 08:10 S53 (R1) re-read (DDHL (4.12)–(4.21)) — confirmed factor by factor; recorded in S53 §4.

## 08:12 Lemma D cross-check (l3f/limit_map.py, finer grids) — NUMERICAL, consistent
* The grid map Δ ↦ Δ′ decreases toward the closed form κ/(5.5 + 12κ) as the grid is refined:
  * Δ = 0.0025: 0.000712 → 0.000625 → 0.000521, with formula 0.000470;
  * Δ = 0.0005: 0.000417 → 0.000208 → 0.000130, with formula 0.000109.
* The grid is a discretised upper bound and the closed form is its continuum limit.

## 08:20 Literature: asymptotic large sieve
* Conrey–Iwaniec–Soundararajan (arXiv 1105.1176) treat the full family of primitive characters with a t-average, and tie Möbius-weighted
  versions to RH.
* Nothing found for thin (one character per modulus) families. That is what (L³) needs (07:36 entry).

## 07:50–08:05 Direct attempts at Q^{1+ε} after the user's redirect (all recorded with the exact point of failure)
1. **The sextic "canonical moment" in the openai baseline.** Moments/Canonical.lean works with complete prime-power correlations, in
   the regime where row completion and the large sieve suffice. It is not a thin-family Lindelöf statement. The quartic route needs
   the thin regime because it aims at θ = 5/6 (L3_REPORT §1). Nothing to transfer.
2. **Suzuki's τ₄(a²) = g̃₄(a)N(a)^{−1/4}.** This makes T a family of theta coefficients at squares. Rankin–Selberg against |θ|² sums
   over all indices, not over the thin family χ_m, so it gives no Lindelöf for the thin family.
3. **Type II with the quadratic coupling.** T(ab) = T(a)T(b)(a/b)₂s(a,b).
   * Cauchy over a turns the coupling into an extra family dimension: the characters (m³a²/·)₄, a family of size QA.
   * That gives Σ_m|R_II|² ≪ QZA + Z²/A under the optimal sparse quartic large sieve. At A = Z^{1/6} this is QZ·Z^{1/6}: better than the
     floor, but still a loss of Z^{1/6}.
   * The loss is exactly the a-phase cancellation that Cauchy discards. Expanding without Cauchy returns to the dual identity: after
     Poisson in m, (a/b)₂ squares to 1.
4. **Any power saving below 3/2.** This needs an asymptotic large sieve for μ over the thin quartic family with *some* power saving δ
   (relative error Z^{−δ} on the dual variance). Quantities checked:
   * individual incomplete-sum bounds (Pólya–Vinogradov, Burgess over Z[i]) give ≈ Z^{−1/3} per pair, i.e. Z^{5/3} in total against
     the Z^{1−δ} needed;
   * GRH, or a quasi-RH at θ, gives pointwise P(h) ≪ Z^{θ} but no variance asymptotic;
   * pointwise subconvexity for ψ_α(m³,½) would need exponent < 1/4, against convexity 3/4.
5. **Second-moment Bochner with mean-value regions** (second.py, plus σ₃-images of the Möbius side). Unconditionally the Möbius side
   has only Re > 1 regions. The convex combinations give u₃ = (1 + 2.5t)/(1+t) ≥ 1 at the centre, i.e. no gain below Q^{3/2}.
Conclusion: no route to Q^{1+ε} (or to Q^{3/2−δ}) found in this session. Every mechanism tried either reduces to the self-dual
identity or needs an asymptotic large sieve for μ over a thin family, which is not known for any thin family.
