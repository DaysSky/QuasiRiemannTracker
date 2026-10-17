# S52 proof notes: the analytic lemmas, stated precisely (DaysSky with Claude assistance)

Notation:
* T(z) = {z ∈ C^n : Re z ∈ U} is the tube over an open U ⊂ R^n.
* For an index set 𝒦 (twists h₀ ∼ K), ℓ²(𝒦) has the normalised norm ‖f‖² = |𝒦|^{−1}Σ_{h∈𝒦}|f_h|².

## Lemma 1 (Bochner with growth; vector-valued) — PROVED
Let U ⊂ R^n be a bounded connected open set, V a Banach space, and F : T(U) → V holomorphic with ‖F(z)‖ ≤ M for all z ∈ T(U).
Then F extends holomorphically to T(conv U), with ‖F‖ ≤ M there.

*Proof.*
1. For each continuous linear functional φ with ‖φ‖ ≤ 1, f = φ∘F is holomorphic on T(U) with |f| ≤ M.
2. By Bochner's tube theorem (Hörmander, Thm 2.5.10), f extends to T(conv U). For |λ| > M, the function 1/(f − λ) is holomorphic on
   T(U), so it also extends. Hence f ≠ λ on T(conv U), and |f| ≤ M there.
3. The extensions of φ∘F, over all φ, define a weakly holomorphic, hence holomorphic, V-valued extension.
4. ‖F‖ = sup_φ |φ∘F| ≤ M. ∎

**Corollary 1a (polynomial growth and affine exponents).** Suppose that on T(U)
  ‖F(z)‖ ≤ C(1+|Im z|)^A · N^{ℓ(Re z)}  for some affine ℓ : R^n → R.
Then the same bound, with (1+|Im z|)^{2A} and the same ℓ, holds on T(conv U).
* Apply Lemma 1 to G(z) = F(z)·Π_j(z_j + M₀)^{−A}·N^{−ℓ(z)}, where ℓ(z) is the complex-affine extension and M₀ > sup|Re z| + 1.
* Taking the infimum over all affine ℓ ≥ b on U gives the convex envelope of b, by LP duality. This is what the LP in hull*.py computes.

## Lemma 2 (block bounds) — PROVED from the cited inputs
Let a(c) be the fix (a) coefficients T*(c), and let the d-weights be g̃(h₀e²;χ_d)D^{−½−y} (|g̃| ≤ 1 for (d, h₀e) = 1; the non-coprime d
are handled as in Lemma 4). For the smooth dyadic block c ∼ C, d ∼ D with D ≤ (N(h₀e²)C)^{1/2+ε}:
* **(B1)** ( |𝒦|^{−1}Σ_{h₀∈𝒦} |Σ_{c∼C}Σ_{d∼D} …|² )^{1/2} ≪ C^{−x}D^{−y}(Ξ₄(KC,D)/K)^{1/2}(CD)^ε.
  Proof: Cauchy over c, then the quartic large sieve for the moduli h₀c, with e² a coefficient twist.
* **(B2)** Same left side ≪ D^{−y}C^{−½−x}(Ξ₄(K,D)/(KD))^{1/2}·D^{1/2}·S*(D,C)^{1/2}(CD)^ε.
  * S*(D,C) = Σ_{d∼D}|Σ_{c∼C}T*(c)χ̄_d(c)W|² is the T*-family second moment with rows D and length C.
  * Proof: write the double sum as Σ_d β_d χ̄_d(h₀) with β_d = g̃(χ_d)D^{−½−y}·(c-block), then apply the large sieve in h₀.
* Here Ξ₄(A,B) ≤ (A + B + (AB)^{2/3})(AB)^ε (Gao–Zhao; S49 Theorem 1 gives more in unbalanced ranges).
* **Inputs for S*(D,C).**
  * **(i) Convexity:** S* ≤ C·D^{5/2}, from |ψ̂(d³,½+it)| ≪ N(d)^{3/4}(1+|t|)^{O(1)} and the Mellin representation of the block.
  * **(ii) Reflection:** if C > N(d)^{3/2}, the c-FE of ψ̂ (fix (a), clean) gives S*(D,C) ≪ (C/C′)S*′(D,C′), with C′ = D³/C and S*′ the
    same bound for the dual family (cusp and twist variants). If C > D^{3+ε} the block is O(C^{−A}).
  * **(iii)** Any previously proved bound S*(Q,Y) ≪ Y^{s(q)} (S49, or an earlier stage of the recursion).

## Lemma 3 (regions) — PROVED from Lemma 2
* The identity-chamber function is F = (normalisers)·Z_{h₀}. Its representation summed over blocks converges absolutely and locally
  uniformly, with the ℓ²(𝒦)-norm ≪ K^{b(x,y;t_K,t_E)}·(1+|Im|)^{O(1)}, on Region_d ∩ {y ≥ 0} and on Region_c ∩ {x ≥ 0}.
* Region_d: every direction (t_C,t_D) = (1,δ), δ ∈ [0,½], has min(B1,B2) < 0. Region_c: the same with (ρ,1), ρ ∈ [0, 3/2].
* b is the maximum of min(B1,B2) over the finitely many blocks not covered by the asymptotic decay.
* The other 11 chambers are pulled back by the termwise functional equations (DDHL (4.28) with the normalisers of §1). Their h₀-factors
  N(h₀e²)^{φ_w(v)} come from the d-reflections.

## Lemma 4 (non-coprime and cusp bookkeeping) — argument given; constants not tracked
* **(a) Local types.** Let ν = h₀e²c. For π^j ∥ ν, the local Gauss sum g(ν;χ_{π^k}) is nonzero only if k = 0, or k = j+1, or
  (k ≤ j and 4 | k). The normalised sizes are |g̃| ≤ N(π)^{j/2} and N(π)^{k/2} respectively.
* **(b) Expansion.** By twisted multiplicativity,
  g̃(ν;χ_{δd′}) = g̃(ν;χ_δ)·g̃(νδ²;χ_{d′})·(unit), for d′ coprime to ν and δ | ν^∞ of allowed local type. Hence
  ψ(ν,w) = Σ_{δ allowed} w_δ N(δ)^{−w} ψ^{(ν)}(νδ², w), with |w_δ| ≤ N(δ)^{1/2}. Here ψ^{(ν)} sums over d′ coprime to ν, and its
  coefficients are χ̄_{d′}(νδ²)g̃(χ_{d′}). For Re w ≥ ½ (y ≥ 0) the δ-sum is a finite product over π | ν of (1 + O(N(π)^{−y}·…)), which is
  ≪ N(ν)^ε.
* **(c) Coupling.** δ depends on c only through which primes of rad(δ) divide c. Fix δ and restrict c to the sublattice rad(δ)·Z[i]. Then
  B1 is a large sieve over a subset of moduli (bounded by positivity) and B2 is unchanged. The δ-sum converges with the extra density
  factor N(rad δ)^{−1}.
* **(d) c-side terms.** Non-squarefree c and c sharing primes with d (where T*(c) or the reflected c-coefficients are supported) are
  handled the same way, with the c-side local types of S51 §20.
* **(e) Cusps.** The 24 components ψ_{j1} of DDHL (4.28) are Dirichlet series over d in fixed 2-adic classes, with 2-adic normalisations.
  Lemma 2 applies to each, uniformly.
* **(f) Other chambers.** Each reflected chamber's series has the same structure with its twist (ℓ_c, ℓ_d), so (a)–(e) apply verbatim.

## Lemma 5 (evaluation) — PROVED from Lemmas 1–3
For h₀ ∈ 𝒦 and the smooth first moment X(h₀e²; Q, Z) with c-twist e:
  ‖X‖_{ℓ²(𝒦)} ≪ Q^{3/4}Z^{½+x_r}K^{b_r} + Z^{½+x}Q^{½+y}K^{b}.
* (x_r, ¼) is any point of the hull, and (x,y) any hull point with 0 < y < ¼.
* Every normaliser argument must be ≥ 0 at both points (fix (a): 3x+y, 2x+y, x+y, y ≥ 0).

## Theorem (conditional on Lemma 4 only)
Σ_{N(m)≤Q}|ψ_α(m³,½+it)|² ≪ Q^{λ+ε}(1+|t|)^{O(1)}, with λ the fixed point of recur2.py. The LP is run on samples inside the regions,
so the λ obtained is a valid upper bound.
