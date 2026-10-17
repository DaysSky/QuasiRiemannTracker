# Theorem ZF — complete proof

DaysSky with Claude assistance, 16 Oct 2026. This file contains the full argument behind l3j/THEOREM_ZF.md (statement, corollaries,
numerics and context are there).

Notation and conventions are as in PROBLEM.md:
* primary c ≡ 1 (mod 2+2i);
* χ_c = (·/c)₄;
* g(r; χ) = Σ_{x mod c} χ(x)e(Im(rx/c)), and g̃ = g/√N(c);
* T(c) = α(c)g̃(χ̄_c), with α(c) = c/|c|;
* λ(c) = χ_c(−1) = (−1)^{(N(c)−1)/4}, and e(c) ∈ {0,1} is defined by λ(c) = (−1)^{e(c)}.

**Hypothesis (G)^{tw}_{θ′}.** This is the squarefree, large-sieve form used throughout the project (PROBLEM.md (G)). For every
ε > 0 and j ≥ 0 there is C such that for all Z, Q ≥ 1 and h ∈ Z[i]∖0:

  Σ_{m∈Z[i]} V(m/√Q) |Σ_{c sqfree primary} T(c)W_j(N(c)/Z)χ̄_c(mh)|² ≤ C·(j+1)^{A+2}·Z·Q^{1+θ′+ε}·N(h)^ε.

**Theorem.** If (G)^{tw}_{θ′} holds with 0 ≤ θ′ < ½, then L(s, (k₀/·)₄λ^e) ≠ 0 for Re s > σ(θ′) = (7+10θ′)/(8+8θ′), for all k₀ ∈ Z[i]∖0
and e ∈ {0,1}.

**Relation to (L³).** By Mellin (Lemma 1), (G)^{tw}_{θ′} follows from the Dirichlet-series form
Σ_{N(m)≤Q}|ψ^{sf}_α(m³h, ½+it)|² ≪ Q^{1+θ′+ε}N(h)^ε(1+|t|)^A, where ψ^{sf} sums over squarefree c only.

* PROBLEM.md's ψ_α sums over all primary c. The two differ by local Euler-type factors at the primes dividing m: a 2×2 system per
  prime, see Lemma 1, step 4.
* (L³)^{tw} for ψ_α implies the ψ^{sf} form wherever those local matrices are boundedly invertible on the line. Checking that is a
  finite local computation, NOT done here.
* So the theorem is stated for (G)^{tw}, the form the project actually uses.

## 0. Fixed choices

* **Row weight.** B = 2. V is the Bessel potential on C ≅ R² with V̂(ξ) = (1+|ξ|²)^{−B}. It is positive, continuous and radial,
  V(x) ≪ e^{−|x|}, and V̂ ≥ 0.
* **Class.** Fix κ ∈ {0,1} and let 𝒞 = 𝒞_κ be the squarefree primary c with e(c) = κ.
* **Narrow annulus.** Fix 0 < η ≤ η₀(B), with η₀ chosen in §3. W is a smooth bump supported in [1, e^η], and W_j(x) := W(x)·(log x/η)^j
  for j ≥ 0. Write u_c := log(N(c)/Z) ∈ [0, η] on the support.
* **Rows.** Q = Z^x with x = 1/(1+θ′) ∈ (2/3, 1].

## 1. T-side bound

For a weight 𝒲 put R_𝒲(m) := Σ_{c∈𝒞} T(c)𝒲(N(c)/Z)χ̄_c(m) and S(𝒲) := Σ_{m∈Z[i]} V(m/√Q)|R_𝒲(m)|².

**Lemma 1.** Under (G)^{tw}_{θ′}, S(W_j) ≪_{η,ε} (j+1)^{A+2}·Z·Q^{1+θ′+ε}. The same holds for the gcd-twisted rows
Σ_m V(m/√Q)|Σ_{a} T(a)W_j(N(ha)/Z)χ̄_a(mh²)|², with Z replaced by Z/N(h), uniformly in h up to N(h)^ε.

*Proof.*
* **Direct part.** Both bounds are direct instances of (G)^{tw}. The twisted rows are the index twist h²; their weight W_j(N(h)·/Z) is
  W_j at scale Z/N(h). The class restriction is handled by step 1 below.
* **From the Dirichlet-series form.** Steps 2–7 record how (G)^{tw} would follow from the ψ^{sf}/ψ_α forms.
* **Caveat on step 4.**
  * For each prime π | m the relation between the states m³h and m³hπ² is a 2×2 system M_π(s).
  * Its entries from e ≡ 0 (mod 4) are Ramanujan-sum terms of size O(1) on Re s = ½, not small. So bounded invertibility of M_π on
    the line is not automatic.
  * This is exactly why the theorem is stated for (G)^{tw}.

Steps:
1. **Class restriction.** 1_{e(c)=κ} = ½(1 + (−1)^κλ(c)), and λ(c)χ̄_c(m) = χ̄_c(−m). So the class-restricted R is ½(R(m) ± R(−m)) for the
   unrestricted family, and it suffices to bound the latter.
2. **Mellin.** R_{W_j}(m) = (2πi)^{−1}∫_{(½)} ψ^{sf}(m³,s)Z^sW̃_j(s)ds. Here ψ^{sf}(r,s) := Σ_{c sqfree primary} α(c)g̃(χ̄_c, r)N(c)^{−s},
   using T(c)χ̄_c(m) = α(c)g̃(χ̄_c, m³).
3. **Cauchy–Schwarz.** |R_{W_j}(m)|² ≤ Z·(∫|W̃_j(½+it)|dt)·∫|ψ^{sf}(m³,½+it)|²|W̃_j(½+it)|dt, and ∫|W̃_j(½+it)|(1+|t|)^A dt ≪_η (j+1)^{A+2}.
4. **Local inversion.** Write a primary c as c_sf·d, with d the powerful part, coprime to c_sf. By CRT and quartic reciprocity:
   * α(c)g̃(χ̄_c, m³) = α(d)g̃(χ̄_d, m³)·σ(c_sf,d)·α(c_sf)g̃(χ̄_{c_sf}, m³d²), with a sign σ ∈ {±1} that depends only on classes mod 16.
   * For π^e ∥ d with e ≥ 2, χ̄_d is induced from χ̄_π^{e mod 4} at π, so g(χ̄_{π^e}, m³) = 0 unless π^{e−1} | m³. When it is nonzero,
     |g̃(χ̄_{π^e}, m³)| ≤ N(π)^{e/2} (Ramanujan-sum case e ≡ 0 mod 4) and ≤ N(π)^{(e−1)/2} otherwise. Hence
     |α(d)g̃(χ̄_d, m³)|·N(d)^{−½} ≤ 1.
   * Moreover χ_{c_sf}(d²) depends on d only through d² modulo 4th powers and through rad(d), the latter via coprimality.
   * So for each prime π | m the relation between ψ^{sf}(m³h, s) and ψ_α(m³h′, s), for h, h′ ranging over a finite set of local states at
     π, is a finite triangular system. Its diagonal is 1 and its entries are bounded on Re s = ½.
   * Inverting prime by prime gives ψ^{sf}(m³,s) = Σ_{h∈ℋ(m)} b_h(m,s)ψ_α(m³h,s) with rad(h) | rad(m), #ℋ(m) ≤ τ(m)^{O(1)} and
     |b_h(m,½+it)| ≤ τ(m)^{O(1)}.
5. **Summing.** Use Cauchy–Schwarz over h ∈ ℋ(m). Group the m by h and write m = rad(h)m′. Then (L³)^{tw} at index twist
   rad(h)³h gives Σ_{m′} ≪ (Q/N(rad h))^{1+θ′+ε}N(h)^ε. Since Σ_r τ(r)^{O(1)}/N(r) ≪ (log Q)^{O(1)}, the sum over h costs only Q^ε.
6. **Tails and growth.**
   * V decays exponentially, so the rows N(m) > Q·(log Q)² contribute O(Z^{−10}).
   * Shifting the contour to Re s = ½ needs ψ_α(r,·) entire with polynomial growth in vertical strips. This holds because the
     α-twisted Kubota series has no poles (Kubota; DDHL (4.28)).
7. **Twisted rows.** The same argument applies verbatim, since the index there is (mh²)³ = m³h⁶. ∎

## 2. Exact dual form for coprime pairs

**Lemma 2.** For c₁ ≠ c₂ in 𝒞 with (c₁,c₂) = 1 and any weights 𝒲₁, 𝒲₂:

  T(c₁)𝒲₁·conj(T(c₂)𝒲₂)·Σ_m V(m/√Q)χ̄_{c₁}(m)χ_{c₂}(m) = Q·Σ_{k∈Z[i]} a_{c₁}(k)·conj a_{c₂}(k)·(1 + QN(k)/(N(c₁)N(c₂)))^{−B},

with a_c(k) := μ(c)χ_c(k)𝒲(N(c)/Z)/√N(c).

*Proof.*
1. **Poisson.** Z[i] is self-dual for e(Im ·). For φ periodic mod q = c₁c₂:

   Σ_m V(m/√Q)φ(m) = (Q/N(q))Σ_k φ̂(k)V̂(√Q·k/q),  with φ̂(k) = Σ_{x mod q} φ(x)e(Im(kx/q)).

2. **CRT.** φ̂(k) = g(k; χ̄_{c₁}χ_{c₂}) = χ̄_{c₁}(c₂)χ_{c₂}(c₁)·g(k;χ̄_{c₁})·g(k;χ_{c₂}).
3. **The two Gauss-sum identities.**
   * T(c)g(k;χ̄_c) = μ(c)χ_c(k)√N(c) (L3_REPORT §10, checked numerically there).
   * g(k;χ) = χ(−1)·conj g(k;χ̄).
4. **Sign.** The leftover sign is χ̄_{c₁}(c₂)χ_{c₂}(c₁)·χ_{c₂}(−1) = (−1)^{e₁e₂}·(−1)^{e₂}, by quartic reciprocity for primary elements.
   * This equals 1 when e₁ = e₂ = κ.
   * Checked by brute force: l3j/signcheck.py, 156 cases, 0 mismatches.
5. **Collecting.** (Q/N(q))·√N(q) = Q/√N(c₁c₂), and V̂ is radial, so V̂(√Qk/q) = (1+QN(k)/N(q))^{−B}. ∎

## 3. Approximate positivity

* For each k put r_k² := QN(k)/Z² and G_k(s) := (1 + r_k²e^{−s})^{−B}. In Lemma 2 the weight equals G_k(u_{c₁} + u_{c₂}).
* Define the formal form over **all** pairs (c₁,c₂) ∈ 𝒞²:

  D_k^{(j)} := Σ_{c₁,c₂} a^{(j)}_{c₁}(k)·conj a^{(j)}_{c₂}(k)·G_k(u_{c₁}+u_{c₂}),  with a^{(j)} built from W_j,

  and the k-independent sums P_j(k) := Σ_{c∈𝒞} μ(c)χ_c(k)W_j(N(c)/Z)/√N(c).

**Lemma 3.** There are η₀(B) > 0 and C such that, for η ≤ η₀ and all k and j,

  G_k(0)|P_j(k)|² ≤ 3D_k^{(j)} + C·Bη·G_k(0)·Σ_{n≥1} γ_n|P_{j+n}(k)|²,

where γ_n ≥ 0, Σ_n γ_n(n+1)^{A+3} ≤ 1, and nothing depends on k.

*Proof.*
1. **Derivative bounds.** log G_k(s) = −B log(1+y), with y = r_k²e^{−s}. The first three s-derivatives of log(1+y) are −y/(1+y),
   y/(1+y)² and −y(1−y)/(1+y)³, all of absolute value ≤ 1 for y > 0.
   * So log G_k(s) = log G_k(0) + g₁s + (g₂/2)s² + R₃(s), with |g₁| ≤ B, |g₂| ≤ B/4 and |R₃(s)| ≤ B|s|³/6, uniformly in k.
2. **Factorisation.** Write b_c := a_c·e^{g₁u_c+(g₂/2)u_c²}. Then D_k = G_k(0)Σ b_{c₁}b̄_{c₂}·(1 + H(u_{c₁},u_{c₂})), with
   H(x,y) = e^{g₂xy + R₃(x+y)} − 1 = Σ_{a+b≥2} h_{ab}x^ay^b, symmetric with real coefficients.
3. **Smallness of H.** By the majorant series, ω_a := Σ_b |h_{ab}|η^{a+b} satisfies Σ_a ω_a(a+1)^{A+3} ≤ C·Bη². Also ω₀ ≤ CBη³, because
   the pure y^b terms come only from R₃.
4. **Lower bound for D_k.** Put B̂_a := Σ_c b_c(u_c/η)^a. Then |Σ h_{ab}B_a·conj B_b| ≤ Σ_a ω_a|B̂_a|², using |xȳ| ≤ ½(|x|²+|y|²), so

   D_k ≥ G_k(0)[(1−ω₀)|B̂₀|² − Σ_{a≥1} ω_a|B̂_a|²].

5. **Removing the exponential factor.**
   * e^{g₁u+(g₂/2)u²} = Σ_n e_n(u/η)^n, with e₀ = 1 and Σ_{n≥1}|e_n|(n+1)^{A+3} ≤ 2Bη for η ≤ 1/(8B).
   * So B̂_a = Σ_n e_nP_{a+n}.
   * By Cauchy–Schwarz, |P₀|² ≤ 2|B̂₀|² + 4Bη·Σ_{n≥1}|e_n||P_n|² and |B̂_a|² ≤ 2Σ_n|e_n||P_{a+n}|².
6. **Combining.** Steps 4 and 5 give the claim for j = 0. Replacing W by W_j only relabels P_n as P_{j+n}. ∎

*Numerical check (l3j/lemma3test.py).*
* **Setup.** 20 000 random cases with η ∈ {0.05, 0.1, 0.2}, 5–60 coefficients and r² ∈ [10⁻³, 10³]. They include zero-mean coefficients
  and coefficients that oscillate in u.
* **Quantity.** The ratio (G(0)|P₀|² − 3D)/(Bη·G(0)·Σ_{n≥1}2^{−n}|P_n|²).
* **Result.** The ratio never exceeded 1.05, so C ≈ 1 suffices in practice.
* **Necessity of the correction term.** D < 0 occurs in 26% of the random zero-mean cases at η = 0.2.

**Corollary 3.** Let 𝒜 := sup_{j≥0} (j+1)^{−A−3}·Σ_k G_k(0)|P_j(k)|², which is finite for fixed Z. For η ≤ η₀(B) small enough that
C·Bη ≤ ½:

  𝒜 ≤ 6·sup_j (j+1)^{−A−3}·Σ_k D_k^{(j)}.

## 4. Upper bound for Σ_k D_k^{(j)}

**Lemma 4.** Under (G)^{tw}_{θ′}, for all j, Σ_k D_k^{(j)} ≪ (j+1)^{A+3}·(Z·Q^{θ′} + Z²/Q)·Z^ε.

*Proof.* The proof is by induction on Z, with the statement taken at once:
* for all row scales Q;
* for all index-twisted families, i.e. coefficients T(c)χ̄_c(h₀) on c ∈ 𝒞 (equivalently c coprime to h₀, with rows m ↦ mh₀), uniformly
  in h₀ up to N(h₀)^ε;
* for both classes κ.

Lemma 1 holds for this whole class, because index twists compose (mh₀ ↦ (mh₀)³ = m³h₀³). Split the pairs by gcd.

1. **Coprime pairs (c₁ ≠ c₂).** By Lemma 2, Q·Σ_k D_k^{cop} = S^{cop}, the T-side contribution of coprime pairs. Inclusion–exclusion
   over common divisors gives

   S^{cop} = Σ_h μ(h)S^{(h|)},  where S^{(h|)} := Σ_m V(m/√Q)|Σ_{c∈𝒞, h|c} T(c)W_j(N(c)/Z)χ̄_c(m)|² ≥ 0.

   * By twisted multiplicativity, T(ha) = T(h)T(a)(h/a)₂s(h,a), and (h/a)₂χ̄_a(m) = χ̄_a(mh²). The class restriction makes s(h,a)
     constant. So S^{(h|)} ≤ the gcd-twisted rows of Lemma 1, which are ≪ (j+1)^{A+2}(Z/N(h))Q^{1+θ′+ε}N(h)^ε.
   * Summing |μ(h)|S^{(h|)} over N(h) ≤ 2Z gives S^{cop} ≪ (j+1)^{A+2}ZQ^{1+θ′}Z^{2ε}.
2. **Pairs with gcd > 1, including c₁ = c₂.** Their formal contribution is Σ_{h>1} −μ(h)·F_h, where F_h is the formal form over pairs
   both divisible by h.
   * a_{ha}(k) = μ(h)χ_h(k)N(h)^{−½}·μ(a)χ_a(k)W_j(N(ha)/Z)N(a)^{−½}, and u_{ha} is still in [0,η]. So F_h(k) = N(h)^{−1}1_{(k,h)=1}·D′_{h,k},
     where D′ is the Step-3 form for the Möbius sums Σ_{(a,h)=1} μ(a)χ_a(k)… of length Z/N(h) with the same kernel G_k.
   * The restriction (a,h) = 1 is the index twist h₀ = h⁴, since 1_{(a,h)=1}χ_a(k) = χ_a(kh⁴). So D′ lies in the class covered by the
     induction, in the class κ − e(h) mod 2.
   * Writing G_k with Z′ = Z/N(h) and Q′ = Q/N(h)² gives r_k² = Q′N(k)/Z′². So D′ belongs to the length-Z′, row-scale-Q′ problem.
   * By the upper-bound version of Lemma 3 (|D′| ≤ G(0)[(1+ω₀)|B̂₀|² + Σ ω_a|B̂_a|²]) and the induction hypothesis at length Z′ < Z
     (for all j, row scale Q′):

     |Σ_k F_h(k)| ≪ (j+1)^{A+3}N(h)^{−1}·(Z′Q′^{θ′} + Z′²/Q′)·Z^ε ≤ (j+1)^{A+3}·(ZQ^{θ′}N(h)^{−2} + (Z²/Q)N(h)^{−1})·Z^ε.

   * When Q′ < 1 the T-side of the shorter problem is empty up to O(Z^{−10}), and the bound is just its diagonal.
   * Summing over N(h) ≤ 2Z costs log Z.
3. **Conclusion.** Σ_k D_k^{(j)} = Q^{−1}S^{cop} + Σ_{h>1} −μ(h)Σ_kF_h(k), which gives the claim. ∎

**Proposition 5.** Under (G)^{tw}_{θ′}, 𝒜 ≪ (Z·Q^{θ′} + Z²/Q)·Z^ε. *(Corollary 3 with Lemma 4.)*

## 5. Repeated frequencies

**Proposition 6.** Fix k₀. Put M := P₀(k₀) evaluated with c ∈ 𝒞_κ, i.e. M = Σ_{c∈𝒞_κ} μ(c)(k₀/c)₄W(N(c)/Z)/√N(c). Then

  |M|² ≪_{k₀} Z^{(3/4)(2−x)+ε} = Z^{(3/4)(1+2θ′)/(1+θ′)+ε}.

*Proof.*
1. **The frequencies.** Let J := (Z²/(QN(k₀)))^{1/4}, and let 𝒫 be the set of primes ϖ with N(ϖ) ∈ [J/2, J]; then |𝒫| ≫ J/log J.
   For k = k₀ϖ⁴ we have r_k² = QN(k₀)N(ϖ)⁴/Z² ≤ 1, so G_k(0) ≥ 2^{−B}.
2. **The repeated sum.** χ_c(k₀ϖ⁴) = χ_c(k₀)·1_{ϖ∤c}, so P₀(k₀ϖ⁴) = M − E_ϖ with |E_ϖ| ≤ Σ_{ϖ|c}N(c)^{−½} ≪ Z^{½}/J.
3. **Case |M| ≥ 2 max|E_ϖ|.** Positivity of every term in 𝒜 gives 2^{−B}(|𝒫|/4)|M|² ≤ 𝒜. So |M|² ≪ (ZQ^{θ′} + Z²/Q)·J^{−1}·Z^ε.
   * With Q = Z^x this is (Z^{1+xθ′} + Z^{2−x})·Z^{−(2−x)/4}.
   * At x = 1/(1+θ′) the two terms are equal, and the exponent is (2−x) − (2−x)/4.
4. **Otherwise** |M| ≪ Z^{½}/J = Z^{x/4}N(k₀)^{¼}. This is smaller, because x/4 ≤ ¼ < (3/8)(1+2θ′)/(1+θ′). ∎

## 6. Zeros

* **The class sums.** 1_{e(c)=κ} = ½(1 + (−1)^κλ(c)). So Proposition 6 for κ = 0 and κ = 1 bounds both Σ_c μ(c)(k₀/c)₄λ(c)^e W(N(c)/Z)
  (e = 0, 1), over squarefree primary c, by Z^{½}·|M| ≪ Z^{σ(θ′)+ε}.
* **From annuli to dyadic ranges.**
  * The bound is uniform over the narrow-annulus bumps W, all of which have support of length η in log N(c).
  * A smooth φ supported in [1/2, 2] is a finite sum of translates W(·/e^{iη}) with smooth coefficients.
  * So Σ_c μ(c)ξ(c)φ(N(c)/X) ≪ X^{σ+ε} for every such φ.
* **Mellin.** φ̃(s)·L(s,ξ)^{−1} = Σ_c μ(c)ξ(c)N(c)^{−s}φ̃(s) is the Mellin transform of a dyadic sum over X. It is therefore holomorphic
  in Re s > σ.
  * Here L(s,ξ) is the Euler product over odd primes. The omitted 2-factor has no zeros.
* **Zeros.** If L(ρ, ξ) = 0 with Re ρ > σ, choose φ with φ̃(ρ) ≠ 0, which is possible since φ̃ is entire and not identically zero.
  Then φ̃/L has a pole at ρ, a contradiction. ∎

## What is used and from where

* PROBLEM.md: the Gauss-sum identities.
* L3_REPORT §10 step 3: T·g = μχ√N, checked numerically.
* Quartic reciprocity and the supplementary laws: standard.
* Kubota (with DDHL (4.28)): ψ_α entire, with polynomial growth.
* signcheck.py: the sign in Lemma 2.

* The proof uses no unproved input beyond the hypothesis (G)^{tw}_{θ′}.
* Passing from the Dirichlet-series form (L³)^{tw} to (G)^{tw} needs the local 2×2 inversion of Lemma 1 (step 4) on Re s = ½. That
  inversion is not done here.
* Places where an independent reader should check details: the majorant-series constants in Lemma 3, and the induction bookkeeping in
  Lemma 4.
