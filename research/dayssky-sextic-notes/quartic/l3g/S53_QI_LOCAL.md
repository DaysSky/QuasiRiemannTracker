# S53. Task #26, local theory: the Chinta–Gunnells G₂ series over Q(i) with n = 4, where μ₈ ⊄ Q(i)

DaysSky with Claude assistance, 16 Oct 2026. Labels: PROVED / NUMERICAL / STANDARD / OPEN.

## 0. Summary
* **What breaks.** Over Q(i), with n = 4 and S = {(1+i), ∞}, the self-symbol is (ϖ,ϖ)_S = (−1/ϖ)₄ =: s_ϖ, which is −1 when N ϖ ≡ 5 mod 8.
  Chinta–Gunnells (CG, JAMS 2010) assume μ_{2n} ⊂ F, so s_ϖ = 1 for them. Their proof uses the hypothesis in two places:
  * the local Weyl action, where γ(i)γ(−i) = 1/q;
  * the slice expansion (their Lemma 5.9), where self-symbols (ϖ^a, ϖ^b)_S are dropped.
* **What replaces it.** Keep the self-symbols and the raw Gauss sums. The resulting local action is CG's Definition 3.1 with Q multiplied
  by a sign:

      Q′_{β,ℓ,i}(x) = s^{μ_{ℓ,i}(β)·(B_i(β)+1)} · Q_{β,ℓ,i}(x),   with B_i(β) = b_i + Σ_{j>i} c(j,i) b_j,

  with raw γ(i) = g(1,ϖ;εⁱ)/q, so that γ(i)γ(−i) = sⁱ/q. Call this the **slice-signed action**.
* **Results.**
  * The slice-signed action is a W-action. Its averaged p-parts N(x;ℓ) are polynomials.
  * Those p-parts satisfy the local functional equation that the TRUE one-variable slices (Kubota series with Hecke characters over
    Q(i)) require.
  * Both statements are PROVED: exact computer algebra in Q(√q, γ(1), x₁, x₂), plus a two-line sign argument for general ℓ.
* **Correction to S51 §19.** S51 §19's Brylinski–Deligne-signed action γ(i)s^{C(i,2)} is also a W-action, but it is the wrong one. It
  gives the d-axis coefficient −g(1,ϖ;ε³) instead of the true Gauss sum, and it fails the slice functional equation in 343 of 554
  classes.
* **The two p-parts differ only by monomial signs** (same support, same |coefficients|). So the exponent analysis of
  S52_CORRECTION_NAIVE_VS_CG.md (E′ convergence, slice comparison) is unchanged.

## 1. The slice expansion with self-symbols — PROVED (Hilbert-symbol algebra; rank-1 check in §3)
Setting: CG §4–5 notation, slice in node i, one prime ϖ ∉ S, a_j = ϖ^{β_j} for j ≠ i, m_i = ϖ^{λ_i}m_i′, l = λ_i + Σ_{j≠i}|c(j,i)|β_j,
t = ‖α_i‖² (1 for the short c-node, 3 for the long d-node; both odd), x = q^{−s_i}.

* **(a) Reduction to a ϖ-coprime Kubota series.** Twisted multiplicativity (CG (4.3)–(4.5)) plus reciprocity (c/b) = (b,c)_S(b/c)
  give E(s) = Σ_k N^{(k)}(x)·D_ϖ(s, u₀ϖ^{K(k)}; θ_k), with K(k) = (l − 2k) mod 4 and
  θ_k(c) = Ψ(…, ϖ^k c, …)(c,ϖ^k)^{−t}Π_{j>i}(c,a_j)^{−2⟨α_i,α_j⟩}.
  * The Hilbert symbols from reciprocity cancel those in θ_k by antisymmetry, (x,y)_S(y,x)_S = 1. So there is no self-symbol here.
  * For node i = 2, the cross term (c₁/c′₂) needs no reciprocity. CG (5.7) twists Ψ_i only by j > i, and the two facts match.
* **(b) CG Lemma 5.4 holds verbatim.** Inverting CG Lemma 5.3 produces two signs:
  * (Ψ̂_η)^_{η₂} = (ϖ^{i+1}, ϖ^{3−i})_S^t·Ψ = s^{(i+1)(3−i)}Ψ;
  * g(u₀ϖ^i, ϖ^{i+1})·g(u₀ϖ^{2−i}, ϖ^{3−i}) = q³s^{i+1}.

  Their product is s^{(i+1)(4−i)} = 1, because (i+1) + (4−i) is odd. So the denominator stays 1 − q³x⁴.
* **(c) Re-indexing (CG Lemma 5.9) gains a self-symbol.** Let k′ = (l+1−k) mod 4, K′ = K(k′) and η = ϖ^{K′+1}. Then

      θ̂_{k′,η} / θ_k = (η, ϖ^{k′})_S^{−t} Π_{j>i}(η, a_j)_S^{−2⟨α_i,α_j⟩} = s^{(K′+1)(k′ + Σ_{j>i} c(j,i)β_j)} =: S(k′).

* **(d) The true slice factor** is therefore

      f′_k(x) = [N^{(k)}(x) − 1_{K′≠3}·S(k′)·q^{K′}g(u₀,ϖ;ε^{t(K′+1)})·x^{K′+1}·N^{(k′)}(x)] / (1 − q³x⁴).

  The rank-1 Kubota functional equation (DDHL (4.28), §4) gives the slice functional equation if and only if (Theorem 4.1′)
  f′_k(x) = (qx)^{l−K(k)}·f′_k(1/(q²x)).
  * The exponent l − K(k) is CG's l+1−(k′−k)_m. The printed (k−k′) in CG Theorem 4.1 is the same up to the convention for k′ (checked
    at s = +1 below).

## 2. The local theorems — PROVED (l3g/relflint.py, thm41sym.py: exact rational functions over FLINT, √q = r and γ(1) = a indeterminate)
* **Theorem A (W-action).**
  * Claim: for every ℓ ∈ Z_{≥0}², the slice-signed action satisfies σ₁² = σ₂² = (σ₁σ₂)⁶ = 1 on all 16 classes.
  * Checked identically in Q(r,a,x₁,x₂) for ℓ ∈ {0..7}², with γ(2) = ±1/r.
  * General ℓ: ℓ ↦ ℓ + 4e_k multiplies P and Q at node k by (qx_k)⁴ and leaves the signs and γ-indices unchanged (μ shifts by 4). With
    M = x^{−4ω_k} (ω_k ∈ Λ for G₂, and −4ω_k ∈ Λ′), f|_{ℓ+4e_k}w = M^{−1}·((Mf)|_ℓ w). So the actions are conjugate, and the
    relations follow from ℓ mod 4. ∎
  * Controls: raw Gauss sums fail 28/48 at s = −1; s = +1 reproduces CG Thm 3.2.
* **Theorem B (polynomial p-parts and Theorem 4.1′).**
  * Claim: N(x;ℓ) = h·D, with h = Σ_w (1|_ℓ w)/Δ(wx), is a polynomial in x for every ℓ, and every slice class satisfies Theorem 4.1′
    with the sign S of §1(c).
  * **Polynomiality.** CG Thm 3.5's proof only uses that |_ℓ is a W-action of the shape (3.14). Its (3.30) allows any constant γ, so
    the sign is irrelevant.
  * **Theorem 4.1′ for all ℓ.** By CG Thm 3.6's proof (invariance of h plus two pointwise identities), it suffices that
    * (I1) P_β − c_f(qx)^{m−δ}Q′_{σ•β} = (qx)^{l+1−δ} and
    * (I2) (qx)^{2(m−δ)}P_{σ•β} − Q′_β(qx)^{m−δ}/c_f = (qx)^{l+1−δ},

    where c_f = S(k′)γ(−tδ) and δ = μ(β) mod 4.
  * **Sign bookkeeping.**
    * In (I1), the sign is S·s^{δ(B(σ•β)+1)}·s^{tδ}, from γ(−tδ)γ(tδ) = s^{tδ}/q. Since S = s^{δB(σ•β)}, this equals s^{δ(1+t)} = 1
      (t is odd).
    * In (I2), the sign is s^{δ(B(β)+1)+δB(σ•β)}. With B(σ•β) ≡ B(β) + δ, this equals s^{δ(1+δ)} = 1.
    * So both identities reduce to CG's own, with γγ = 1/q. ∎
  * **Direct exact check.** For ℓ ∈ {0..7}²: N is polynomial, and Theorem 4.1′ holds in 4191 of 4191 slice classes (469 of 469 for ℓ ≤ 3,
    with γ(2) = −1/r too).
  * **Controls.** s = +1 reproduces CG (all pass). The BD-signed action of S51 §19 gives polynomial N but fails 343/554 classes. Two other
    sign candidates (s^{μB}, and no sign) fail.
* **Comparison.** For ℓ = (0,0) the slice-signed and BD-signed p-parts have the same 14 monomials. The coefficient ratios are ±1:
  (0,1), (1,2), (4,1), (5,2) and (9,6) flip, and the rest agree. The same holds at ℓ = (1,0), (0,1), (1,1) (l3g/compare.py).
  * So the q-exponents in S52_CORRECTION (13 non-constant monomials, E′ thresholds) are unchanged.
  * H_{x₂¹} = (raw Gauss sum)·q^0 still has no x₁-dependence, so the c-slices at squarefree d still equal the naive series.

## 3. Rank-1 check of the self-symbol — PROVED (l3g/rank1sym.py, exact in (r,a))
* With the true Gauss-sum local factor at ϖ^l ∥ index, N(x) = Σ_k g(ϖ^l, ϖ^k; ε^t)x^k, the factor f′ of §1(d) satisfies
  f′(x) = (qx)^{l−K}f′(1/(q²x)) for l = 0..12 and t = 1, 3 at s = ±1.
* Without S it fails for every even l at s = −1.
* Kubota's theorem says the true series satisfies this, so the check confirms the derivation of S, including its sign convention.

## 4. Global step (CG §5–6 over Q(i)) — PROVED from DDHL (4.28) and Remark 4.2 (= [Dia04] (2.35)–(2.48), [Suz83])
* **Inputs used by CG Thm 5.8 / 6.1:**
  * Lemma 5.3/5.4: §1(b), verbatim.
  * Lemma 5.9: §1(c), with the extra S.
  * Thm 4.1: Theorem B.
  * Prop 5.6/5.7 (Ψ-bookkeeping): below.
  * Bochner/Phragmén–Lindelöf: STANDARD.
* **(R1) Rank-1 FE in operator form — from the literature.**
  * Statement. For odd b ∈ Z[i], φ ∈ M_t(Ω) and ∞-type ℓ:

        D*(s,b;φ;ℓ) = N(b)^{1−s}α(b)^ℓ·D*(2−s,b;T_bφ;−ℓ).

    Here D(s,b;φ;ℓ) = Σ_c g(b,c;ε^t)φ(c)α(c)^ℓN(c)^{−s}, and D* includes the index-free factors G_∞ and ζ(4s−3).
    T_b ∈ End M_t(Ω) ⊗ C(2^{−s}) depends on b only through its class in Z₂[i]^×/(Z₂[i]^×)⁴.
  * Source. DDHL (4.28) holds for every ν with A_{ji} independent of ν. For (γ_i, λ) = 1, ψ_{i1} is ψ_β(s,ν,ℓ) up to a constant
    (4.15). Every ψ_{j1} is a combination of ψ_β(s,(−1)^aλ^{2b}ν,ℓ) whose coefficients depend on ν only through ν′ mod 4,
    (i^aλ^b/ν′)₄ and e(·ν/λ⁴) ((4.16)–(4.19)); all of these are invariant under ν ↦ νz⁴.
  * Re-read check (08:10Z). Each 2-adic factor in DDHL (4.15)–(4.19) is invariant under ν ↦ νz⁴ for odd z. For z primary,
    z⁴ ≡ 1 mod λ⁷, which gives the following.
    * ě(−γ′αν/λ⁴) changes by e(Tr(λZ[i])) = 1.
    * (i^aλ^b/ν′z⁴)₄ = (i^aλ^b/ν′)₄.
    * ν′ mod 4 is unchanged, and S(κ,ν) depends only on (j,k) = (unit, ord_λ).
    * ψ₁ and ψ_{1+λ³} are themselves cusp components, (4.12)–(4.13).
    * The extra indices λ^{2b}ν are S-unit twists, with factors 2^{−bs}.
    So (R1) holds as stated.
  * Spanning. For S-units w, D(s,bw;φ) = D(s,b;φ·(w/·)^{−t}), and (w/c) = (c,w)_S^{±1}. The characters (·,i)_S and (·,λ)_S
    generate the dual of F_S^×/Ω, because the pairing with Ω/F_S^{×4} = ⟨i,λ⟩ ≅ (Z/4)² is perfect. So the ψ_β(s,wν) span every
    D(s,ν;φ), and T_ν is well defined.
* **(R2) Injectivity — PROVED.** If D(s,b;ψ;ℓ) ≡ 0 for some ψ ∈ M_t(Ω) ⊗ C(2^{−s}), then ψ = 0.
  * First separate the powers of 2^{−s}; the c run over odd ideals.
  * The coefficient at a split prime norm p ∤ 2N(b) is g(b,π)ψ(π)α(π)^ℓ + g(b,π̄)ψ(π̄)α(π̄)^ℓ, with π, π̄ the two ideals of norm p.
  * If ψ(C) ≠ 0, then ψ(C̄) ≠ 0 and the ratio R(π) = g(b,π)α(π)^ℓ/(g(b,π̄)α(π̄)^ℓ) is constant on the primes π ∈ C.
  * Since g(b,π;ε^t)⁴ = g(1,π;ε^t)⁴ = (π³π̄)^{±1} (DDHL (3.10)), R⁴ = α(π)^{4(1+2ℓ)} for t = 1, or α(π)^{4(2ℓ−1)} for t = 3. These
    are nontrivial Grössencharacters (ℓ ∈ Z).
  * By Hecke, their values on degree-1 primes in a ray class are equidistributed. Contradiction. ∎
* **(R3) The slice identity at ϖ⁴ — PROVED (every coefficient checked by hand).** For ϖ ∤ u₀ and x = q^{−s}:

      D(s,u₀ϖ⁴;φ) = (1+q⁴x⁴)·D(s,u₀;φ) − g(u₀,ϖ;ε^t)·x·D(s,u₀ϖ²;θ₁φ),   θ₁φ(c) = φ(ϖc)(c,ϖ)_S^{−t}.

  * Checked coefficient by coefficient at c′, ϖc′, …, ϖ⁵c′, using CG Prop 2.1 (1)–(3), which hold for μ₄ ⊂ F.
  * At ϖ⁴c′ the left side is (q⁴−q³)g(u₀,c′)φ(c′). The right side gives q⁴ − q³·s², where one s is the self-symbol (ϖ,ϖ)_S^{3t} in θ₁ and
    the other comes from g(1,ϖ;ε^t)g(1,ϖ;ε^{−t}) = s^t q.
  * This independently confirms §1(c).
  * The local factors 1+q⁴x⁴ and −g(u₀,ϖ)x are §1(d)'s f′₀ and f′₁ for l = 4, and they satisfy Theorem 4.1′.
* **(E1) Equivariance — PROVED.**
  * **Step 1: apply the FE both ways.** Apply (R1) to the left side of (R3), then re-expand by (R3) at 2−s with −ℓ. Apply (R1) to each
    term on the right side and use Theorem 4.1′ (the α(ϖ)-powers cancel, as in the twist bookkeeping). This gives, with s′ = 2−s and
    x′ = q^{−s′},

        f₀(x′)·D(s′,u₀;Δ₀) + f₁(x′)·D(s′,u₀ϖ²;Δ₁) ≡ 0,
        Δ₀ = T_{u₀ϖ⁴}φ − T_{u₀}φ,   Δ₁ = θ₁T_{u₀ϖ⁴}φ − T_{u₀ϖ²}θ₁φ.

  * **Step 2: Δ₀ = 0.** At split primes π ∉ {ϖ, ϖ̄} only the constant term f₀(0) = 1 contributes. So (R2) gives Δ₀ = 0, i.e.
    T_{u₀ϖ⁴} = T_{u₀}.
  * **Step 3: Δ₁ = 0.** Then f₁ = −g(u₀,ϖ)x′ ≠ 0 forces D(s′,u₀ϖ²;Δ₁) ≡ 0, and (R2) gives Δ₁ = 0.
  * **Base case.** Together: T_{u₀ϖ²}θ₁ = θ₁T_{u₀}. This is (E1) for (b,y) = (u₀,ϖ), with ϖ ∤ b.
  * **All b.** For b divisible by ϖ, choose ϖ′ in the same F_S-class with ϖ′ ∤ b. Then θ_{ϖ′} = θ_ϖ, and T_{bϖ′^{−2}} = T_{bϖ^{−2}} by
    (R1)'s class dependence.
  * **All y.** For y = ϖ^k, compose: θ_{y₁y₂} = (y₂,y₁)^t·θ_{y₂}θ_{y₁}, and the constants agree on both sides.
* **(E2) Ω-covariance — PROVED.**
  * For S-units u, D(s,bu;φ·(·,u)^t) = D(s,b;φ) identically (reciprocity, with (c/u) = 1). So T_{bu}(φ·(·,u)^t) = (T_bφ)·(·,u)^t.
  * For u ∈ F_S^{×4}, (·,u)_S = 1, and T_{bu} = T_b by (R1).
* **Prop 5.6/5.7.** Define σ_iΨ by Φ(a) := [T_{m_ib_i(a)}Ψ_i^{(a)}](1).
  * (E1) gives Φ_i^{(a)} = T_{m_ib_i(a)}Ψ_i^{(a)} for every a, via Ψ_i^{(…a_ic…)} = (c,a_i)^tΠ_{j>i}(c,a_j)^{2⟨α_i,α_j⟩}(Ψ_i^{(a)})_c.
  * (E2) gives Φ ∈ M(Ω,Φ).
  * So CG Thm 5.8 holds over Q(i): E*(s_i,a;Ψ) = |Am_i|^{1−s_i}α(Am_i)^ℓ·E*(2−s_i,a;σ_iΨ). Then CG Thm 6.1 (FE and continuation by
    Bochner) holds as written.
* **The α-twist is consistent — PROVED (bookkeeping).** In the slice, b_k = u₀ϖ^{K(k)} gives α(b_k)^ℓ = α(u₀)^ℓα(ϖ)^{ℓK(k)}.
  Theorem 4.1′ in the twisted variable α(ϖ)^ℓx contributes α(ϖ)^{ℓ(l−K(k))}. The product is α(m_iA)^ℓ, independent of k: the twist rule
  of S51 §17.
* **Status of #26.** PROVED, modulo:
  * the published rank-1 inputs DDHL (4.28) and Remark 4.2 (from [Dia04], [Suz83]), read as in (R1);
  * Hecke's equidistribution theorem;
  * the STANDARD Bochner/Phragmén–Lindelöf step.

  No Brubaker–Bump μ₈ hypothesis is used.
* **Literature.**
  * Chinta–Offen (metaplectic Casselman–Shalika for GL_r) work with μ_n ⊂ F only. Their scattering coefficients carry a sign
    (−1,ϖ)^{(quadratic in λ)} of the same shape as S (their (7.9)–(7.10)); they assume (−1,ϖ) = 1 only when matching CG.
  * Chen (arXiv 2411.13143) needs 2n | q−1 only for (ϖ,ϖ) = 1 (his Prop 2.2.1(iii)) and for the Gauss-sum relation g_k g_{−k} = q
    (Prop 2.3.2(v)). Theorems A and B are the G₂, n = 4 replacements.

## 5. Consequences
* The G₂ series over Q(i) that S52 needs is the CG series with slice-signed p-parts. Its one-variable slices in every chamber are exact
  Kubota series (Theorem B plus §4). So S52 Lemma 4(f) ("other chambers verbatim") holds for it by the functional equation
  Z*(s;Ψ) = |m_i|^{1−s_i}Z*(σ_i s; σ_iΨ).
* The transfer back to the actual (L³) series is unchanged: naive = Z × E′ × normalisers^{−1} at (−1/12, ¼) and (x, −3x), with
  −1/12 < x < 0, where E′ converges absolutely. The exponents are identical by §2's comparison.
* **Status of (L³) ≪ Q^{3/2+ε}.** #26 is now PROVED modulo the published rank-1 FE (DDHL (4.28) and Remark 4.2) and STANDARD
  Bochner/Phragmén–Lindelöf (§4). With S52_LIMIT_PROOF and S52_FULLSUM (which use the framework (F1)–(F8)), this gives the full (L³) sum
  ≪ Q^{3/2+ε}. Its remaining items are S52_PROOF Lemma 4's non-coprime and cusp constants, labelled there as argued but not tracked.
* Unconditional state unchanged: Q^{13/8} (S49). (L³) at Q^{1+ε}: OPEN (BELOW_3_2.md).

## 6. S52 Lemma 4 for the CG series at the points the proof uses (l3g/lemma4_exp.py, lemma4_all.py)
* **Where the proof evaluates.** S52_LIMIT_PROOF only needs block bounds at the chamber points (x′,y′) = (⅓+δ′, 0) of chambers s_c and
  dcdc (Lemmas A–C). The transfer back to the actual series happens at (−1/12, ¼) and on (x,−3x), with −1/12 < x < 0, where E′
  converges (S52_CORRECTION).
* **Bochner needs no connectivity here.**
  * The normalised CG function becomes entire once multiplied by the linear forms of its finitely many polar lines (W-translates of the
    theta lines, CG Thm 6.1). Those lines carry S52's polar terms (S52 §16).
  * That product has finite order (STANDARD: Phragmén–Lindelöf with the functional equations).
  * On the tube over a convex set, sup_t log|F(σ+it)| is then a convex function of σ (three lines), up to polynomial factors in |Im|.
  * So the bound at a convex combination of evaluation points is the convex combination of their bounds, whatever happens between them.
* **Corrections at generic primes (ℓ = 0) — NUMERICAL, exact-arithmetic exponents.** The non-first-order monomials of the
  slice-signed p-part have worst normalised exponent e − k₁(1+x′) − k₂(1+y′):
  * −1.333 at (⅓, 0), from (1,2);
  * −1.5 at (½, 0).

  Both are < −1, so the correction Euler product converges absolutely. Lemma 4(b)–(d) then go through for the CG series as for the naive
  one.
* **Primes dividing the twist (h₀e², with h₀ fourth-power-free).**
  * For every s = ±1 and ℓ ∈ {0..7}², all p-part coefficients satisfy e(k) ≤ k₁ + k₂ (the maximum of e − k₁ − k₂ is 0, attained up to a
    constant ≤ 4).
  * So at every chamber point with x′, y′ ≥ 0 each local factor is O_ℓ(1), and the product over p | h₀ is ≪ N(h₀)^ε.
  * Primes dividing only e cost at most N(e)^{O(1)} = Y^{O(ε)}.
* **Status.**
  * Lemma 4(a)–(e) for the CG series at the points used: verified, with exponents read from exact arithmetic at √q = 10007 and 10009.
    The exponents are half-integers, so rounding is exact.
  * Lemma 4(f) (other chambers): CG Thm 6.1 over Q(i) (§4).
  * So the S52 chain for (L³) ≪ Q^{3/2+ε} has no open item beyond published inputs:
    * DDHL (4.28) and Remark 4.2;
    * Gao–Zhao;
    * BGL;
    * Hecke;
    * Bochner.
  * It does still need an independent end-to-end audit before any upstream PR.
