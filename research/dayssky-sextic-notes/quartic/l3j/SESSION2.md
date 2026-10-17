# (L³) — second 5-hour block, 16 Oct 2026, from 08:20Z (DaysSky with Claude assistance)

Target: Σ_{N(m)≤Q} |ψ_α(m³,½+it)|² ≪ Q^{1+ε}. Notation as in l3f/BELOW_3_2.md:
* Z = c-length, Q = Z^{2/3}, K = Z²/Q = Z^{4/3};
* T-side: S(Q,Z) = Σ_m |Σ_c T(c)χ̄_c(m)W|²;
* μ-side: P(k) = Σ_c μ(c)ξ_k(c)W(N(c)/Z), with S − c₀QZ = (Q/Z)Σ_{N(k)≲K}(|P(k)|² − C_kZ).

Labels: PROVED, EXPONENTS (bookkeeping only), HEURISTIC, OPEN. Every entry ends with the exact point where the route stops.

## 1. Heath-Brown decomposition of μ on the μ-side (08:25–08:45)

* **Setup.** Write μ on N(c) ≤ 2Z with Heath-Brown's identity. The pieces are products of short Möbius variables a_i (N(a_i) ≤ Z^{1/J})
  and smooth variables b_l. Since ξ_k is a character, each piece factors as P_conf(k) = Π A_i(k)·Π B_l(k), after Mellin-separating the
  joint cutoff.
* **(a) Type I, small a: no off-diagonal (EXPONENTS).**
  * Take P_I(k) = Σ_{a~A}α_a ξ_k(a)Σ_b ξ_k(b)w(N(b)/B), with AB = Z.
  * Poisson in b (conductor N(k) ≈ K) gives P_I(k) ≈ (B/√K)·ε_k·Σ_{a, N(h)≲K/B} α_a ξ_k(a) ψ̄_k(h), plus the h = 0 term when ψ_k
    is trivial.
  * In |P_I|² the root number ε_k drops out.
  * An off-diagonal tuple (a₁,h₁,a₂,h₂) gives a k-character of modulus ≤ A⁴K²/Z². This is below K^{1−ε} iff A ≤ Z^{1/6−ε}, and
    then the smooth k-sum is negligible.
  * So for A ≤ Z^{1/6−ε} the variance of P_I is its diagonal-tuple part, up to O(Z^{−100}).
* **(b) The diagonal-tuple part is not the c-diagonal.**
  * For k = k₀j⁴ the b-character is (essentially) ξ_{k₀}. When k₀ is a unit, the Poisson zero mode h = 0 survives:
    |P_I(j⁴)|² ≈ B²|Σ_{(a,j)=1}α_a|².
  * There are #j ≈ K^{1/4} = Z^{1/3} such frequencies. Each Type I piece therefore carries a main term of size up to Z^{1/3}·Z², above the
    allowed Z².
  * Summed over all Heath-Brown pieces, these main terms reassemble into Σ_{N(j)≤K^{1/4}} |Σ_{(c,j)=1} μ(c)W|², the 4th-power
    frequencies of the full P (§2). They are small only if Σ_{c~Z}μ(c)W ≪ Z^{5/6}.
* **(c) Cross terms carry quartic theta coefficients at non-square indices.**
  * Between a piece with L₁ Poisson'd smooth variables and one with L₂ ≠ L₂′, the factor ε_k^{L₁−L₂} survives.
  * ε_k = g̃(ξ_k) is a quartic Gauss sum as a function of k, and it is not α-twisted. So Σ_k ε_k·(characters of modulus n) is a long
    Kubota sum in k, of length K. Its main term is τ₄(n)K^{3/4}, a quartic theta coefficient at the general index n.
  * Those coefficients are unknown off squares (Suzuki, Eckhardt–Patterson; numerically, at primes they vanish iff N(π) ≡ 5 mod 8,
    L3_ATTEMPT S29).
  * **Stops here.** The cross terms need explicit quartic theta coefficients at non-square indices. In the sextic (OpenAI) case the
    analogous coefficients are explicit (Patterson), which is exactly what the quartic family lacks.
* **(d) Type I₂ and balanced Type II stay as hard as (L³).**
  * Example: c = b₁b₂ with both b_i smooth and ≈ Z^{1/2} < K^{1/2}. Then no b_l can be Poisson'd usefully, the k-modulus of the
    off-diagonal is Z² against length K, and Poisson in k returns to the T-side.
  * Heath-Brown combinatorics cannot avoid this configuration.

## 2. 4th-power frequencies, zero-free regions, and a correction to NOTES §8G (08:45–09:20)

* **Exact accounting (PROVED algebra).** In the μ-side identity, the frequencies k = k₀j⁴ (N(j) ≤ (K/N(k₀))^{1/4}) all carry the same
  Möbius sum Σ_{(c,j)=1}μ(c)ξ_{k₀}(c)W.
  * For k₀ = 1 this is a ζ_{Q(i)}-type sum, twisted by norm-class characters mod 8 (Dirichlet L-functions mod 8 after factoring).
  * For other small k₀ it is a quartic Hecke L-function over Q(i). OpenAI's quasi-RH does not cover these (PROBLEM.md).
* **Naive requirement (HEURISTIC).** If nothing compensates, (G) at row-exponent M forces β(ξ_{k₀}) ≤ (6+M)/8 for the real parts β of
  zeros:
  * M = 1 gives 7/8;
  * M = 2/3 gives 5/6;
  * at M = 2/3 with loss Q^{θ′}, it gives 5/6 + θ′/3.
  So on this count, any bound below Q^{3/2} would already need a zero-free strip for quartic Hecke L-functions over Q(i).
* **The naive count is not reliable.**
  * At M = ½ it demands 13/16. But (G) at M = ½ is equivalent, by Kubota reflection (L3_REPORT §4), to (G) at M = 1 at scale Z^{1/2},
    whose count demands only 7/8.
  * The same statement cannot have two thresholds. So for zeros with 13/16 < β ≤ 7/8, the other frequencies must compensate at
    M = ½, and nothing rules out the same at M = 2/3.
* **Positivity fails, so the (L³) ⟹ PO step is not established.**
  * Turning (L³) into a positive k-sum (Theorem Z's hypothesis PO) needs the dual weight V̂(√Q·k/(c₁c₂)) to be positive
    semi-definite in (c₁,c₂) for each k.
  * Write the weight as G(log N(c₁) + log N(c₂)). It is a Hankel kernel, PSD only if G is exponentially convex, hence log-convex.
  * A log-convex function cannot be flat and then decaying. But a flat weight up to N(k) ≈ K is exactly what the K^{1/4} repetition
    of Theorem Z needs.
  * Power weights V̂ = r^{−2λ} are PSD and summable for λ > 1. But they give the frequencies the weights N(k)^{−λ}, so the j⁴
    repetition shrinks to O(1) effective copies, and the conclusion becomes trivial.
* **Correction.**
  * NOTES §8G's "Corollary: (L³) ⟹ 7/8 zero-free region" and L3_REPORT §5's "(L³) ⟹ consequences" were outlines whose positivity step
    fails for wide c-ranges.
  * Theorem Z itself (PO ⟹ zero-free) is unaffected.
* **Resolved at 09:00Z, l3j/THEOREM_ZF.md.**
  * Restrict c to a narrow annulus N(c) ∈ [Z, Ze^η]. On that annulus the weight G_k(u₁+u₂) is a product times 1 + O(η²).
  * The O(η²) error is a combination of |μ-side sums|² with other weights, so it is absorbed by taking a max over a finite family of
    weights.
  * Optimising the row scale Q = Z^{1/(1+θ′)} gives: (L³) with loss Q^{θ′}, θ′ < ½, implies a zero-free half-plane
    Re s > (7+10θ′)/(8+8θ′) for every quartic Hecke L(s,(k₀/·)₄λ^e) over Q(i).
  * In particular (L³) gives 7/8, and any bound Q^{3/2−δ} gives a strip of width about δ/6. Two parts are STANDARD and not written out:
    the gcd > 1 pairs and the non-squarefree c.
* **Consequence.** Below Q^{3/2}, every route must prove a quasi-RH for that family. This is a rigorous explanation of the 3/2 wall.
  The naive count at M = 2/3 (5/6) is still unjustified, for the reflection reason above.

## 3. Amplification over all characters mod m (09:20)

* **Construction.** Embed ψ_m = (·/m)₄ in the full character group mod m, with an amplifier w(χ) = |Σ_{a∈𝒜} χψ̄_m(a)|²/|𝒜|.
* **The full family is balanced (EXPONENTS).** Q² characters against AFE length about Q². But positivity alone loses the factor Q.
* **The amplified sum** is
  Σ_m φ(m) Σ_{a,a′} ψ_m(a′ā) Σ_{ca′≡c′a (m)} T(c)T̄(c′).
  * Off the diagonal ca′ = c′a it becomes shifted convolutions Σ_{ca′−c′a=h} T(c)T̄(c′)·F(h) of quartic Gauss sums in the modulus
    variable.
* **Stops here.** There is no spectral or Kubota control of correlations of Gauss sums at shifted moduli; this is a Chowla-type
  problem for Gauss sums.

## 4. t-average and Mellin positivity (09:25)

* Averaging |ψ_α(m³,½+it)|² against ĥ(t) ≥ 0 makes the c-kernel h(log N(c₁)/N(c₂)) PSD.
* In the triple series Z_T(s,s̄,w), the w-functional equation gives a Hermitian μ-side series Σ_k N(k)^{w−1}|Σ_c μξ_k N(c)^{−s−w+½}|²
  for real w. Its coefficients are nonnegative.
* **Stops here.** The functional equation sends the region where (L³) gives absolute convergence (Re w > 1) to b = 1 − w < 0, where
  the μ-series diverges. Landau positivity acts only at an abscissa of convergence, so it gives nothing there.

## 5. GRH with exp(short prime sums) for 1/L (09:30) — EXPONENTS

* Under GRH, 1/L(½+σ₀+it, ξ_k) = exp(−P_x) with x a power of log, and the exponential can be expanded into a Dirichlet polynomial of
  length K^{1/2}. This gives Σ_k|1/L|² = diagonal·(1 + O(K^{−σ₀/4})).
* Moving P(k) to the line Re s = ½ + σ₀ costs Z^{2σ₀}. The net error is Z^{7/3 + 5σ₀/3}, against the allowed Z².
* So this route fails for every σ₀ > 0, which confirms that GRH by these means gives Q^{3/2} and no more.

## 6. Hybrid quartic × quadratic sparse large sieve for Type II (09:35) — EXPONENTS

* Type II on the T-side (with Cauchy over the short variable) needs Σ_{N(m)≤Q} Σ_{a~A} |Σ_{b~B} β_b (m³a²/b)₄|² ≪ (QA + B)^{1+ε}‖β‖².
  This would give QZA + Z²/A, better than the floor for 1 < A < Z^{1/3} and best (QZ^{7/6}) at A = Z^{1/6}.
* By duality, the off-diagonal is a product of an incomplete quartic sum in m (length Q, modulus about B²) and an incomplete
  quadratic sum in a (length A, modulus about B²). Here B² > Q and B² > A, so neither is complete.
* Heath-Brown's quadratic sieve in a ignores the m-orthogonality and gives QZ(A + B). Gao–Zhao or S49 in m gives back the floor.
* Even granted, balanced Type II (A ≈ B ≈ Z^{1/2}) gives QZ^{3/2}, and T has no μ-factor to decompose into favourable pieces
  (§1(d)). **Stops here.**

## 7. Numerical confirmation of the THEOREM_ZF mechanism (09:05–09:20) — NUMERICAL
* **The test.** Take the control β = μT, whose Möbius side is μ², the extreme "pole" case. At Z = Q = 500, 2000, 4000 (all m ∈ Z[i],
  Gaussian row weight):
  * observed S/diag − 1: 0.058, 0.159, 0.194;
  * exact prediction from the frequencies k = u·j⁴: 0.099, 0.173, 0.222.
  It grows at the predicted rate, with no compensation at M = 1. For T itself the excess stays ≤ 0.03.
* **Side finding.** l3/gen_rows.c's "untwisted control" is ᾱT, not μT; it has no pole, hence no inflation. j4var.py reproduces both
  gen_rows ratios exactly.

## 8. Can the quasi-RH that THEOREM_ZF requires be supplied? (09:20) — HEURISTIC
* **What is needed.** (L³) implies zero-freeness in Re s > 7/8 for L(s, (k₀/·)₄λ^e) over Q(i). Base change from K = Q(ζ₁₂) reduces this
  to Hecke L-functions over K, since L(s, ξ∘N_{K/Q(i)}) = L(s,ξ)·L(s,ξω).
* **OpenAI's sextic argument over K needs:**
  * K ⊇ μ₆. This holds.
  * Sextic Gauss sums turning into cubic Gauss sums under Poisson. This is algebraic and works over any K ⊇ μ₆.
  * Explicit cubic theta coefficients over K. Kazhdan–Patterson uniqueness covers n = r+1 = 3 for GL₂ (for suitable covers); Patterson's
    explicit formula is over Q(ω), and the general-K formula is not checked here.
  * A quadratic large sieve over K: Goldmakher–Louvel, Math. Proc. Camb. Phil. Soc. 2012 (arXiv 1112.1642).
  * The 7/8 bookkeeping. In norm terms it looks degree-independent, but this is not checked.
* **Conclusion.** The quasi-RH that (L³) needs looks reachable in principle by a field-generalised OpenAI argument, which would be a
  large project.
* **It would not suffice.** Under full GRH the known arguments still give only Q^{3/2} (L3_REPORT §10b, §5 above), because the
  4th-power-free frequencies need a ratios-type (beyond-GRH) input. So (L³) needs both:
  1. a quasi-RH for the quartic Hecke family over Q(i) (necessary, by THEOREM_ZF);
  2. power-saving control of the family average of 1/(L·L̄) beyond GRH.
  Hypothesis (R) of L3_REPORT §10 supplies 2 and is sufficient. It is stated for 1/L on Re s ≥ ½ + ε, so it already contains zero-free
  information on average.

## 9. Bochner hulls for the second-moment series: GRH input and the correct functional-equation structure (09:15–09:25)

* **First attempt (l3j/tits_grh.py).** I added a GRH region for the Möbius side to l3h/tits.py's model. The hull then reaches the (L³)
  point exactly. But the certificate (tits_cert.py) is a single point σ₃(apex), so the result is an artifact.
* **What is wrong with the model.** tits.py, and BELOW_3_2 §7, treat the Poisson map σ₃ = P as a functional equation of one function.
  In fact P relates two different series:
  * Z_T(u), the T-side series;
  * Z_μ(v) = Σ_k N(k)^{−½−v₃}/(L(½+v₁,ξ_k)L(½+v₂,ξ̄_k)), the Möbius-side series.
  The functional equations are σ₁, σ₂ (Kubota) on Z_T and J on Z_μ, the latter with a class twist λ(k). The T-side group is
  ⟨σ₁, σ₂, J_T = P⁻¹JP⟩, the group of L3_REPORT §10c, not (2,6,6).
* **Second attempt (l3j/hull_groupoid.py).** With this correct structure and exact functional equations granted, the hull is unbounded
  (every u₃ is reached), with or without GRH. That is also an artifact, for two reasons.
  1. The diagonal pole at w = 1 has infinitely many images under the infinite group, and these polar hyperplanes accumulate. Bochner
     needs the function made holomorphic by multiplying out the poles, which is impossible here.
  2. For non-squarefree m and k the conductors are smaller than N(m)³ and N(k). The correction factors this forces converge only in
     bounded regions, so the functional equations are not exact far from the convergence regions.
* **Conclusion.**
  * Convex-hull computations that grant exact functional equations do not settle this problem in either direction. That includes
    BELOW_3_2 §7's u₃ ≥ 2, which used the wrong group.
  * A meaningful version must track the polar hyperplanes and the domains where the corrections converge. That is exactly the
    unsolved "infinite-type multiple Dirichlet series" problem over number fields.
  * No conditional or unconditional gain comes from this route.

## Status at 09:05Z
No proof of (L³). New this block:
* §2 and THEOREM_ZF:
  * The (L³) ⟹ zero-free implication is now argued with positivity, and strengthened. Any bound below Q^{3/2} forces a zero-free strip
    for every quartic Hecke L-function over Q(i).
  * So a proof of (L³) has to contain a new quasi-Riemann hypothesis (7/8) for that family. OpenAI's theorem does not cover it.
* The naive 4th-power counts at M < 1 are inconsistent with Kubota reflection, so they are not used.
* §1(c) identifies a second missing input: quartic theta coefficients at non-square indices, entering through root numbers in the
  dual variable.

## Status at 09:40Z (end of this line of attack)
* THEOREM_ZF (+ THEOREM_ZF_PROOF, lemma3test, j4pred/j4obs, signcheck).
  * Assume the squarefree (G) form of (L³), uniform in index twists, with any power saving below Q^{3/2}. Then every quartic Hecke
    L-function over Q(i) is zero-free in Re s > (7+10θ′)/(8+8θ′); (L³) itself gives 7/8.
  * The optimal-shape square-case large sieve for T (L3_REPORT §12 sub-problem 1) gives 7/8 as well.
* What a proof of (L³) must therefore contain:
  1. a new quasi-Riemann hypothesis for the quartic Hecke family over Q(i), not covered by OpenAI's theorem;
  2. beyond-GRH control of the 4th-power-free frequencies. GRH alone, by every known argument, stops at Q^{3/2} (§5, L3_REPORT §10b).
* Open technical gap, recorded: the passage from ψ_α (all c) to the squarefree (G) form (local 2×2 inversion, THEOREM_ZF_PROOF Lemma 1).
* (L³) at Q^{1+ε}: OPEN.
