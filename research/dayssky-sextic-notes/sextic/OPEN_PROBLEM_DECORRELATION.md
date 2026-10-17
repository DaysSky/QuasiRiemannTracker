# The decorrelation estimate behind the next frontier move

DaysSky, with Claude assistance. 17 Oct 2026. Frontier (unconditional): θ = 0.874836 (Gao, GTimG/QuasiRiemannTracker PR #62).

## Setting

- K = Q(√−3), O = Z[ω], with sextic residue symbols χ_n(u) = (u/n)₆ for primary n.
- Rows: sixth-power-free u with N(u) ≍ U.
- Plain sums: S_i(u) = X^{−1/2} Σ_n W_i(N(n)/X) χ_n(u), with X = U^{0.417} (the record's binding plain length). Write P(u) = S₁(u)S₂(u), a product of length D = X² = U^{0.834}.
- Slot moduli: primes m with N(m) = Z = U^z, z ∈ [0.04, 0.08] at the binding point.

## The estimate (averaged thin family)

    Σ_{m ~ Z prime} Σ_{N(v) ≤ U/Z²} |P(v m²)|²  ≪  U^ε · Z^{1+c} · (U/Z²)        (★)

That is, natural size up to a loss Z^c.

- c = 0 gives plain slot cost 2. Under the certified ladder (CONDITIONAL_THEOREMS.md, row "plain slot weight 2") that is θ ≤ 0.8734, a gain of about 1.4·10⁻³.
- Any c < 0.36 beats the paper's 6κ ≈ 4.5 at the binding point.
- Why the average suffices: the off-diagonal in §92 of DH_SEXTIC.md is (U/(DZ))·Σ_mΣ_h|B_{hm⁴}|². Cauchy over m inside the h-sum is already in place.

## Equivalent forms (all PROVED equivalences, DH_SEXTIC.md §98 addenda and §104)

1. **Twist form.** χ_n(vm²) = χ_n(v)·ξ_m(n), where ξ_m(n) = χ_n(m)² is the cubic character of conductor m. (★) is the plain fourth moment of the full family over v ≤ U/Z², with coefficients twisted by ξ_m, averaged over m.
2. **Decorrelation form.** Let T(ψ) = E_v|Σ_N c_N ψ(N) χ_N(v)|², where c = W₁∗W₂. Then T(ψ) = Diag + Off(ψ), with Off(ψ) = Σ_t F(t)ψ(t) and

       F(t) = Σ_{N≠N′, N≡tN′ (mod m)} c_N c̄_{N′} Σ_v χ_N χ̄_{N′}(v).

   (★) holds if F(t) = o(Diag·Z^{c−1}) for t ≢ 1, or more weakly if the Fourier mass of F at the characters of order dividing 6 is small.
3. **Twisted fourth moment.** Σ_v ∫|L(½+it, χ_v ξ_m)|⁴ dt, with conductor dependence N(m)^{c}.
4. **Hadamard form.** c*(A∘B)c, with A the row Gram matrix and B the slot Gram matrix (rank ≤ Z).

## Numerical evidence (mmnum/)

- thin_1e6.log and thin_3e6.log: the ratio thin/full is 0.77–1.25 up to N(m) = U^{0.355}. The fitted loss exponent is 0.05–0.07 at U = 3·10⁶ and falls as U grows.
- twistall_1e6.log: for every character ψ mod m (N(m) = 7, 13, 19), T(ψ) equals the untwisted value to ±0.002. That includes the foreign characters of order 4, 9, 12 and 18. So F(t) ≈ 0 for t ≠ 1, and nothing is special about the thin twists in the data.

## Why known methods stop (PROVED, with references to DH_SEXTIC.md)

- **General coefficients cannot work.** Coefficients aligned with one dual row produce sixth-power copies, so the general-coefficient operator norm of (★) exceeds natural by W^{1/6}/Z = U^{0.06–0.10}, with W = D²V^{−1}. This is the shape of de Faveri's optimal sieve (arXiv 2610.04045). Any proof must use the smoothness of c.
- **Lemma 18.1's smoothness recursion pays N(m)^1 (the twist rule).** The twist acts as a row shift at every level (χ_a(m²)G(a,h) = G(a, hm⁴)). Only the terminal Hecke reflection, which is unavoidable because A = 0.834 > 5M/6, needs the conductor N(m).
- **Parity.** The slot exponent is 2i, always even, so the twist is cubic. Exponent 1 (rows vm) would be free by merging. De Faveri's cubic Prop. 8.1 merges because −2 ≡ 1 (mod 3). In the sextic case the quadratic component ε = χ³ survives every Poisson step.
- **Other routes, each circular:** Parseval over characters mod m, the conjugate grouping, "thin = full − coprime" asymptotics, lower bounds through test polynomials, Schur/Hadamard, a structured de Faveri recursion, hybrid t-averaging, Halász/Huxley.
- **Literature (Oct 2026), no input:**
  - de Faveri 2610.04045: optimal sieve; its Theorem 1.2 is cubic and lives at Q = √(AB).
  - de Faveri–Dunn–Hoffstein 2607.07911: the cubic sieve is sharp.
  - Diaconu–Ion–Paşol–Popa 2607.27131: second moments only, and needs 2r-th roots of unity.
  - Dunn 2403.13151: cusp forms only.

## What a proof needs

A version of Lemma 18.1, the plain fourth moment of the sextic family (natural for all lengths), for coefficients that are smooth on residue classes to a small moving modulus m, with the cost in N(m) below N(m)^{0.36}. For Dirichlet characters the analogue is the uniform twisted fourth moment (Blomer–Fouvry–Kowalski–Michel–Milićević, via Kloosterman/spectral methods). For the sextic family that would need a metaplectic Kuznetsov formula on the 6-fold cover with Kloosterman-type bilinear bounds. That has not been developed.

**Spectral shape of a proof (HEURISTIC).** The second transform of the twisted Gauss-row norm produces sums over moduli u, v of Gauss–Gauss correlations F(u, v; j)·χ̄_u(m⁴)χ_v(m⁴): metaplectic Kloosterman-type sums at level about m. A Kuznetsov formula on the cubic cover at level m would split them three ways:
- **Residual spectrum.** This is the theta function, twisted by ξ_m. Its size is damped by N(m)^{−1/6}·N(m)/(N(m)+1) (measured: mugauss_t1/twist2_1e5.log, to 2%). That is the source of the Gauss-sum bias and of the de Faveri–Dunn–Hoffstein lower bound.
- **Continuous spectrum.** Eisenstein series, expected square-root size.
- **Cuspidal spectrum.** This needs a spectral large sieve uniform in the level m. The spectral density grows like N(m), while the normalised Fourier coefficients shrink like N(m)^{−1/2}, so a Deshouillers–Iwaniec-type bound would be level-neutral, as in the Dirichlet case (BFKMM).

If such a level-uniform metaplectic spectral large sieve held, the cuspidal and continuous parts would be level-neutral. The residual part is a twisted theta term, which in the data is smaller than its untwisted counterpart, and the untwisted counterpart is already inside Lemma 18.1's natural bound. So the twisted norm would plausibly come out natural up to N(m)^{o(1)}. That last step is a heuristic count, not a computation. Status: the cubic metaplectic Kuznetsov formula exists at fixed level (Livné–Patterson, Proskurin, Louvel, Dunn), but no level-uniform version is known. This is the most concrete research route to (★) found so far.

## Second open input (inverse slot), same flavour

Marked RD: orthogonality between distinct sublattice shifts hm⁴ and hm′⁴ of the Möbius dual (PLAN_RD.md). c_I = 7/4 would give 0.87425. The source of c_I = 2 is the reflection conductor (DH_SEXTIC.md §32). Gao's affine caps 13/400 and 3/200 are these same two costs (VariableSupply.lean).

**Literature check for the spectral route (17 Oct, web search).** The Bruggeman–Kuznetsov formula for the cubic cover over Q(ω) is in Louvel, arXiv 1107.0197 (residual spectrum = Kubota's theta term).
It is used at the lowest level Γ₂ = ⟨SL₂(Z), Γ₁(3)⟩ in Dunn, arXiv 2403.13151 (Algebra & Number Theory 2025), whose bounds depend on the progression modulus as N(v)⁴. No level-aspect (Γ₁(3q))
uniform version, and no level-uniform metaplectic spectral large sieve, was found. So the spectral route needs that theory built first.

**Correction and refinement of the spectral route (07:36Z).** Lemma 13.4 shows that the second transform produces characters, F(D²a, E²b; j) ∝ χ_a(j)χ_b(−j), not Kloosterman sums. So the
recursion stays inside character families, and Kuznetsov enters only through smooth sums of Gauss sums over the modulus. With the twist absorbed, Σ_a A(a)χ̄_a(m⁴)G(a, h) = Σ_a A(a)G(a, hm⁴),
the twisted Gauss-row norm is a mean square of Kubota-type sums along the sublattice of numerators hm⁴. For a single smooth A this is
  Σ_{N(h)≤H}|ψ₆(hm⁴, s)|²,
i.e. de Faveri–Dunn–Hoffstein's Lindelöf-on-average bound (their (2.9)) restricted to a sublattice. That is a level-m⁴ Rankin–Selberg problem for the metaplectic Eisenstein series, which they treat by
Zagier regularisation at level 1. Classical analogues (Σ_{n≤H, q|n}|a_f(n)|² = cH/q + O(H^θ q^β)) have power-saving errors with polynomial level dependence. With H ≈ U^{0.67} and q = m⁴ ≈ Z⁴,
θ = 3/5 and β = 1 would cover z ≤ 0.033; better exponents would cover more. But the record's A = W₁∗W₂ is a convolution, so the relevant object is the A₂-type sextic double Dirichlet series
Σ G(a₁a₂, h)… (Brubaker–Bump–Friedberg Weyl group multiple Dirichlet series), and its mean square along sublattices is unstudied. Programme:
(i) sublattice Rankin–Selberg for the n = 6 metaplectic Eisenstein series, with explicit level dependence;
(ii) the same for the A₂ double series;
(iii) insert into §92 step 4.
This is multi-session research.

**Error budget for step (i) (HEURISTIC, 07:37Z).** If the level-m Rankin–Selberg gives Σ_{N(h)≤H}|ψ_χ(h, ½+it)|² = (natural main term ≍ H) + O(H^{1/2+ε}·N(m⁴)^{c′}), then with H ≈ U^{0.67} the error fits
inside the allowed Z^{0.36} for z ≤ 0.084/c′. So c′ ≤ 1 covers the binding slot range z ∈ [0.04, 0.08]. The record's plain witness is a product W₁∗W₂, which the fourth-moment count needs: a single
plain sum of length 0.834 has in-bin excess at most δ·0.166 and is a poor witness. So the A₂ double series version is the one actually required.

## Sextic Kubota numerics for the spectral route (kub6/, 17 Oct, after 07:40Z)
- **Implementation (NUMERICAL, checked against brute force).** Sextic Gauss sums g(h, a) = Σ_{x mod a}χ_a(x)e(Tr(hx/(a√−3))) over Z[ω]. Checks: |g|² = N(a) at prime and composite moduli,
  g(h, π) = conj((h/π)₆)·g(1, π), and twisted multiplicativity g(h, ab) = (b/a)₆(a/b)₆·g(h, a)g(h, b). The factor (b/a)₆(a/b)₆ is (a/b)₆² times a sign: sextic reciprocity has a quadratic sign.
- **Sublattice mean squares (kub6/sublat.py).** M(m) = mean_{N(h)≤H}|Σ_a w(a)g̃(hm⁴, a)|² against the untwisted M₀ over moduli coprime to m:
  | regime | weights | M₀/Σ\|w\|² | M(m)/M₀, N(m) = 7…61 | random-twist control |
  |---|---|---|---|---|
  | rows-dominant (Y = H = 3000) | single | 0.97–0.99 | 0.95–0.99 | 0.96–0.99 |
  | column-dominant (Y = 4·10⁴, H = 400) | single | 1.08 | 0.84–0.87 | 0.85–0.90 |
  | column-dominant (Y = 4·10⁴, H = 400) | product W₁∗W₂ | 0.92 | 0.91–0.96 | 0.92–0.98 |
  So the twisted (sublattice) Gauss-row norm is natural in the record's column-dominant regime, for the product (A₂) weights too. The untwisted single-weight excess 1.08 is the theta bias, and
  twisting removes it.
- **Single weights: an analytic bound (PROVED modulo standard Mellin/AFE bookkeeping).** For smooth w,
  Σ_a w(N(a)/Y)g̃(h, a) = (theta residue at s = 2/3) + Y^{1/2}∫ψ₆(h, ½+it)ŵ dt.
  The h-mean square of ψ₆(h, ½+it) is Lindelöf on average, by the sextic sieve at the balanced point: Θ₆(H, H^{1/2}) ≪ H^{1+ε} under BGL or de Faveri. So the untwisted single-weight norm is
  (HY)^{1+ε} plus the theta term, i.e. natural. For the sublattice, ψ₆(hm⁴, s) has AFE length (HZ⁴)^{1/2}, and de Faveri's Θ₆(H, H^{1/2}Z²) = H + H^{1/2}Z² + HZ^{2/3} + H^{3/4}Z^{5/3} gives loss
  Z^{2/3} for Z ≤ H^{1/4}. That beats the twist rule's Z, but not the Z^{0.36} the record needs, and it is for single weights only.
- **Why products are the crux.** B(k) = Σ_{n₁,n₂}W W g̃(k, n₁n₂) = Σ W W·(n₂/n₁)₆(n₁/n₂)₆·g̃(k, n₁)g̃(k, n₂). The cubic coupling between the two plain variables stops B from factoring, and a
  pointwise cubic-sieve bound on the coupling loses X^{2/3}. The coupling is invisible on the diagonal (|c′_N|² = |W∗W(N)|²); it acts only in the off-diagonal, which is where (★) lives.

**Twisted-pair route checked (PROVED bookkeeping; closed).** Target (§83–84 reduction): T(m, m′) = Σ_u W|S₁S₂(u)|²ψ_u(m)ψ̄_u(m′) ≪ U·N(mm′)^{−θ′} with θ′ > 0.318 (0.325 against Λ = 3.96).
- The zero frequency is a shifted diagonal Σ_k c_{m′k}c̄_{mk}·U ≍ U/Z = U·N(mm′)^{−1/2}, which already has θ′ = 1/2. So everything rests on the nonzero frequencies.
- **Tempting step:** choose W with Ŵ ≥ 0, localise the kernel in fine boxes of N, and apply Cauchy–Schwarz in h, hoping the self-terms are the untwisted off-diagonal Σ_u W|P|² − Diag. If they
  were, a power-saving asymptotic Off = O(U^{1−η}) with η > 0.64z would give input 2.
- **Why it fails:** the self-term Σ_hŴ|C_m^B(h)|² is built from G(Nm, h) = (m/N)₆(N/m)₆G(N, h)G(m, h). It equals N(m)·(the thin-twisted Gauss-row norm), not the Poisson dual of Σ_u W|P_m|².
  For the common factor m on both sides, Poisson modulo Nm·N′m keeps only h ≡ 0 (mod m) and does not factor as G·Ḡ. So the Cauchy self-terms are the thin families again.

**Convexity is enough (HEURISTIC conductor count).** Model: Σ_{N(h)≤H}|F(hm⁴)|² = (natural main ≍ H) + O(H^{1/2}·cond^{1/4}), with cond = N(m)^k for a degree-k Rankin–Selberg L-function at level m and
H = U^{0.668+2z} the dual range. The error stays inside Z^c for z ≤ 0.334/(k/4 − c − 1). Results:
- single Kubota series (GL₂ cover, k = 4): every z;
- A₂ double series (GL₃ cover, Brubaker–Bump–Friedberg, k = 9): z ≤ 0.267 (c = 0) or 0.375 (c = 0.36).
The binding range is z = 0.04–0.08, so convexity-quality Rankin–Selberg on the metaplectic GL₂/GL₃ covers with polynomial level dependence would make (★) natural, giving plain slot cost 2 and
θ ≤ 0.8734.
**Why the paper's elementary induction cannot mimic this (structural).** Rankin–Selberg reflects only the error term, so the level enters as cond^{1/4} beside a main term untouched at w = 1.
Lemma 18.1's induction is an ε-loss (Lindelöf-on-average) induction with no main/error split. Each child's width drop σ is paid by its excess σ, which reflection absorbs at the child's full
conductor, so the twist's N(m) multiplies the main size. An asymptotic version of Lemma 18.1 (main term plus a power-saving error at every level) would let N(m) fall on the error only. That would
be the elementary replacement for the spectral route, and the paper's tight induction does not provide it.

**Classical analogue located (CITED).** Petrow–Young, arXiv 1908.10346: Lindelöf-on-average fourth moment of Dirichlet L-functions along a coset of the subgroup of characters mod d (q*|d), via
Conrey–Iwaniec/Kuznetsov, leading to the Weyl bound. For Dirichlet characters a fixed twist ξ of a full family is free, because the large sieve is optimal for general coefficients
(Σ_{q≤V}Σ_χ|Σa_nχξ(n)|² ≤ (V² + N)‖a‖²). The sextic family's general-coefficient sieve carries the copies/theta loss, so (★) needs either an automorphic (metaplectic Kuznetsov/Rankin–Selberg)
argument, as in Petrow–Young, or structure à la Lemma 18.1. Also seen: Gao–Zhao, shifted moments of cubic and quartic L-functions (Monatsh. Math. 2026), which is GRH-conditional and not an input.
**Slots on other legs (PROVED exponent counts; no gain).**
- Long leg: copies preserve slots (Q_{ua⁶} = Q_u), but the slot costs 5/6 per unit against a saving of 2q = 0.466, so long(t) + z(5/6 − 2q) only gets worse.
- Inverse leg: pointwise deletion of slot length beyond (1 − r)/2 costs κ per unit against a saving of 2q, a net κ − 2q = 0.28 > 0, so it also gets worse.
