# Research plan: the rows-dominant asymptotic (RD) behind inputs 1 and 3

DaysSky, with Claude assistance. 17 Oct 2026. Reductions: DH_SEXTIC.md §26, §32, §69, §70, §88, §90; values: CONDITIONAL_THEOREMS.md.

## Target
For the canonical family of OpenAI Lemma 17.2 (rows k, N(k) ≤ Z^M; squarefree columns n with coefficients ᾱγ₂(n)·d(n), where γ₂ is the normalised cubic
Gauss sum and d is a product-form slot mark of length z₀), equivalently for the dual sieve of §69:

  Σ_{N(h)≤K} |Σ_{N(n)≤D} c_n χ̄_n(h) W(N(n)/D)|² = K·‖c‖²·𝔖 + O(U^ε D²),   c_n = ε(n mod 4)(n/|n|)conj(g̃₃(n)),

with K = D²/U, in the rows-dominant range D < K ≤ D^{1+c₀}. The main term is exact (𝔖 is the coprimality factor) and the error is a power below it.

## What each range is worth (certified, CONDITIONAL_THEOREMS.md)
| c₀ (LM slope 5/(6(1+c₀))) | θ ≤ |
|---|---|
| 1/24 (slope 4/5) | 0.87462 |
| 1/9 (slope 3/4) | 0.87425 |
| 0.19 (slope 7/10) | 0.87368 |
With marks (M up to F₀ − (c_I − 1)z₀): c_I = 3/2 ⟹ 0.87325; c_I ≤ 6/5 ⟹ 0.8728. Together with plain weight 2 ⟹ 0.8708.
If RD_{c₀} holds for the marked family with the same relative extension, M ≤ (1 + c₀)(F₀ − z₀), then at the binding r ≈ 0.75, z₀ ≈ 0.125 it gives
c_I ≈ 2 − c₀(F₀ − z₀)/z₀ ≈ 2 − 5.8c₀. So c₀ = 1/9 gives slope 3/4 and c_I ≈ 1.36 together; the pair (c_I = 3/2, slope 3/4) is a conservative representative
(certified: θ ≤ 0.8728 at h = 0.7957, CONDITIONAL_THEOREMS.md row 5′).

## Classical-strength checkpoint (17 Oct, 04:55Z)
- **Marked range needed.** At the binding r = 3/4 the marked dual has length J = N²Z²/U = U^{2r+2z−1}. For c_I = 2 this is J = N, exactly self-dual. For c_I = 7/4 it is J = U^{r+(1−r)/7}, so
  N^{1+c₀} with c₀ = (1−r)/(7r) ≈ 0.048. For c_I = 1 it is J = U, so c₀ = (1−r)/r = 1/3. This agrees with c_I ≈ 2 − 5.8c₀ above.
- **BFI strength.** Bombieri–Friedlander–Iwaniec reach level x^{4/7} for well-factorable weights, where the classical barrier is x^{1/2}: a factor 8/7 in the level exponent. The same factor on the
  large-sieve barrier q ≤ U/N for the twisted Möbius mean square (FOUR_INPUTS input 3) gives 2z ≤ (8/7)(1 − r), i.e. c_I = 7/4. Its certified value is in CONDITIONAL_THEOREMS.md (row 2″).
- So **RD with c₀ ≈ 1/24 for both the unmarked and the marked family** gives LM slope 4/5 and c_I ≈ 7/4 together. That is the classical-strength target: **certified θ ≤ 0.87405** (CONDITIONAL_THEOREMS.md row 2‴). At Maynard's strength (x^{3/5}, factor 6/5): c_I = 5/3, marked c₀ ≈ 0.067, giving θ ≤ 0.8740 alone and 0.8738 with slope 4/5 (rows 2ᵐ, 2⁗). At Lichtman's (x^{66/107}): c_I = 107/66, giving 0.87382 and 0.87362 (rows 2ˡ, 2ˡ′).

**Why BFI does not transfer directly (§94).** The dispersion needs a smooth variable. In a Heath-Brown decomposition of μ the smooth factor carries cubic Gauss sums, whose smooth sums have the Patterson main term M^{5/6}. So the theta main terms of the Type I pieces must be shown to cancel across the decomposition first. That is step 2 below. Numerically (§94 addenda) the theta coefficient on multiples of d is Π_{p|d}(N(p)+1)^{−1}, so Möbius kills it through 1/ζ_{Q(ω)}(1 + w). Step 2 is therefore a uniform-in-d power-saving error term for the twisted Kubota series Σ_m g̃(dm)N(m)^{−s} past the theta pole. Its pointwise shadow, a power saving below X^{5/6} for Σμg̃χ, already needs the Type III sums (A = B = C = X^{1/3}) that Dunn–Radziwiłł call self-dual and open (§95). The theta rows themselves are harmless for RD.

## Why the paper's tools stop (proved)
- **Self-duality (§70, re-derived 17 Oct).** Poisson in h maps the off-diagonal to exactly (K/D)·Σ_{j≤U}|M_j(D)|², the original LM sum. Since c_n = μ(n)g̃₆(n), the
  dual coefficient is μ(n)|g̃₆(n)|²χ̄_n(j). No Poisson step can help.
- **Reflection lengthens the dual past the self-dual point (§14.2, §26).** The canonical recursion keeps M ≤ F₀ − z₀ and gives only the diagonal bound beyond it.
- **Large values and pointwise bounds.** Pointwise, the small-conductor sums Σ_n μ(n)g̃₆(n)χ̄_n(k₀) need D^{0.773} against convexity's D^{0.875} (§79). Lindelöf-on-average bounds for Kubota's
  series carry an ε and give no exact main term.
- **General-coefficient sieves cannot do it.** The cubic and quartic large sieves are not orthogonal (de Faveri–Dunn–Hoffstein), and the extremal configurations are sixth-power copies.
  RD must use that c_n is bias-free (Möbius kills the Patterson pole; numerics in §69 and §76).

## Proposed route (research programme)
1. **Dispersion with well-factorable moduli (BFI-type).** In the marked form the off-diagonal is Σ_{p,p′}a_pā_{p′}T(p, p′), with T(p, p′) twisted Möbius mean squares over
   progressions mod pp′ (§88). Split the moduli and apply Cauchy–Schwarz so that one prime moves outside, as in Bombieri–Friedlander–Iwaniec. Then apply Poisson to reach sums of
   metaplectic Kloosterman sums with smooth moduli.
2. **Kuznetsov on the cubic cover.** Proskurin and Livné–Patterson give spectral expansions of sums of cubic Kloosterman/Gauss sums over moduli. The cuspidal part needs a
   spectral large sieve (Dunn, arXiv 2403.13151). The residual part is the cubic theta. Möbius should remove it (no bias, §69), but this has to be proved
   (a "Möbius orthogonality to theta" statement). For prime-supported sums the bias is real: Dunn–Radziwiłł (arXiv 2109.07463) prove Patterson's
   Σ_{p≤X} g̃(p) ~ cX^{5/6}/log X. For μ-weighted sums the numerics show none (|Σμg̃| ≤ 0.63√X up to 2·10⁵, §69). There is no Euler product (twisted
   multiplicativity), so the usual 1/L argument is not available.
3. **Spectral gap.** A fixed gap for cubic metaplectic Maass forms on Kubota's Γ(3) ⊂ SL₂(Z[ω]) away from the theta eigenvalue, uniform in the level (pp′)², is needed for the power saving.
   None was found in the literature (17 Oct search).
4. **Checkpoint.** c₀ = 1/24 (LM slope 4/5) already beats the frontier (0.87462, certified); c₀ = 1/9 gives 0.87425.

## New direction (§96): a density-estimate loop
DR's corrected cubic sieve (1.13) reduces the cubic sieve for μ-supported sequences to zero counts of cubic Hecke L-functions. A density estimate N_cub(5/6; H) ≪ H^{2/3−η} for the cubic family
beats Heath-Brown's (AB)^{2/3} for μ. Heath-Brown-derived density estimates meet the threshold exactly at σ = 5/6 (no gain), so the input must come from a detector stronger than Montgomery's, such as the OpenAI-type architecture run on
the cubic family. It would improve the critical cubic sieve for μ. But density estimates alone cannot give RD (the §96 addendum: zero counts give LM slope (2σ − 1) + E/(s − 1)), so this loop is a side road, not the route to RD.

## Evidence
The dual sieve sits on the coprimality plateau 0.94 for K up to 32D at P = 3·10³, 10⁴ and 3·10⁴ (up to 386 939 rows), within 1% of a random control, while the unweighted Gauss sequence climbs to 1.18 (§76 addendum). The primal LM is flat up to U = 10⁶ (§41, §57). Marked slot costs look natural (§31). The Möbius fourth moment is
flat too (§88 addendum 2), which by §88 would make slots free outright.

## The single inequality to prove for c_I = 7/4 (θ ≤ 0.87425), stated at the binding point
With U the family size, N = U^r (r = 0.745), slot primes m, m′ of norm ≍ Z = U^z (z = (1 − r)/(7/4) ≈ 0.146), H = N²Z²/U = N^{1+c₀} (c₀ = (1 − r)/(7r) ≈ 0.049), and
B(k) = Σ_{N(n)~N} μ(n)g̃₆(k, n)W(N(n)/N), the needed statement is, for distinct slot primes (on average over the slot weights a_m with |a_m| ≤ 1):
  Σ_{m≠m′} a_m ā_{m′} Σ_{N(h)≤H} g̃₆(h, m)·conj(g̃₆(h, m′))·B(hm⁴)·conj(B(hm′⁴)) ≪ U^{−η}·Z·H·N   for some η > 0.
The right side is U^{−η} times the trivial (diagonal) bound, and the left side is the off-diagonal of the marked dual (§92, §98). The diagonal m = m′ and all large-sieve upper bounds are already available.
Numerically the marked moment is at natural size up to r + 2z = 1.4 (§31). At the binding r = 0.745 itself (sampled rows, U = 10⁹, mmnum/slotbig_1e9.log), ρ = 1.00–1.06 at r + 2z = 1.145 and 1.245, beyond the c_I = 2 line and far beyond what this needs (r + (7/4)z = 1).

## New leads for the trilinear step (CITED, 17 Oct 06:07Z search)
- **Non-abelian amplification for bilinear forms with Kloosterman sums** (arXiv 2511.08445) and Milićević–Qin–Wu (arXiv 2511.07550, Nov 2025) give power-saving bounds for bilinear and trilinear forms
  in Kloosterman sums in ranges where Weil plus Cauchy–Schwarz gives nothing. Our open step is the metaplectic analogue: a trilinear form in (h, n, m) with cubic Gauss sums g̃₆(h, m), g̃₆(hm⁴, n).
  A metaplectic version of their amplification is the most specific new tool to try.
- **Dunn, "Metaplectic cusp forms and the large sieve"** (arXiv 2403.13151, ANT 19:9, 2025): power saving for cubic metaplectic cusp-form coefficients over primes (bilinear estimate via circle method and Voronoi).
  This is the cuspidal half of the problem. Möbius kills the residual (theta) half in the main term (§94), so combining the two is the concrete programme for PLAN_RD step 2.
- **Fit check (abstract read 06:09Z).** 2511.08445 bounds Type II sums of Kloosterman sums with composite modulus c = pp′ (p ≍ p′), length √c, saving c^{−1/12}, via Fourier analysis on SL₂(Z/cZ). Our slot moduli
  are exactly products mm′ of two same-size primes. But in the marked dual the mm′-dependence is the abelian character χ̄_mχ_{m′}(h) times B(hm⁴)·conj(B(hm′⁴)). Poisson in h returns the original family (§102),
  and no Kloosterman sums appear. Their method would enter only after a dispersion step in the slot moduli that produces metaplectic Kloosterman sums of modulus mm′. The bridge to build is therefore
  "BFI-type dispersion in the slot moduli, then a metaplectic analogue of SL₂(Z/cZ) non-abelian amplification".
- **Bridge found (§103).** Factor the row variable u = u₁u₂ (hyperbola decomposition; ψ_u(n) is multiplicative in u). The progression u ≡ b (mm′) then becomes u₂ ≡ b·ū₁, and Poisson in u₂ splits by CRT into
  Gauss sums mod nn′ (the Möbius dual) and Kloosterman-type sums mod c = mm′. After Cauchy–Schwarz this is a Type II form in Kloosterman sums of composite modulus with two same-size prime factors,
  the exact setting of 2511.08445 / 2511.07550. Checks (a)–(c) in §103 decide whether this gives c_I < 2 unconditionally.
  First pass (§103): the (n, n′) sum recombines into a coupling weight |B(h·u₁*)|², so the form is a twisted Type II Kloosterman form. The sixth-power row factorisation avoids the coupling but fails on lengths.
  **Target estimate:** Σ_{u₁~U₁, h~H} α(u₁)|B(h·u₁*)|²·S(h, b·ū₁; mm′) ≪ (mm′)^{−η}·(trivial), on average over slot primes m ≍ m′ ≍ Z.
