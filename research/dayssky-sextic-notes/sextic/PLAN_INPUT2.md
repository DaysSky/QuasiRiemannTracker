# Research plan: the plain slot cost (input 2) via a twisted fourth moment on the cubic cover

DaysSky, with Claude assistance. 17 Oct 2026. Supporting reductions: DH_SEXTIC.md §31, §78–§84; FOUR_INPUTS.md.

## Target and value
**Twisted plain fourth moment (TPF).** For slot primes m ≠ m′ of norm ≤ Z = U^z, with z ≤ z₀ and n₁ + n₂ + 2z ≤ 1:

  T(m, m′) = Σ_{N(u)≍U} |S_u(U^{n₁})S_u(U^{n₂})|² ψ_u(m)ψ̄_u(m′) ≪ U^{1+ε}·N(mm′)^{−θ′}.

(S_u is the centrally normalised plain block and ψ_u(n) = (u/n)₆. Rows are all sixth-power-free u, or a smooth lattice weight with V̂ ≥ 0.)

**What each θ′ buys (MODEL, §83).**
- θ′ > 0.32: slots in the plain leg beat the paper's 6κ (and > 0.325 beats Theorem 2.2's Λ = 3.96).
- θ′ = 1/2: the plain slot weight becomes 2. On the current record pipeline (Gao's constants, §83 addendum 4) that is θ = 0.873382, i.e. −1.45·10⁻³
  against the frontier 0.874836. Weight 3 gives 0.874313 and natural weight 1 gives 0.870524.
- The same input feeds the HIGH side through ϑ ≤ q/c_P (§81).

## Why existing tools stop (§83–§84)
1. Cauchy–Schwarz on the dual cross-correlation reproduces 6: 2 from Poisson normalisation and 4 from the sublattice h·m⁴.
2. Without Cauchy, T has no main term in either representation, and Poisson maps it to itself (self-duality).
3. Lemma 18.1's twist rule (bound Z^{m+q}) is conductor-sized for the mildly thin subfamily {u ≡ 0 mod m²}, so it loses exactly Z.
4. Multiple Dirichlet series give affine D̃₄: continuation to a half-space only.

## Correction (02:08Z): the right analogue is the quadratic family, where this is open
The sextic family is thin: Poisson in u gives Gauss sums, not a congruence a ≡ b (mod q). So the shifted-convolution and Kloosterman mechanism of the
Dirichlet-family proof below does not arise. The matching family is quadratic characters {χ_d}: also thin, Poisson-dual to Gauss sums, with an affine D̃₄
fourth moment. There the first unconditional fourth-moment asymptotic (Shen–Stucky, arXiv 2402.01497, 2024, after Li's third moment) has error
O(X(log X)^{6+ε}) against a main term X(log X)^{10}. That is a log saving, and the authors say "one would need additional ideas to improve the error
term". A twisted version with a power saving in the twist modulus, which is what input 2 needs (θ′ > 0.32), is beyond the current state of the art even
for quadratic characters. The steps below remain a possible route, but they are a research programme, not an adaptation of known proofs.

## The thick-family analogue (for comparison only)
For primitive Dirichlet characters mod q: Σ_χ |L(1/2, χ)|⁴ χ(m)χ̄(m′) = main(m, m′) + O(q^{1−η}) with a power saving for mm′ ≤ q^{θ₀}. This is due
to Young (fourth moment), Blomer–Fouvry–Kowalski–Michel–Milićević (twisted, with the off-diagonal through Kloosterman sums and Kuznetsov) and
later refinements. The off-diagonal becomes shifted divisor sums Σ d(n)d(n + h), treated by the δ-method or Kuznetsov with Weil's bound. The power
saving comes from the spectral gap.

## Proposed route on the cubic/sextic cover
1. **Opening.** Expand |S₁S₂|² and execute the u-sum by Poisson on the lattice. The frequencies become sextic Gauss sums G(χ_{N₁}χ̄_{N₂}, h) with
   N₁ = nm and N₂ = n′m′ (§83 formula). Separate classes mod 36 so the reciprocity units factor.
2. **Metaplectic Kloosterman form.** Write Σ_{n,n′} c_n c̄_{n′} g̃6(h, nm)conj(g̃6(h, n′m′)) through twisted multiplicativity as sums of sextic Gauss sums
   in the modulus variable. Then apply Voronoi or reflection in n, using the plain coefficients' reflectability (Lemma 18.1's mechanism). The aim is
   sums of metaplectic Kloosterman sums K₆(a, b; c) = Σ_{d mod c} (c/d)₆ e((ad + bd̄)/c) with c running over a range of size about D/Z.
3. **Kuznetsov on the 6-fold (or cubic) cover.** Proskurin's general Kloosterman-sum Kuznetsov formula, specialised to Γ₁(36) ⊂ SL₂(Z[ω]) with the
   Kubota multiplier. The spectral side has metaplectic Maass cusp forms, the Eisenstein continuous spectrum, and the residual spectrum: the cubic/sextic
   theta functions, which carry the Patterson-type main terms.
4. **Residual (theta) terms.** Show they cancel or are lower order for the twist ψ_u(m)ψ̄_u(m′) with m ≠ m′. Heuristic (§83): shifted Θ rows k·m⁴ are
   generic, and Θ-copies contribute U^{z−1/6} relative, so they are harmless for z < 1/6.
5. **Cuspidal and continuous terms.** Bound them with a spectral large sieve on the cover (Dunn, arXiv 2403.13151, proves such bounds for cubic metaplectic
   cusp forms) and with the spectral gap of the metaplectic Laplacian away from the theta eigenvalue. Polynomial dependence on the level (mm′)⁴ gives
   θ′ = 1/2 for Z ≤ U^{z₀} with some z₀ > 0. Any z₀ ≥ 0.04 (the binding slot widths) already moves the frontier.

## Known ingredients
- Kubota's metaplectic Eisenstein series and Patterson's cubic theta coefficients, with explicit cusp expansions (Dunn–Radziwiłł, arXiv 2109.07463 §5
  and App. A, used by the OpenAI paper).
- Proskurin's Kuznetsov formula for general Kloosterman sums; Dunn's metaplectic spectral large sieve (arXiv 2403.13151).
- The A₂ Weyl group multiple Dirichlet series (Brubaker–Bump–Friedberg) for the plain-product dual coefficients (§84). These are useful for the
  main-term bookkeeping, though they do not give the power saving.

## Checkpoints
- Numerics (§83 addendum): |T(m, m′)| sits at the row-noise level, with θ_eff ≥ 0.6 in every feasible range up to U = 3·10⁶.
- A first milestone that already counts: θ′ > 0.32 at slot width z = 0.04 for n₁ + n₂ = 0.834.

## Literature status of the spectral steps (17 Oct, 02:00Z search)
- **Step 3 (Kuznetsov on the cover).** Louvel's dissertation, "Twisted Kloosterman sums and cubic exponential sums" (d-nb.info/101056661X), develops cubic
  metaplectic forms on SL₂(Z[ω]) classically and representation-theoretically, with two spectral summation formulas for cubic Kloosterman sums.
  Dunn (arXiv 2403.13151) sets up Γ = SL₂(Z[ω]), the level-3 subgroup and the cubic Kubota character, and proves a spectral large sieve for cubic metaplectic
  cusp forms (step 5).
- **Which cover.** For plain products the dual coefficients are a_n·g̃6(n), and on squarefree n, g̃6(n) = μ(n)·ε(n mod 4)(n/|n|)conj(g̃3(n)) (§77). In cubic
  terms they carry μ; in sextic terms they are clean. So input 2 lives naturally on the 6-fold cover (sextic Kubota multiplier). Proskurin-type
  general-n Kuznetsov formulas cover it in principle. Its residual spectrum is the sextic theta function, whose Fourier coefficients are not known
  explicitly. Step 4 needs only size bounds there: §83 estimates the Θ contribution at U^{z−1/6} relative.
- **Spectral gap.** I found no explicit spectral gap for cubic or sextic metaplectic Maass forms on SL₂(Z[ω]) beyond the residual theta eigenvalue. The power
  saving in step 5 needs a gap (any fixed one), uniformly in the level (mm′)⁴.
