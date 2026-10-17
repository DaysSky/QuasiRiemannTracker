# The four inputs that would move the quasi-Riemann frontier substantially

DaysSky, with Claude assistance. 17 Oct 2026. Frontier: θ = 0.874836 (Gao, GTimG/QuasiRiemannTracker PR #62), built on the OpenAI
sextic architecture over Q(ω). Details and proofs of every reduction: DH_SEXTIC.md §26–§92. Certified conditional bounds: CONDITIONAL_THEOREMS.md. Research plans: PLAN_RD.md (inputs 1 and 3), PLAN_INPUT2.md (input 2). Labels as in DH_NOTE.md.

Every unconditional lever left inside the architecture is worth at most about 10⁻⁷ (§74). Every substantial lever reduces to one of
the following four statements. The values are measured on the record pipeline (MODEL), or the Cycle61 grid model where marked.

## 1. LM: the long Möbius mean square past the conductor
**Statement.** For s ∈ [1.10, 1.25] and some η > 0: Σ*_{N(u)≍U} |Σ_n μ(n)(u/n)₆W(N(n)/U^s)|² ≪ U^{1 + (5/6)(s−1) − η + ε}·U^s.
- **Value.** dθ/dη ≈ −0.044. Certified in exact rationals: η = 0.001 ⟹ θ ≤ 0.87480 (§73). In slope form (Σ_u|M_u(U^s)|² ≪ U^{1+α(s−1)+ε}·U^s on [1, 3/2],
  copies give α = 5/6): α = 0.80 / 0.75 / 0.70 gives 0.874589 / 0.874186 / ≈ 0.87366 (§89). **Certified in exact rationals (§90): α = 3/4 ⟹ θ ≤ 0.87425,
  α = 7/10 ⟹ θ ≤ 0.87368.** The full LM gives 13/15 (DH_THEOREM.md).
- **Equivalent dual form (PROVED, §69).** With K = D²/U and c_n = ε(n mod 4)(n/|n|)conj(g̃₃(n)) (ε periodic, §77):
  Σ_{N(h)≤K}|Σ_n c_n χ̄_n(h)W|² = K‖c‖²·(local factor) + O(U^ε D²).
- **Why standard tools stop.**
  - Self-duality makes sixth-power copies a fixed point (§3, §70, §75).
  - The row twist (h/n)₂ blocks Rankin–Selberg (§69, §77).
  - Convexity in the dual leaves small-conductor Möbius sums needing D^{0.773} (§79).
  - Large values fail (§9, §71).
  - The thick-family analogue (power-saving Barban–Davenport–Halberstam/Montgomery–Hooley below √x) is open (§77).
- **Evidence (NUMERICAL).** The dual sieve sits on the exact coprimality plateau 0.94 up to K = 32D (§76). The primal plateau holds up
  to U = 10⁶ (§41, §57).

## 2. Plain slot cost c_P: from 6κ toward 1 (critical length)
**Statement.** Σ_u |S_u(Z^{n₁})S_u(Z^{n₂})Q_u|² ≪ Z^{1+ε} for n₁ + n₂ + c·z ≤ 1 with c < 6κ (Lemma 18.1 has c = 6κ; Cycle61
Theorem 2.2 has c = Λ ≈ 3.96 on a window).
- **Value (record pipeline with Gao's constants, §83 addendum 4, re-checked with all conditions in §89).** Slot weight 3 / 2 / 1 gives θ = 0.874313 / 0.873382 / 0.870524, i.e. −5.2·10⁻⁴ /
  −1.45·10⁻³ / −4.3·10⁻³ against the frontier.
  **Certified in exact rationals (§90): plain weight 2 ⟹ θ ≤ 0.8734.**
- **Mechanism (PROVED reading, §80).** The induction amplifies dual rows by h ↦ hp⁶: gain ℓ, row cost 6ℓ, so children carry an excess
  6(w + ℓ) that is paid by deleting slots at κ per unit. Only p⁶ fixes every Gauss sum. Other shifts bring in column phases χ_u(p)^k,
  whose average is a short prime character sum (GRH-type).
- **Sharpest form (§83).** It suffices that the twisted fourth moments T(m, m′) = Σ_u|S₁S₂(u)|²ψ_u(m)ψ̄_u(m′) save N(mm′)^{−θ′} with θ′ > 0.32 (not
  1/2). The fixed-family rows are harmless for slot length z < 1/6. Poisson plus the sieve gives θ′ = 0, because the twisted dual is past the
  conductor by exactly z.
- **Final form (§83 addenda).** Cost 2 is equivalent to the plain fourth moment over the subfamily u ≡ 0 mod m² at its natural size (U/Z²) rather than its
  conductor size (U/Z, which is what Lemma 18.1's twist rule gives). This is a thin-family Lindelöf-on-average with thinness only Z = U^{0.04–0.08}. The
  matching family is quadratic characters, not all Dirichlet characters. There even the untwisted fourth moment has only a log-saving asymptotic
  (Shen–Stucky, arXiv 2402.01497), so a power-saving twisted version is open in the simplest analogue too.
- **Self-duality (§84).** For m ≠ m′, T has no main term in either representation (primal diagonal 0; dual row-diagonal O(Z) by orthogonality), and Poisson
  returns T with coefficient 1. The data put |T| at row-noise level (θ_eff ≥ 0.6, §83 addendum).
- **Averaging over slot primes (PROVED reduction, §92).** Cauchy–Schwarz over the slot primes inside the dual row sum, then positivity on the union of the sublattices m⁴O
  (which overlap boundedly), gives cost 5 with no zero-free input. That strictly improves Lemma 18.1's κ = 1 case (6), but not 6κ₀ ≈ 4.5. It cannot be combined with deletion or
  with live slots to give 5κ, and 5κ₀ would be worth only 0.874661 anyway, because the inverse leg takes over.
- **Closest classical analogue (CITED, 17 Oct search).** Twisted fourth moments of Dirichlet L-functions: Blomer–Fouvry–Kowalski–Michel–Milićević (arXiv 1411.4467, Amer. J. Math. 2017) prove the fourth moment for prime q with error q^{−1/32+ε},
  and Zacharias evaluates a twisted fourth moment. Input 2 needs twists by moduli U^{0.08–0.16} with a saving θ′ > 0.32 per unit of twist conductor. The twist ranges in those papers were not checked here.
- **Evidence (NUMERICAL, §31).** Natural cost 1 holds to within U^{10⁻³} from U = 10⁵ to 3·10⁶, far outside the proved ranges.

## 3. Inverse slot cost c_I: from 2 toward 1 (critical length)
**Statement.** The marked Möbius moment Σ_u|M_u(U^r)Q_u|² ≪ U^{1+ε}·(norms) for r + c_I·z ≤ 1 with c_I < 2.
- **Value.** dθ/dc_I ≈ 2.2·10⁻³ near 2. c_I = 39/20 ⟹ θ ≤ 0.87474 is certified (§73). On the record pipeline (§89, window-aware, slot supply ℓ/h enforced):
  c_I = 1.8 / 1.5 / 1.2 gives 0.874353 / 0.873204 / 0.872762. Below ≈ 1.2 the slot supply ℓ/h ≈ 0.21 binds and the value saturates at 0.872762.
  **Certified in exact rationals (§90): c_I = 3/2 ⟹ θ ≤ 0.87325; c_I = 6/5 ⟹ θ ≤ 0.8728**, and hence any c_I ≤ 6/5 gives the same through this proof.
- **Mechanism (PROVED reading, §32, §65, §68).** The mark q | A enters the cubic-theta reflection conductor through twisted
  multiplicativity (τ(qA′) carries a cubic character mod q). The reflected length is multiplied by q² (+2z), while Lemma 14.2's sieve over
  marks saves at most z. Both sieves in (14.7) are at their diagonals.

- **Barrier form and classical analogue (PROVED reformulation, §88; analogy CITED).** By §88 the marked moment is Σ_{p,p′}a_pā_{p′}T(p, p′), where
  T(p, p′) = Σ_u|M_u(N)|²χ_{pp′}(u) is the Möbius mean square twisted by a character of the row variable modulo q = pp′. The large sieve on the progressions u mod q
  (size U/q, length N) is sharp only for q ≤ U/N, which is exactly r + 2z ≤ 1. Cost 1 needs q up to (U/N)², i.e. beating the large-sieve barrier for moduli
  q = pp′ that are well-factorable products of slot primes with product-form weights a_pā_{p′}. That is the shape of Bombieri–Friedlander–Iwaniec (level x^{4/7}
  for well-factorable weights, via dispersion and Deshouillers–Iwaniec Kloosterman bounds), where the classical barrier is √x. The sextic version would need
  metaplectic Kloosterman/Kuznetsov input on the cubic cover, the same spectral programme as PLAN_INPUT2.md.
  **BFI strength, certified (§93):** a level-exponent gain of 8/7 over this barrier (the BFI factor) is c_I = 7/4, which gives θ ≤ 0.87425 (window-aware top row, all side conditions pass). **Maynard strength** (x^{3/5}, arXiv 2006.07088; factor 6/5) is c_I = 5/3, which gives **θ ≤ 0.8740**, or 0.8738 together with LM slope 4/5. **Lichtman strength** (x^{66/107}, arXiv 2309.08522) is c_I = 107/66, which gives **θ ≤ 0.87382**, or 0.87362 with slope 4/5.
- **Dual Type II form (PROVED reformulation).** After Poisson in u and the shift identity, the marked dual row is R(j) = Σ_{p,n}a_pb_n·g̃₃(pn)·χ̄_{pn}(j):
  a Type II sum of twisted cubic Gauss sums with the variables coupled by (n/p)₃, in mean square over j ≤ J = (NZ)²/U. Cost 1 is the bilinear large sieve
  (J + NZ)·‖a‖²‖b‖² for it. Heath-Brown–Patterson-type bilinear estimates for cubic Gauss sums handle such couplings through Heath-Brown's cubic sieve. The paper's
  Lemma 14.2 already uses that sieve for this coupling, and it is at its diagonal in the binding regime (§68). Going past it needs input on the Gauss-sum side (spectral).

## 4. Sixth moment of the HIGH-side numerator
**Statement.** Σ_u|L(1/2 + ϑ + it, ψ̄*_u)|⁶ ≪ U^{1+λ} with λ < (1 − R)/2 ≈ 0.20, or equivalently for the central approximate-functional-equation
block of length U^{1/2}.
- **Value (corrected §91).** Ignoring slots it improves the HIGH exponent by d(1 − R − 2λ)/12, up to −7·10⁻³ at λ = 0 (§81). The fourth route keeps a slot gain ϑ/2
  that a plain sixth route lacks, so the requirement is λ_m < (1 − R)/2 − 3ϑ(1 − 2m) on every block m ≤ 1/2, and the full gain is roughly −2·10⁻³ (MODEL, rough).
- **Why it stops.** The central block has no slot room. Three plain factors of length U^{1/2} cannot be shortened by reflection (total 3/2,
  supercritical), so copies give λ = 5/12. The fourth moment plus a Burgess cap gives λ = 3/8 (§82 addendum). The moment Σ_u Π_{i≤6} L(s_i, ψ_u^{±1}) has an indefinite Weyl group, so there is no
  multiple-Dirichlet-series continuation (§68).
- **Against the ζ analogue (CITED, §87).** The best unconditional sixth moment of ζ, ∫|ζ(1/2+it)|⁶ ≪ T^{5/4+ε}, has loss 1/4. λ < 0.20 is beyond it;
  even a twelfth moment ≪ U^{2+ε} for the family would give λ = 1/4 by the same Hölder step.
- **Evidence (NUMERICAL, §82).** With the Θ row u = 1 removed, m₆/m₂³ of the central block equals a random multiplicative control to within 1.5% up
  to U = 3·10⁶ (1.09·10⁷ rows). No power loss is visible.

## Literature check, 17 Oct 02:45Z (CITED, abstracts read as data)
- Diaconu–Ion–Paşol–Popa (arXiv 2607.27131, Jul 2026): first and second twisted moments of r-th order Hecke L-functions over fields containing μ_{2r},
  r ≥ 3. For the sextic family this needs μ₁₂ ⊄ Q(ω); it covers the cubic family over Q(ω). Second moments only, so it does not reach the fourth-moment
  shape of input 2.
- de Faveri–Dunn–Hoffstein (arXiv 2607.07911): the cubic and quartic large sieves are not perfectly orthogonal (Gauss-sum bias). They conjecture the same
  for every order n ≥ 3. So the general-coefficient sextic large sieve cannot have the (U + N) shape: inputs 1–3 must use the arithmetic of μ and of the plain
  coefficients.
- McGoldrick (arXiv 2610.17144) and Sang (arXiv 2610.15291, θ = 0.8748867, superseded by the record) supply none of the four inputs.

## The single statement behind inputs 1 and 3 (PROVED reductions §26, §69, §88; certified values §90)
**(RD_{c₀}) Rows-dominant asymptotic for the canonical family.** For Lemma 17.2's family (rows k of norm ≤ Z^M, squarefree columns with the bias-free
coefficients ᾱγ₂(n)·d(n)), in the rows-dominant range F₀ < M ≤ (1 + c₀)F₀:
  E = (diagonal main term)·(1 + O(Z^{−η})),
i.e. the off-diagonal is smaller than the diagonal by a power. This is the regime where only the diagonal bound Z^M is proved.
- Without marks it gives the LM slope 5/(6(1 + c₀)) (§26). c₀ = 0.042 beats the record. c₀ = 1/9 is slope 3/4, which ⟹ θ ≤ 0.87425 (certified). c₀ = 0.19 is slope 0.70,
  which ⟹ θ ≤ 0.87368 (certified).
- With product-form marks of length z₀, the same estimate for M up to F₀ − (c_I − 1)z₀ (the proved range is M ≤ F₀ − z₀; §26, §32) gives inverse cost c_I. c_I = 6/5 ⟹ θ ≤ 0.8728 (certified).
- Its dual form is the ε-free sextic large sieve with exact leading term for c_n = ε(n mod 4)(n/|n|)conj(g̃₃(n)) (§69). Numerically it sits on the coprimality plateau (§76).
- It is an asymptotic, not an upper bound. In the dual form the main term K‖c‖² must be exact and the error O(U^εD²), a relative saving U/D. So Lindelöf-on-average bounds
  for Kubota's Dirichlet series (as used by de Faveri–Dunn–Hoffstein, 2607.07911) do not supply it.
- What it needs beyond the paper: Möbius-weighted cubic Gauss sums with no bias at the copy rows k = k₀a⁶, i.e. cancellation in Σ_n μ(n)g̃₆(n)χ̄_n(k₀)
  at small conductors k₀ (§79: pointwise this needs D^{0.773}). On average it is the analogue of breaking the large-sieve barrier for well-factorable moduli
  (Bombieri–Friedlander–Iwaniec), with metaplectic Kloosterman/Kuznetsov input on the cubic cover in place of Deshouillers–Iwaniec.

## Common roots
- **Copy efficiency 1/6.** Sixth powers are the only multipliers fixing every column character. This is forced by the cubic theta series:
  the Poisson frequencies of the averaged base sum are u·a⁶ (OpenAI overview; §78). It sets α′ = 5/6, c_P = 6κ, Theorem 2.2's K_a, and through
  ϑ ≤ q/c_P part of the HIGH side.
- **Theta reflection squares a mark's conductor.** This sets c_I = 2.
- **Supercritical moments.** LM and the sixth moment are both beyond the conductor; for LM, self-duality makes copies a fixed point.
