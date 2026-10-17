# Four statements that would move the quasi-Riemann frontier (draft for discussion; not posted)

DaysSky, with Claude assistance.

The current all-moduli Dirichlet zero-free boundary on the QRH tracker is θ = 0.874836 (Gao, PR #62). It is built on OpenAI's sextic architecture over
Q(√−3) and its refinements (Collet, Levy, timbalud, DaysSky, Gao). A systematic sweep found no unconditional input that moves θ by more than about
10⁻⁷. Every substantial improvement needs one of the following four statements. Each comes with its value on the record pipeline and the
mechanism in the existing proofs that it would replace.

Notation: ψ_u(n) = (u/n)₆ on Z[ω], u sixth-power-free with N(u) ≍ U. S_u(X) and M_u(X) are the centrally normalised plain and Möbius sums of length X.

**1. Möbius mean square just past the conductor.** For s ∈ [1.10, 1.25] and some η > 0:
Σ_u |M_u(U^s)|² ≪ U^{1+(5/6)(s−1)−η+ε}·U^s. Sixth-power copies give η = 0.
- Value: θ drops by about 0.044η. η = 0.001 already gives θ ≤ 0.87480, certified in exact rationals. In slope form (copies give 5/6), slope 0.75
  gives θ ≈ 0.87419. The full statement gives 13/15.
- Equivalent dual form: an ε-free sextic large sieve with exact leading term, rows-dominant, for the bias-free sequence ε(n mod 4)(n/|n|)conj(g̃₃(n)).
- Analogue: a power-saving Barban–Davenport–Halberstam variance below √x.

**2. Twisted plain fourth moments (plain slot cost).** For slot primes m ≠ m′ of norm ≤ U^z (z ≈ 0.04–0.08) and n₁ + n₂ + 2z ≤ 1:
Σ_u |S_u(U^{n₁})S_u(U^{n₂})|² ψ_u(m)ψ̄_u(m′) ≪ U^{1+ε}·N(mm′)^{−θ′}.
- Value: θ′ > 0.32 beats the current slot weight 6κ. θ′ = 1/2 gives θ ≈ 0.87338 (−1.45·10⁻³).
- Equivalent form: the plain fourth moment over the subfamilies u ≡ 0 mod m² at natural rather than conductor size, i.e. thin by only U^z.
- Analogue: a power-saving twisted fourth moment for quadratic characters. Even the untwisted one has only a log saving (Shen–Stucky 2024).

**3. Inverse slot cost below 2.** The marked Möbius moment Σ_u |M_u(U^r)Q_u|² ≪ U^{1+ε} for r + c·z ≤ 1 with c < 2 (Lemma 17.1 has c = 2).
- Value: c = 1.95 ⟹ θ ≤ 0.87474 (certified); c = 1.5 / 1.2 gives θ ≈ 0.87320 / 0.87276 on the record pipeline. Below ≈ 1.2 the slot supply ℓ/h binds.
- Mechanism to beat: cubic-theta reflection squares a mark's conductor.

**4. Sixth moments of the long AFE blocks.** Σ_u |S_u(U^m)|⁶ ≪ U^{1+λ_m} with λ_m < (1 − R)/2 − 3ϑ(1 − 2m) for 1/6 < m ≤ 1/2, i.e. below about 0.205 at the central block (copies give 5/12; the fourth moment with a Burgess cap gives 3/8).
- Value: about −2·10⁻³ through the HIGH-side numerator. That needs small loss on every approximate-functional-equation block of length U^m with 1/6 < m ≤ 1/2; the fourth-moment route
  keeps a slot gain the sixth-moment route lacks. For comparison, the best unconditional sixth moment of ζ has loss 1/4.

**Certified conditional values (exact rationals, relative to the radical record's proof; CONDITIONAL_THEOREMS.md).**
- LM slope 4/5 / 3/4 / 7/10 ⟹ θ ≤ 0.87462 / 0.87425 / 0.87368.
- Inverse cost 3/2 / ≤ 6/5 ⟹ 0.87325 / 0.8728.
- Plain weight 2 ⟹ 0.8734.
- Inverse cost 6/5, slope 7/10 and plain weight 2 together ⟹ 0.8708; at natural strength (slope 1/2, plain weight 1) ⟹ 0.8692.
- DH for the family ⟹ 13/15.
Statements 1 and 3 both follow from one rows-dominant asymptotic for the canonical family of OpenAI Lemma 17.2 (PLAN_RD.md).

**Evidence.** All four hold numerically with no visible loss. LM and its dual sieve sit on the exact coprimality plateau up to U = 10⁶. Slot costs are natural
to within U^{10⁻³}. The twisted moments sit at the row-noise level. The sixth moment matches a random multiplicative model up to U = 3·10⁶.

Details, proofs of the reductions, scripts and certificates: available on request (working notes DH_SEXTIC.md §1–§89, FOUR_INPUTS.md, PLAN_INPUT2.md).
