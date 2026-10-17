# The quartic Lindelöf-on-average problem (L³): status report

DaysSky with Claude assistance, 2026-10-15. Labels: **PROVED** (argument given or standard), **NUMERICAL** (checked by computation, no proof), **HEURISTIC**, **OPEN**.

**Summary.**

> **16 Oct, 09:40Z: (L³) contains a quasi-Riemann hypothesis (l3j/THEOREM_ZF.md).**
> * **Statement.** Suppose the squarefree (G)-form of (L³) holds with loss Q^{θ′} for some θ′ < ½, i.e. any power saving below Q^{3/2},
>   uniform in index twists h and over all scales. Passing from ψ_α over all c to this form needs a local inversion that is not done here. Then every quartic Hecke L-function L(s,(k₀/·)₄λ^e) over Q(i) is zero-free in
>   Re s > (7+10θ′)/(8+8θ′).
>   * (L³) itself gives 7/8.
>   * Q^{3/2−δ} gives a strip of width about δ/6.
> * **Mechanism.** The repeated dual frequencies k₀ϖ⁴ carry the Möbius sums of each target character. Positivity comes from a
>   narrow-annulus expansion of the dual weight. The row scale is optimised.
> * **Checks.** The mechanism is confirmed numerically on the μT control family (l3j/j4pred.py). The sign identity is checked by brute
>   force (signcheck.py).
> * **Consequences.**
>   * The Q^{3/2} wall is a quasi-RH wall.
>   * The quartic route is OpenAI's argument with (L³) in the role of its free cubic-theta direct estimate (THEOREM_ZF, last section).
>   * This corrects NOTES §8G and §5 below: the earlier "(L³) ⟹ 7/8" was an outline whose positivity step failed for wide c-ranges.
>   * It also corrects BELOW_3_2 §7: the (2,6,6) model treats Poisson as an FE of one function, which is wrong (l3j/SESSION2.md §9).

* **Status.** (L³) is not proved and not disproved.
* **Numerics.** It holds with constant 1.00 ± 0.01 over every scale tested, up to Z = 6·10⁶ (§8).
* **Conditional proof.** (L³) follows from an explicit ratios-type hypothesis for the quartic Hecke family (§10).
* **Why it is hard.** Unconditionally, every Poisson, reflection or large-sieve route stops at the self-dual point M = 2/3 (§4, §6). The problem is the metaplectic twin of an open GL3 thin-family second moment (§5).
* **Not enough on its own.** Even a proof of (L³) would not reach 5/6 alone: quartic prime-slot moments carry quartic-theta main terms (§7).
* **Same wall for the current endpoint.** The proven sextic endpoint's best remaining lever is the same kind of self-dual barrier (§11).
* **New this round (§10b).** On GRH the T-family obeys the optimal large-sieve shape, S ≪ (QZ + Z²)Z^ε. That gives (G) for all M ≥ 1 and the (L³) sum ≪ Q^{3/2+ε}, against Q^{7/4} unconditionally. The remaining gap to (L³) is exactly Q^{1/2}, and GRH cannot close it by this argument.
* **Confirmed this round (§10c).** The multiple-Dirichlet-series attack on the second moment has an infinite, hyperbolic group of functional equations, so the finite-group continuation mechanism is unavailable. This re-derives an earlier subagent finding (NOTES.md §8).
* **New this round (§10d).** The (L³) rows carry no m-dependent size. In a control, an Euler product of the same conductor shows Selberg-CLT correlations (0.47); the T-family shows none (|corr| ≤ 0.04). So ψ_α hides no Euler product, and GRH methods that act through Euler products cannot apply.
> **Status at 23:45Z, 15 Oct (newest results; details in §13–§15):**
> * (L³) sum Σ_{N(m)≤Q}|ψ_α(m³,½)|²:
>   * **Q^{13/8+ε} PROVED unconditionally** (S49). This is the large-sieve limit.
>   * **Q^{3/2+ε} CONDITIONAL** (S50 + S51; refined at 00:08Z: Q^{1.51}–Q^{1.54} on the CG extension alone, Q^{3/2} needs one more L-factor; see §20 note below). It needs Chinta–Gunnells' G₂ multiple Dirichlet series (n = 4, symbol ε³) extended to an
>     archimedean ℓ = 1 twist and to Q(i) ⊉ μ₈, plus standard convexity bounds.
>   * **Q^{1+ε} OPEN**, and structurally hyperbolic: Coxeter group (2,6,6).
> * Why: Cauchy–Schwarz trades the hyperbolic second moment for the finite G₂ first moment at a cost of Q^{1/2}. The 13/8 → 3/2 step is
>   exactly the T–theta correlation 𝒞(h₀) = Σ_c T(c)τ₄(h₀c), whose residue series has its functional-equation centre at 3/8.
> * Unbalanced quartic large sieve (S49, Theorem 1): proves de Faveri–Dunn–Hoffstein Conjecture 1.3 for max(A,B) ≥ min(A,B)^{6/5}.
> * θ: no change from any of this.
> * **Update 23:58Z (S51 §17–§18).** The conditional step is now narrower.
>   * DDHL (4.28), the Kubota functional equation over Q(i) with n = 4 and ∞-type ℓ, confirms the twist rule. It also fixes the
>     W(G₂)-invariant Gamma normaliser: root twists |ℓ| = 0, 1, 1, 2, 3, 3, and only the d-root is untwisted, so the only pole is the theta pole.
>   * The ℓ = 1 extension is therefore bookkeeping already supported by the literature.
>   * The μ₈ ⊄ Q(i) issue reduces to the quadratic sign (c,−1)_S = χ_c(−1). That sign is a class character mod 4 and preserves CG's
>     space M(Ω) (PROVED).
>   * What remains open is one finite computation: the braid relation for CG's action with that sign kept, on Q₂(i)^×/4th powers (order 256).
>   * **00:03Z (§19).** Exact computer check: CG's local G₂ action with raw Q(i) Gauss sums fails, including σᵢ² ≠ 1, at primes with
>     N ϖ ≡ 5 mod 8. With the Brylinski–Deligne sign χ_ϖ(−1)^{C(i,2)} it satisfies every relation. Globally the sign is χ_c(−1) on the
>     c-node, a class character absorbed into Ψ. Left before Q^{3/2} is unconditional: the global FE write-up, the comparison of
>     non-coprime p-parts, and Bochner plus Phragmén–Lindelöf.
>   * **00:08Z (§20), correction.** The p-part comparison shows that CG's normalising L-factors put L(·, α⁴) and L(·, α¹²) on their
>     critical line at the residue centre.
>     * The α⁴ factor is purely in c and is removed by convolution (PROVED bookkeeping).
>     * The α¹² factor (long root (3,1)) is not. So the CG-conditional result is (L³) ≪ Q^{1.51}–Q^{1.54}, and Q^{3/2} needs that factor handled
>       (for example GRH for one Hecke L-function).
>     * All variants are ≤ Q^{1.54}, below the proved 13/8.
>   * **01:00Z — S52 (l3f/S52_MEANVALUE_BOCHNER.md): an unconditional route below 13/8, with no Chinta–Gunnells and no GRH.**
>     * Method: Bochner's tube theorem on mean-value regions (large-sieve LoA in d, S49 (L³)-type mean squares in c, blockwise
>       convexity) over all 12 G₂ chambers, plus an envelope-of-holomorphy sup argument for polynomial growth and h₀-exponents.
>     * Result: P4-type bounds for X(h₀) for small twists; through S50 and one bootstrap, (L³) ≪ Q^{1.5625}, or Q^{1.5573} with fix (a).
>       These are revised (S52 §10) after correctly accounting for the j⁴-coprimality in S50; the first pass gave 1.5417 / 1.5365.
>     * **01:28Z (S52 §11–§12).** Block-optimal mean-value regions, RMS twists (vector-valued Bochner with the true quartic large sieve in
>       h₀), and an S(q) recursion that feeds S52's bounds back into the regions. Fixed point first reported as Q^{1.5156}; corrected in S52 §14
>       (generic twist range at q > 2/3). Superseded: the theta polar term was missing (S52 §16); the valid value is Q^{1.5573}, and the block-region recomputation is running.
>     * Full sum (S52 §8–§9): reflect each non-squarefree row box through the c-FE (fix (a)); boxes then sit strictly below the main box,
>       and S52 is continuous in (Q,Y). So the full sum is ≪ Q^{1.5573+ε} (after §10), modulo STANDARD bookkeeping (non-coprime corrections, cusp
>       vectors, T↔T* transfer).
>       Large twists Z^{0.24} ≤ N(h₀) ≤ Z^{1/3} are why it stops above 3/2; the Cauchy route cannot go below 3/2.
>     * **02:19Z (S52 §17), supersedes the line above.** Audit fixes: proved bootstrap start, no table clamping, step-up reads; the polar term
>       was already included.
>       * Corrected chain: 1.5104 → 1.5055 at grid step 0.01 (valid); recur10 with linear reads gives 1.5002.
>       * **S52 §18 / S52_LIMIT_PROOF.md: Q^{3/2+ε} given the framework.** It uses an explicit hull point, an exact cone lemma and
>         closed-form bookkeeping on the q-continuum.
>       * **CORRECTION (S52_CORRECTION_NAIVE_VS_CG.md): conditional on task #26 (CG G₂ WMDS over Q(i)).** The binding hull point v₂ = (⅓, −1)
>         lies far outside E′'s region (worst exponent +1, |Bx₂| > 1), so the naive series cannot reach it. The argument holds for the CG
>         series, with transfer back at (−1/12, ¼), where E′ converges absolutely. Unconditional: Q^{13/8} (S49).
>       * **04:35Z, l3g/S53_QI_LOCAL.md: #26 narrowed.**
>         * The local theory over Q(i) is PROVED by exact computer algebra. The correct action is the "slice-signed" one: raw Gauss sums,
>           with Q multiplied by s^{μ(B_i+1)}.
>         * S51 §19's BD sign was the wrong choice. The exponents are unaffected.
>         * 04:55Z: the Ψ-bookkeeping of CG Prop 5.6/5.7 over Q(i) is PROVED, given Kubota's rank-1 FE in Brubaker–Bump form.
>         * #26 is now CONDITIONAL only on that rank-1 FE for Q(i), n = 4 (BB06 assume μ₈; Dia04 / DDHL (4.28) give it in cusp-vector
>           form).
>         * 05:30Z: that rank-1 input is no longer needed in Brubaker–Bump form.
>           * The equivariance follows from DDHL (4.28) + Remark 4.2 by a consistency argument: the ϖ⁴ slice identity, checked by hand,
>             plus injectivity via g(1,π)⁴ = π³π̄ and Hecke equidistribution.
>           * #26: PROVED modulo published DDHL/Dia04 and STANDARD Bochner.
>           * (L³) ≪ Q^{3/2+ε} modulo S52 Lemma 4 constants. Q^{1+ε} still OPEN.
>       * **03:09Z, l3f/S52_FULLSUM.md: the full sum (all four box types) at Q^{3/2+ε}, PROVED given the framework.**
>         * Boxes reduce to twisted main-box problems. The twist conductor F costs 7/12 at the hull point (EXACT F-cocycle) plus
>           ≤ 5/6 + κ/2 in the cone lemma (by reduction to the untwisted one).
>         * The box sum converges for every twist cost c < 2.75; the cost obtained is ≤ 1.57.
>         * Remaining open items: Lemma 4 constants and the T↔T* transfer.
>       * **Below Q^{3/2} (l3f/BELOW_3_2.md).** The problem is self-dual on the range Z < N(k) ≤ Z^{4/3}. Every input used so far is invariant
>         under that duality, so (L³) needs a new one (candidates: metaplectic Kuznetsov, a ratios-type Möbius variance, or μ-Type II
>         structure).
>       * Mechanism: the Region_d boundary point (⅓, 0) in two non-polar chambers, closed by the (L³)-shaped block D = C^{2/3}; contraction
>         s − 2 ↦ ≈0.16(s − 2).
>       * Q^{3/2} is the Cauchy floor and equals the GRH value of §10b. (L³) (Q^{1+ε}) stays OPEN, and θ is unchanged.
>     * §21: a leftover p-part term (5,2) also blocks the centre (X₁ > −1/10). So even GRH for L(α¹²) gives Q^{1.506}; the exact
>       Q^{3/2} needs a matched family.


## 1. Why it matters

The OpenAI quasi-Riemann proof reaches 7/8 using a sextic family over Q(√−3). The current tracker refinements reach 0.8748389.

* **Exponent bookkeeping.** For a family of order k the architecture gives σ₀ = 1 − 1/(2k) − 3ℓ/(2k), where ℓ is the compensation length.
* **The quartic family over Q(i).** For k = 4 this is 7/8 − 3ℓ/8, equivalently σ₀ = 1 − M/4 where M is the row-range exponent of (G). Through the same architecture, (G) at M = 2/3 would give **θ = 5/6**; the quartic ceiling is 4/5.
* **The cost.** It needs one new input, (G) below. It also needs quartic analogues of the other sextic moment lemmas (§7).
* **The sextic ceiling.** The sextic architecture cannot go below 13/15 even with ideal inputs.

## 2. The statement

Notation:

* K = Q(i), χ_c = (·/c)₄ for primary squarefree c, and α(c) = c/|c|.
* g̃(χ) = N(c)^{−1/2} Σ_x χ(x) e(Im(x/c)), and T(c) = α(c) g̃(χ̄_c).
* W is a smooth bump on [1/2, 1].

**(G) at M = 2/3** (OPEN). For Q = Z^{2/3}:

  S(Q,Z) = Σ_{N(m)≤Q} |R(m)|² ≪ Q Z^{1+ε},  where  R(m) = Σ_c T(c) χ̄_c(m) W(N(c)/Z).

**(L³)** is the equivalent Dirichlet-series form:

  Σ_{N(m)≤Q} |ψ_α(m³, ½+it)|² ≪ Q^{1+ε}(1+|t|)^A,  where  ψ_α(r,s) = Σ_c α(c) g̃(χ̄_c, r) N(c)^{−s}.

## 3. Proved structure

1. **T(c)² = μ(c)α(c)** for every squarefree primary c (PROVED in outline: the sign of the normalized quadratic Gauss sum at primes is taken from the mod-4 evaluation; also NUMERICAL, see below).
   * Primes: Jacobi sums give J(χ,χ) = −χ(−1)π, and g(χ)² = J(χ,χ)g(χ²) with the quadratic Gauss sum.
   * Composites: T² is multiplicative, because the twisted-multiplicativity cocycle is ±1.
   * Checked exactly on all 248 squarefree c with N(c) ≤ 700, on 78 350 primes up to norm 10⁶, and (via the multiplicative pipeline) on all 347 518 squarefree c with N(c) ≤ 10⁶, with the sign matching μ(c) on 20 000 random c.
2. **T(c) = μ(c)·conj g̃(χ̄_c)** (PROVED, from 1 and |T| = 1). At primes this reads α(π)g̃(χ̄_π) = −conj g̃(χ̄_π) exactly.
3. **Twisted multiplicativity** (standard, from CRT and quartic reciprocity; NUMERICAL check 118/118 pairs): T(ab) = T(a)T(b)(a/b)₂ s(a,b), with s = −1 iff N(a) ≡ N(b) ≡ 5 (mod 8). Equivalently, g̃(χ̄_{ab}, m) = s·g̃(χ̄_a, m b²)·g̃(χ̄_b, m).
4. **Rows as Möbius-weighted Kubota sums** (PROVED): R(m) = conj Σ_c μ(c) g̃(χ̄_c, m) W(N(c)/Z).
   * The Möbius weight raises the conductor from N(m) to N(m)³. That is why the family is thin: Q members of conductor Q³.
5. **Elliptic-curve form** (PROVED; point counts confirm the normalization, 125/125). α(c)χ̄_c(m) is the normalized Grössencharacter of the quartic twist E_m : y² = x³ − m x. So R(m) = Σ_c g̃(χ̄_c) ψ̃_{E_m}(c) W.
6. **Conjugate pairs** (PROVED): T(c̄) = conj T(c), hence |R(m̄)| = |R(m)|.

## 4. Equivalent forms and the fixed point

* **Poisson in m** (PROVED algebra) maps (G) to the Möbius mean square at H = Z²/Q = Z^{4/3}, D = Z:
  Σ_{N(k)≤H} |Σ_c μ(c)χ_k(c) N(c)^{−1/2} W|² = H·C + O(D).
  It uses T(c)g(χ̄_c,k) = μ(c)χ_c(k)√N(c) and the CRT factorization of g(χ̄_cχ_{c'}, k).
* **Poisson back** returns (G) with the same parameters. At M = 2/3 the composite map is Fourier inversion, so no contraction is available.
* **Kubota reflection** maps (G) at M to (G) at M/(3M − 1). The fixed point is M = 2/3. With an optimal dual large sieve this proves (G) only for M ≤ 1/2; Gao–Zhao's quartic large sieve gives only M ≤ 3/8. (Caveat: Gao–Zhao is stated for squarefree rows. Rows m = m₁m₂²m₃³j⁴ need separate square-class bookkeeping, and the j⁴-multiplicity alone forces Δ ≥ Q^{1/4}Z; see S43 and S47. The ranges quoted here are at the level of exponents.)
* **Parity form.** Let K_z(m) = Σ z^{ω(c)} g̃(χ̄_c, m)W. Then R̄ = 2K₊ − K₁. So (G) says the even-ω part carries exactly half of the parity-blind sum, up to square-root error. This is Type-II information of square-root strength, behind Selberg's parity barrier.

## 5. Where (L³) sits

* **Ratios conjecture ⟹ (G) for all M ∈ (0,1] ⟹ (L³).** Take the Conrey–Farmer–Zirnbauer ratios conjecture for the family {L(s,(k/·)₄)}, with denominators only. Its error H^{1/2+ε} beats the needed O(H^{3/4}).
* **GRH alone does not suffice by any known argument.** The natural GRH argument, pointwise bounds for 1/L, loses H^{1/2} on the off-diagonal, and ψ_α has no Euler product. Whether GRH implies (L³) is not known. What GRH does give is proved in §10b: the optimal large-sieve shape, i.e. (L³) up to a factor Q^{1/2}.
* **Ratio of conductor to family size.** In the Kubota-index picture, ψ(r, ½) has conductor N(r). Over all indices N(r) ≤ X (ratio 1) Lindelöf on average is a large-sieve statement. Over squares r = a² (ratio 2, quadratic characters) David–Dunn–Hamieh–Lin prove A^{3/2} and expect A^{1+ε} "with additional work" (arXiv 2306.11875, (1.11) and Remark 9.7). Over cubes r = m³ (ratio 3, quartic characters) it is (L³). Ratio 3 is beyond any large sieve, like the sixth moment of quadratic Dirichlet L-functions, which is open unconditionally and known on GRH only through Euler products (Soundararajan's method), which ψ_α lacks.
* **(L³) ⟹ consequences.**
  * A zero-free half-plane Re s > 7/8 + θ′/2 for the whole quartic Hecke family over Q(i). This is "Theorem Z" (NOTES §8), re-derived in outline this session.
  * The argument: the k = k₀j⁴ with N(j) ≤ (H/N(k₀))^{1/4} all give the same character χ_{k₀} on c coprime to j. Positivity then gives |Σ_{c~D} μ(c)χ_{k₀}(c)|² ≤ PO(H,D)/#j, i.e. a bound D^{λ/2 + max(1−η/8, ½+3η/8)}, so 1/L(s, χ_{k₀}) is analytic beyond that exponent. Outline only: the coprimality to j changes each finite sum by shorter sums, which need an inductive treatment.
  * The exponent is optimized at η = 1, the square case H = D, which (L³) uniform over scales implies; at η = 4/3 it gives only Re s > 1.
  * It also implies Weyl-strength subconvexity for each ψ_α(m³,·).
* **GL3 analogue** (OPEN). A degree count makes ψ_α a degree-3 object, so (L³) is the metaplectic analogue of Σ_{χ mod q} |L(½, π⊗χ)|² ≪ q^{1+ε} for GL3 π. That is open even for π = 1⊞1⊞1, the sixth moment for a single modulus. The known sixth-moment results (Huxley's upper bound; the Chandee–Li–Matomäki–Radziwiłł asymptotic, arXiv 2409.01457) all average over moduli q ≤ Q.

## 6. Routes tried and why each closes

| Route | Outcome |
|---|---|
| Poisson in m / in k | exact duality, fixed point at M = 2/3 |
| Kubota reflection | self-dual at M = 2/3 |
| AFE for the Möbius side | 1/L has no functional equation; contour shifts cross zeros |
| Bilinear (Type II) via the cocycle | creates index-shifted rows of conductor N(m)³N(b)²; Cauchy–Schwarz loses B |
| Large sieve (any shape) | needs dual length ≤ rows, i.e. M ≤ 1/2; the optimal shape is unknown, and false for general coefficients (de Faveri–Dunn–Hoffstein) |
| All-index Lindelöf on average restricted to cubes | Q³ |
| Cube detection by sieve or larger sieve | Q^{3/2}; or circular (needs a fourth moment) |
| Multiple Dirichlet series | the w-functional equation sends T(c) to μ(c); only (Z/2)² acts, so no continuation past convexity |
| Prime-sign (double-cover) structure | signs carry only the X^{3/4} quartic-theta bias, which μ kills |
| Recursion (Heath-Brown style) | loss λ = 1/6 at H = D^{4/3}; floor 1/12 |

**State of the art for Σ_{N(m)≤Q} |ψ_α(m³,½)|²:**

| Bound | Source |
|---|---|
| Q^{15/8} | Gao–Zhao quartic large sieve |
| Q^{7/4} | Heath-Brown-type recursion (subagent C, λ = 1/6); these methods cannot go below Q^{13/8} (λ ≥ 1/12) |
| **Q^{13/8}** | **PROVED unconditionally (S49, l3c/S49_PROOF.md): sharpened BGL Lemma 3.4 + Poisson in the column variable + Kubota convexity for j⁴-rows. This is the large-sieve limit (dFDH Thm 1.2).** |
| Q^{1.51}–Q^{1.54} (dual-term shift δ ∈ [0, 1/6]) | CONDITIONAL on standard-type facts for one G₂ Weyl-group multiple Dirichlet series (Chinta–Gunnells 2010 G₂, n = 4, symbol ε³, with an archimedean ℓ = 1 twist on the short node; S51 §17–§20) |
| Q^{1.506} | as the row above, plus control of 1/L(½+it, α¹²) (the long-root (3,1) normaliser), e.g. GRH for that one Hecke L-function (S51 §20–§21) |
| Q^{3/2} | as the row above, plus matching the (5,2) p-part (a non-normaliser Euler factor that converges only for X₁ > −1/10; S51 §21) |
| Q^{1.543} | CONDITIONAL on (P4): square-root cancellation in Σ_c T(c)τ₄(h₀c) (T against quartic theta coefficients), S50; numerically supported |
| Q^{3/2} | the optimal large-sieve shape; unproved for quartic symbols in general; PROVED on GRH for T-coefficients (§10b) |
| Q^{1+ε} | target (L³); numerically true with constant ≈ 1 |

## 7. Even a proof of (L³) is not the whole route

The sextic proof's Reflection/, Descent/, CubicSieve/ and Cusp/ layers (about 800 Lean modules in openai-baseline) use the identification "cubic Gauss sums = cubic theta coefficients" (Patterson). The quartic route needs quartic analogues of several lemmas:

* the canonical moment, which is (L³) itself;
* the marked inverse moments with prime slots;
* the plain moments;
* the Gram low side.

The marked inverse moment analogue is again a fixed-point statement of (L³) type.

Prime slots add a further obstruction (NUMERICAL, S32). With prime-only columns the α-twisted mean square is biased and growing: S/diag = 1.23, 1.29, 1.44 and 1.67 at Z = 10⁵, 3·10⁵, 10⁶ and 3·10⁶, while the untwisted prime control stays at 1.00. The excess grows like Z^{1/2}/log Z (predicted 0.71 at 3·10⁶, observed 0.67), and the bias sits entirely in rows with N(m) ≡ 1 mod 8. The reason is that at primes T(π) = −conj g̃(χ̄_π), so prime sums are parity-blind and carry quartic-theta main terms. Quartic prime-slot moments would therefore need explicit control of the quartic theta coefficients, which are unknown at non-square indices; in the sextic setting the analogous terms are tractable because cubic theta coefficients are explicit (Patterson).

The smallest θ gain from the quartic route (HEURISTIC interpolation) needs (G) just inside the open range: M slightly above 1/2, or equivalently, by reflection, slightly below 1 (computed at M ≈ 0.999). Even there it needs two things:

* an optimal-shape quartic large sieve for T-coefficients;
* a further power saving of about Z^{6·10⁻⁴} beyond it.

## 8. Numerics (NUMERICAL; l3/gen_rows.c)

| Z | Q = Z^{2/3} | rows | S/diag (α-twisted) | S/diag (untwisted) | row variance (Exp: 1) | max row |
|---|---|---|---|---|---|---|
| 1e5 | 2154 | 847 | 1.005 | 0.996 | 1.09 | 6.3 |
| 3e5 | 4481 | 1759 | 1.023 | 1.088 | 0.94 | 7.3 |
| 1e6 | 9999 | 3928 | 0.996 | 1.008 | 1.01 | 7.6 |
| 3e6 | 20800 | 8172 | 1.042 | 1.064 | 1.13 | 13.1 |
| 6e6 | 33019 | 12961 | 1.010 | 1.045 | 1.06 | 9.1 |

* **Bias test.** Across 36 dyadic windows (Z from 1.2·10³ to 2.5·10⁶), the row-weighted mean of S/diag is 1.008 ± 0.007, with the s.e. understated by up to √3 because windows overlap. At Z = 2.5·10⁶ it is 0.990. Across the row-exponent scan M ∈ {0.55, 0.6, 2/3, 0.75, 0.85} at Z = 3·10⁵ it is 1.01–1.05; the reflection pair M = 0.6 and 0.75 at Z = 10⁶ gives 1.012 and 1.004, and M = 0.85 at Z = 10⁶ (49 439 rows) gives 0.998. In the square case M = 1 (Q = Z = 3·10⁵, 117 802 rows) it is 1.006, so even the optimal large-sieve shape, unproved for quartic symbols, holds for T-coefficients numerically. So S ~ diag with constant 1.00 ± 0.01.
* **Twisted rows.** |R(m)|²/diag is Exp(1)-distributed, i.e. random-phase behaviour. This is consistent with (L³), and even with a pointwise bound of order log Q. The excess at 3·10⁶ (z ≈ 2.6 with the paired sample size) did not persist: at 6·10⁶ the ratio is 1.010 ± 0.013, with exactly exponential tails and a typical maximum.
* **Untwisted control.** It has theta-bias rows at N(m) = 17³ = 4913, 29², 13², 135² = 18225 and 5⁶ = 15625, i.e. where m³ has a square class or a large fourth-power part.
* **Theta terms in parity-blind sums.** Σ_c g̃(χ̄_c, π)W(N(c)/Z) at a prime index π has a visible Z^{3/4} term (|coef| ≈ 0.025) for N(π) = 17 and 41, which are ≡ 1 mod 8. It is zero within ±0.002 for N(π) = 5, 13, 29, 37 and 53, which are ≡ 5 mod 8. Observed rule: the quartic-theta residue at π vanishes iff χ_π(−1) = −1 (S29).
  * Confirmed on all 160 split primes of norm ≤ 1000, in three non-overlapping groups. N ≡ 5 mod 8 gives at most 0.0042. N ≡ 1 mod 8 with 2 a biquadratic residue gives 0.006–0.010. N ≡ 1 mod 8 with 2 not a biquadratic residue gives 0.015–0.027.
  * Within N ≡ 1 mod 8 the size has two clusters. It is large (≈ 0.02–0.027) when 2 is not a biquadratic residue mod p, and small (≈ 0.009) when p = x² + 64y² (Gauss's criterion for 2 being one).
  * Interpretation: for primary π, N(π) ≡ 5 mod 8 ⟺ π ≡ 3+2i mod 4. Residues at ±iπ are nonzero. So the vanishing is a 2-adic congruence condition on the index, the usual kind of local condition for theta coefficients.
  * These are numerical observations about the quartic theta function. The Möbius-weighted family used for (L³) has no such term: at all 160 prime indices it stays at noise level, max 0.0046 (S35).
* **Prime sums.** Σ_{N(π)≤10⁶} g̃(π) shows the X^{3/4} bias: class 1 mod 4 gives +1838 (z = 9.3), class 3+2i gives −760.
* **Moments and phases.** E[x²] ≈ 2 and E[x³] ≈ 6 (Gaussian values) at every Z. Row phases are uniform. Twisted and untwisted rows of the same m are uncorrelated.
* **Weight.** A sharp cutoff (W = 1 on (Z/2, Z]) gives the same: 1.010 at Z = 3·10⁵ and 0.984 at Z = 10⁶. So does a wider bump on (Z/4, Z]: 0.984 at Z = 10⁶.
* **Reproducibility.** Programs and data are in `l3/` (see `l3/README.md`), and every result is logged in `L3_ATTEMPT.md` (S1–S37). The C pipeline matches brute-force Gauss sums to 1e-12 for every squarefree primary c with N(c) ≤ 15 500, 5383 elements (S17).

## 9. Bottom line

(L³) is neither proved nor disproved. Numerically it holds with constant ≈ 1. Proving it requires square-root, parity-sensitive cancellation for a thin degree-3 metaplectic family. That is at least as hard as Lindelöf on average for GL3 twists over a single modulus, which is open. Every method based on Poisson, reflection, the large sieve or multiple Dirichlet series stops at the self-dual point M = 2/3. What can be proved now is the conditional statement in §10: a ratios-type hypothesis for the quartic Hecke family implies (L³).

## 10. A conditional proof: ratios-type hypothesis ⟹ (G)

This section proves (G), and with it (L³), from one explicit hypothesis about the quartic Hecke family. The hypothesis is a standard-form instance of the Conrey–Farmer–Zirnbauer ratios conjecture. It is **unproved**. Steps marked "standard" are routine analytic number theory whose details are not written out here.

**Family.** For k ∈ Z[i]∖0 and a class character ν of (Z[i]/16)^×, let ξ_{k,ν}(c) = (k/c)₄·ν(c), as a function of odd c. Via quartic reciprocity this is a Hecke character of c with conductor dividing 16·rad(k).

**Hypothesis (R).** Fix ε > 0, a pair of classes (ν, ν′), and a smooth V supported in [1/2, 2]. For some fixed B, uniformly for ε ≤ Re u, Re v ≤ 1/4 and all imaginary parts,

  Σ_{k ∈ Z[i]∖0} φ(k) V(N(k)/H) · [L(½+u, ξ_{k,ν}) · L(½+v, conj ξ_{k,ν′})]^{−1} = H·Ṽ·A_{ν,ν′,φ}(u,v) + O(H^{1/2+ε}(1+|u|+|v|)^B)

for every function φ of k mod 16. (The Poisson dual runs over all nonzero k, not only primary ones; φ absorbs the 2-adic phase e(Im(k/(2+2i))) and the associate classes.)

Here A_{ν,ν′,φ}(u,v) is the "diagonal" Euler product of the recipe: the family average of ξ_{k,ν}(c)·conj ξ_{k,ν′}(c′) for squarefree c, c′ vanishes unless c = c′. For a Möbius-only (denominators-only) ratio the recipe has no swap terms. Cross pairs ν ≠ ν′ are needed (step 3).

**Theorem (conditional).** Assume (R). Then for Q = Z^M with 0 < M ≤ 1,

  S_w(Q,Z) = Σ_m w(N(m)/Q)|R(m)|² = D(Q,Z) + O(Q·Z^{(2−M)/2+ε}),

where D(Q,Z) ≍ QZ is the diagonal. In particular, at M = 2/3 the error is QZ^{2/3+ε}, so S ~ D with a power saving and (G) holds.

**Proof.**

1. *Open the square.* S_w = Σ_{c,c′} T(c)·conj T(c′)·W W′ · Σ_m w(N(m)/Q) χ̄_c(m) χ_{c′}(m).
   * The terms with c = c′ give D(Q,Z).
   * Smooth majorization of a sharp cutoff is harmless because all terms of S are nonnegative.
2. *Poisson in m* (standard, over the lattice coset of primary elements). For c ≠ c′ the character ψ = χ̄_c χ_{c′} is nontrivial mod cc′, so the dual frequency k = 0 vanishes. The dual sum has length H ≍ N(cc′)/Q ≍ Z²/Q = Z^{2−M}. The coset produces finitely many twists ν mod 16.
   * Detail. Primary m form the coset 1 + (2+2i)Z[i], and Z[i] is self-dual for e(Im(·)). Poisson over q(2+2i)Z[i], with q = cc′, and CRT give
     Σ_{m≡1 (2+2i)} ψ(m) f(m) = N(q(2+2i))^{−1} Σ_{k∈Z[i]} ψ̄(u_q) g(ψ,k) e(Im(k/(2+2i))) f̂(k/(q(2+2i))),
     where u_q = (2+2i)^{−1} mod q. The finite character-sum factorization behind this was checked by brute force (48 cases, error 7e-14). So ψ̄(u_q) = χ̄_c(2+2i)χ_{c′}(2+2i), which by the supplementary laws depends only on c, c′ mod 16 (the ν-twists), and the 2-adic phase depends only on k mod (2+2i) (absorbed into V_k).
   * Non-coprime pairs (g = gcd(c,c′) ≠ 1): the character is imprimitive. Its Gauss sums vanish unless g | k, and the contribution is bounded the same way with smaller moduli (standard).
3. *Gauss-sum algebra* (PROVED). By CRT (write x = x₁c′ + x₂c), g(χ̄_c χ_{c′}, k) = χ̄_c(c′) χ_{c′}(c) · g(χ̄_c, k) · g(χ_{c′}, k). This was also checked numerically on 41 coprime pairs, with no extra sign. Combining with T(c)·g(χ̄_c,k) = μ(c)χ_c(k)√N(c) and conj T(c′)·g(χ_{c′},k) = μ(c′)χ_{c′}(−1)χ̄_{c′}(k)√N(c′), the reciprocity factors become ±1 by class mod 16. Both identities were checked on 160 random (c,k), with error below 7e-12, including k not coprime to c.
   * Splitting the sign. The quartic-reciprocity sign s(c,c′) = (−1)^{((Nc−1)/4)((Nc′−1)/4)} is not a product f(c)g(c′). With λ(c) = (−1)^{(Nc−1)/4} it equals ½(1 + λ(c) + λ(c′) − λ(c)λ(c′)), a sum of four products. The other class factors (χ_{c′}(−1), ψ̄(u_q)) are products.
   * Result: a finite sum of cross terms,
     Off_S = Q · Σ_{(ν,ν′)} a_{ν,ν′} Off_M^{(ν,ν′)},  with  Off_M^{(ν,ν′)} = Σ_k V_k Σ_{c≠c′} μ(c)μ(c′) ξ_{k,ν}(c) conj ξ_{k,ν′}(c′) (N(c)N(c′))^{−1/2} W̃ W̃′.
   The weights V_k come from ŵ(kQ/cc′). They separate from (c, c′) after a Mellin transform in N(cc′)/Z², with rapidly decaying transforms (standard).
4. *Mellin in c and c′.* Σ_c μ(c) ξ(c) N(c)^{−1/2−u} W(N(c)/Z) = (2πi)^{−1} ∫ Ŵ(u′) Z^{u′} L(½+u+u′, ξ)^{−1} du′ (up to local factors at 2), and likewise for c′. So the full mean square Σ_k V_k |Σ_c …|² equals a double integral, on Re = ε, of Ŵ(u)Ŵ(v)Z^{u+v} times the family sum in (R).
5. *Apply (R) to each pair (ν, ν′).* The main term H·Ṽ·A_{ν,ν′}(u,v), integrated, reproduces exactly the c = c′ terms: A_{ν,ν′} is the diagonal series. Hence Off_M^{(ν,ν′)} = O(H^{1/2+ε}Z^{2ε}). The integrals converge because Ŵ decays faster than any power of |Im|, which beats the polynomial factor in (R).
6. *Conclude.* Off_S ≪ Q·H^{1/2+ε} = Q·Z^{(2−M)/2+ε}. At M = 2/3 this is QZ^{2/3+ε} = o(QZ). ∎

**Remarks.**

* **Numerics agree.** Numerically S/diag = 1.008 ± 0.007 over 36 windows (§8), consistent with a power-saving off-diagonal. The random-model size of Off_S is also Q·H^{1/2}, the same as the error in (R).
* **What (R) does.** It replaces exactly the missing input: square-root control of the Möbius off-diagonal, uniformly in the shifts. Any proof of (L³) must supply an equivalent statement for this family. By §4 and §6, Poisson, reflection and large-sieve methods return to the same fixed point.
* **Status of (R).** It is not known to follow from GRH. It is predicted by the CFZ recipe, the same heuristic that predicts all known moment asymptotics.

## 10b. What GRH gives: the optimal large-sieve shape (PROVED, conditional on GRH)

**Theorem G1.** Assume GRH for ζ_{Q(i)} and for the Hecke L-functions of the finite-order characters ξ(c) = (r/c)₄ν(c) of Q(i), with r ∈ Z[i]∖0 and ν a class character mod 16. For smooth w, W supported in [½, 2] and [½, 1], and all Q, Z ≥ 2:

  S_w(Q,Z) = Σ_m w(N(m)/Q)|Σ_c T(c)χ̄_c(m)W(N(c)/Z)|² ≪_ε (QZ)^ε (QZ + Z²).

**Proof** (details in L3_ATTEMPT S38).

1. Open the square. Write g = gcd(c,c′), c = gc₁, c′ = gc₁′. The m-character is χ̄_{c₁}χ_{c₁′}·1_{(·,g)=1}, modulo the squarefree L = gc₁c₁′.
2. Poisson over the primary coset, with CRT, gives Gauss sums g(χ̄_{c₁},k)·g(χ_{c₁′},k)·r_g(k), where r_g is the Ramanujan sum, and the weight Ŵ_w(QN(k)/(8N(L))).
3. Apply T(c₁)g(χ̄_{c₁},k) = μ(c₁)χ_{c₁}(k)√N(c₁), its conjugate analogue for c₁′, and twisted multiplicativity in g. Each (g,k) term becomes (Q/8N(g))·r_g(k)·Σ μμ′ ξ(c₁)ξ′(c₁′)N(c₁c₁′)^{−1/2}WW′Ŵ_w, with Hecke characters ξ, ξ′ of the hypothesis type.
4. Mellin-separate Ŵ_w on Re z = 1 + ε. GRH gives Σ_{N(c)≤y, (c,g)=1} μ(c)ξ(c) ≪ y^{1/2}(yN(kg))^ε, so each Möbius sum is ≪ (Z/N(g))^{1+ε}Z^ε.
5. Σ_k N(k)^{−1−ε}|r_g(k)| ≪ τ(g). So the g-part is ≪ τ(g)Z^{2+ε}/N(g). Summing over g gives Z^{2+ε}. The k = 0 term is the diagonal, ≪ QZ. ∎

**Consequences.**
* (G) holds for every M ≥ 1. In particular sub-problem 1 of §12 (the square case) holds on GRH.
* With the approximate functional equation of ψ_α: Σ_{N(m)≤Q}|ψ_α(m³,½+it)|² ≪ Q^{3/2+ε}(1+|t|)^A on GRH.

**Context.**
* This is the quartic analogue, for T-coefficients, of the GRH dispersion estimates of Dunn–Radziwiłł (Annals 200 (2024), Prop. 9.1–9.2). There the extremal sequences for Heath-Brown's cubic large sieve are the Gauss sums themselves, and a theta correction term is needed. Here the Möbius factor in T = μ·conj g̃ removes the theta term, and the shape is optimal.

**Why it stops at Q^{3/2}.** (L³) needs the k-variance Σ_k(|P(k)|² − C) ≪ Z = H^{3/4}. GRH controls each |P(k)| by Z^ε, but not their correlations across k.

## 10c. The multiple-Dirichlet-series group is infinite (PROVED; l3b/mds_group.py)

Take the series Z_T(s₁,s₂,w) = Σ_m N(m)^{−w}ψ_α(m³,s₁)·conj ψ_α(m³,s̄₂), in coordinates x_j = s_j − ½, y = w − ½. Its functional equations act linearly:

| Map | Source | Action |
|---|---|---|
| K₁, K₂ | Kubota FE in s₁ or s₂ (conductor N(m)³) | x_j ↦ −x_j, y ↦ y + 3x_j |
| P | Poisson in m (T-form ↔ Möbius form Σ_k N(k)^{−½−b}/(L(½+a₁,ξ_k)L(½+a₂,ξ̄_k))) | (x₁,x₂,y) ↦ (x₁+y, x₂+y, −y) |
| J | joint FE of the two 1/L factors | (a₁,a₂,b) ↦ (−a₂, −a₁, b − a₁ − a₂) |

* In T-coordinates J becomes J_T = [[−1,−2,−4],[−2,−1,−4],[1,1,3]].
* K₁, K₂ and J_T are reflections whose invariant form has eigenvalues −5/2 and (3 ± √19)/2, so it is indefinite.
* K₁J_T has an eigenvalue −9.899, so it has infinite order. The group grows exponentially: 229 elements of word length ≤ 8.

The Bochner continuation used for cubic moments (a finite Weyl group) is therefore unavailable. This matches the sixth moment of quadratic Dirichlet L-functions, whose group is also hyperbolic. It re-derives, with different generators, an earlier finding recorded in NOTES.md §8 (a hyperbolic Coxeter group (2,6,6), together with the remark that GRH gives exactly Q^{3/2}). §10b is the written-out proof of that remark.

## 10d. No hidden Euler product (NUMERICAL; l3/gen_corr.c, L3_ATTEMPT S41)

The control is an Euler product of the same conductor N(m)³: the partial sums of L(s, ξ_m)³, with coefficients d₃(c)χ̄_c(m). At fixed Q = 9999, over four disjoint dyadic c-windows and 3482 squarefree rows:
* **Control.** It shows Selberg-CLT behaviour: E[x²]/E[x]² ≈ 4, E[x³]/E[x]³ ≈ 30, and cross-window correlation 0.40–0.52.
* **T-family.** It is exactly Gaussian: E[x²]/E[x]² = 2.00, E[x³]/E[x]³ = 6.0, and cross-window correlations −0.04 to +0.02 (s.e. 0.017).

So |ψ_α(m³,½+it)|² has no m-dependent size, unlike any Euler product or Waldspurger-type square root of one. GRH methods that run through Euler products (Soundararajan, Harper) therefore have nothing to act on. A proof of (L³) has to come from cancellation in sums over moduli.

## 10e. Further attempts (round 3) and the precise obstruction

| Attempt | Outcome |
|---|---|
| Zero-density for the quartic family ⟹ optimal shape (S42) | needs A(σ)(1−σ) < 15/8 − 3σ/2. The optimal quartic LS would just give it (A ≈ 2.2 at σ = 0.6); Gao–Zhao gives a trivial density up to σ ≈ 0.7. Closed. |
| Small-conductor rows (S43) | on the Möbius side they need a quasi-RH for quartic Hecke L-functions over Q(i) (implied by (L³) via Theorem Z); on the T-side they are harmless (ψ_α is entire). |
| Finite-type multiple Dirichlet series (S44) | Σ_m g̃(χ̄_m)ψ(m³,s₁)ψ(m³,s₂)N(m)^{−w} is an A3 Weyl group series (finite group), but (L³) is the conjugate-weight-1 version, of hyperbolic type (2,6,6). |
| Phase locking ψ ↔ g̃(χ̄_m) (S45) | would reduce (L³) to the A3 series; numerically absent (all 18 coupling statistics at noise level, 3482 rows). Closed. |

**Function-field contrast (S48).** Over F_q[T] the analogue of (L³) is tame, for two reasons:
* quartic Gauss sums are explicit: g̃(χ_c) = u_{deg c}·conj (c′/c)₄, by the residue theorem and Hasse–Davenport;
* μ(c) = (−1)^{deg c}(c′/c)₂ (Pellet).

So T(c) = v_n·χ_c(−1)(c′/c)₄ exactly, and R(m) is a multiplicative character sum of Disc(c) and Res(c,m) over affine space: Weil–Deligne territory, and exactly evaluable (Anderson–Evans) when m has two roots. Both identities were verified numerically over F_5 and F_13. Over Z[i] neither identity exists. There the arguments of quartic Gauss sums are the Kummer–Patterson problem, so (L³) is a statement about the joint distribution of Gauss-sum arguments and residue symbols.

**The obstruction in one sentence (S46, PROVED).** For quadratic characters, Poisson maps Möbius sums to Möbius sums. So the same-exponent "quadratic (L³)" follows from GRH applied on the side with fewer rows than length. For quartic characters, Poisson maps Möbius sums to Gauss-sum (T) sums, and the side with fewer rows than length is the T-side. Its objects ψ_α(m³,·) have no Euler product (§10d), so GRH says nothing about them. (L³) is therefore exactly a GRH-strength statement for non-Eulerian Kubota series.

## 11. The proven (sextic) endpoint hits the same wall

The verified submission (θ = 0.8748389) is limited by the detector row at (δ, p) ≈ (0.467, 1/2). In the float model, dθ/dR ≈ 0.135 for the detector count exponent R at that point. The three detector branches have sensitivities, per unit saving δ·η:

| Branch | dθ/dη |
|---|---|
| inverse | −0.043 |
| long | −0.020 |
| affine | −0.011 |

* **The inverse branch.** Its raw exponent is 1 − δr − 2(δx)z, with slot capacity z = (1 − r)/2. That is the maximum allowed by the marked inverse moment's constraint r + 2z ≤ 1 (Cycle25/Assembly/Count/Budget.lean); the second constraint 2r + 8z ≤ 3 is slack. Raising the capacity to (1 − r)/2 + τ would give Δθ ≈ −0.043τ, e.g. θ ≈ 0.87441 at τ = 0.01.
* **Why that is hard.** r + 2z ≤ 1 is the cubic-theta reflection threshold. Going past it is the supercritical marked or mixed moment (MM′, NOTES §8). By the same Poisson count, MM′ has the self-dual obstruction that stops (L³): the dual of a length-U^s column sum has U^{2s−1} rows, so it is shorter only when s < 1.
* **The affine constant B is capped twice over.**
  * The marked high reserve leaves only about 4·10⁻⁶ of headroom.
  * The slot bound (affine capacity < 37/1000) has margin 3·10⁻⁹. Relaxing it all the way to the moment's box edge 1/25 gains only 4·10⁻⁸; the bmax check then binds.

**Summary of this section.** Any substantial further gain, whether through the quartic route or inside the sextic architecture, needs a moment estimate past a self-dual (fixed-point) threshold. The quartic route needs (L³); the sextic route needs a supercritical marked inverse moment. Both are open problems of the same type.

## 12. Concrete sub-problems, roughly from most to least approachable

1. **Optimal-shape quartic large sieve for T-coefficients** at Q = Z (M = 1). Now PROVED on GRH (§10b); open unconditionally:
   Σ_{N(m)≤Z} |Σ_c T(c)χ̄_c(m)W|² ≪ Z^{2+ε}.
   * Numerically true with constant 1.006 (§8).
   * Known only with the extra (QZ)^{3/4}.
   * Through Theorem Z it would give a 7/8 zero-free region for the whole quartic Hecke family over Q(i). That would be new for general Hecke characters of Q(i), though not a better θ for Dirichlet L-functions.
2. **Ratios-type hypothesis (R) in a restricted range**, e.g. Re s = ½ + σ₀ with error H^{3/4 − 3σ₀/2} (§10, S20). Any σ₀ with that error suffices for (G).
3. **A Voronoi formula in the modulus c for quartic Gauss sums** (a GL3-type reflection via the Kazhdan–Patterson theta on the 4-fold cover of GL3). It is not in the literature found here. Even with it, (L³) is the self-dual case of a GL3 thin-family second moment.
4. **(L³) itself, or the supercritical marked inverse moment (MM′) in the sextic setting.** Either would move θ substantially; both are fixed-point problems of the same kind.

## 13. S49: unbalanced quartic large sieve and (L³) at Q^{13/8} — PROVED (l3c/S49_PROOF.md); the first bullets below are the 22:15Z candidate statement, kept for the record
- BGL14's own proof (Lemma 3.3 Poisson step, Lemma 3.4 with ℓ ∈ {1,2,3} for n = 4, plus exact j⁴ multiplicity and a classical large-sieve box)
  appears to give Ξ₄(A,B) ≪ (AB)^ε(A + A^{1/3}B + A^{3/4}B^{1/2}) for B ≤ A ≤ B². This beats BGL's (AB)^{2/3} for every A > B and equals the dFDH
  lower bound for A ≥ B^{6/5}, i.e. dFDH Conjecture 1.3 in that range. It gives nothing at A ≍ B.
- Spot check of a mixed box by hand (x₁ = x₂ = 0.3, η = 0.9): the classical-LS box c4 = 1.2 is what keeps it under the target 1.267, so c4 is essential.
- Possible (L³) consequence: unconditional Q^{13/8+ε}. That is the floor for any large-sieve method, since dFDH's lower bound at (Q, Q^{3/2}) is
  Q^{13/8}. Open items (i)–(iii) are listed in L3_ATTEMPT S49.
- Resolution (22:33Z). Theorem 1 (Ξ₄ ≪ A + B + A^{3/4}B^{1/2} + A^{1/2}B^{3/4} + min(A^{1/3}B, AB^{1/3})) is PROVED from BGL Lemmas 3.1–3.3 plus the new
  Lemma A. Theorem 2 ((L³) sum ≪ Q^{13/8+ε}) is PROVED from Lemma A, Lemma B (column Poisson; its shared-prime case is STANDARD) and Kubota
  convexity. No θ consequence (see S49 resolution in L3_ATTEMPT).

## 14. S50: what stands between Q^{13/8} and Q^{3/2} (l3d/S50_REDUCTION.md)
- Cauchy route plus Poisson in c plus a double Kubota reflection shows that S49's limit Q^{13/8} is exactly the T–theta correlation
  𝒞(h₀) = Σ_c T(c)τ₄(h₀c) at full Cauchy–Schwarz strength. This matches dFDH's extremal (theta-biased) example.
- Any power saving in 𝒞(h₀) for bounded h₀ improves 13/8; square-root size gives Q^{1.543}.
- Numerics (|X|/√S = O(1) over ten twists) indicate square-root size.
- 𝒞 is a quartic Heath-Brown–Patterson-type sign-equidistribution problem, in the opposite (non-WMDS) orientation. Suzuki's known values of
  τ₄ do not reach the index exponent 1 that the duality forces.

## 15. S51: the obstruction series is a G₂ Weyl-group multiple Dirichlet series (l3d/S51_G2.md)
- The generating series of the S50 first moment, Σ_{c,d} α(c)g̃(χ̄_c)·g̃(h₀;χ_d)·χ̄_d(c)N(c)^{−s}N(d)^{−w}, has two Kubota functional
  equations generating W(G₂) (dihedral of order 12).
- It matches the Chinta–Gunnells G₂ series exactly (n = 4, residue symbol ε³, c the short root, d the long root), apart from the
  archimedean α-twist on c.
- Consequences:
  - The theta residue r(s) has its functional-equation centre at s = 3/8, which is square-root size for (P4).
  - Convexity alone gives the (L³) sum ≪ Q^{3/2+ε}, the GRH-level bound. This is conditional only on extending CG to the α-twist and
    Q(i) ⊉ μ₈ and on standard growth bounds.
- Classification by index exponent r (family ν = m^r):
  - first moments: A2 / B2 / G2, all finite;
  - second moments: (2,3,3) finite, (2,4,4) affine, (2,6,6) hyperbolic.
  (L³) is r = 3, the hyperbolic case, which explains why it is out of reach. Cauchy–Schwarz trades it for the finite G₂ first moment at
  the cost of Q^{1/2}.
- Correction to S44: its A3 identification used the index-m series; our family has index m³.

## 16. S52: an unconditional route below 13/8 by mean-value Bochner (l3f/S52_MEANVALUE_BOCHNER.md, S52_PROOF.md) — 16 Oct, 00:50–01:35Z
* **Theorem (modulo S52_PROOF Lemma 4, which has an argument but untracked constants, and the box / T↔T* transfer of S52 §8–§9).**
  Σ_{N(m)≤Q}|ψ_α(m³,½+it)|² ≪ Q^{3/2+ε}(1+|t|)^{O(1)} for every ε > 0 (S52 §18, S52_LIMIT_PROOF.md). No GRH.
  * **CORRECTION (03:38Z, l3f/S52_CORRECTION_NAIVE_VS_CG.md):** this is CONDITIONAL on task #26 (the Chinta–Gunnells G₂ WMDS over Q(i)
    with the BD sign and the α-twist). The naive (L³) double series cannot be continued to the chamber-dcdc hull points (E′ diverges
    there, far outside its region), so the argument must run on the CG series and transfer back at the evaluation point (−1/12, ¼),
    where E′ converges absolutely (checked).
  * The same condition applies to every S52 bound using non-identity chambers (including 1.5573).
  * **04:35Z (l3g/S53_QI_LOCAL.md):** #26's local half is PROVED. The correct sign is the slice sign s^{μ(B_i+1)}, not S51 §19's BD sign.
    What remains is the Ψ-bookkeeping of CG Prop 5.7 over Q(i), labelled STANDARD.
  * Unconditional: Q^{13/8} (S49).
  * 1.5055 comes from the bootstrap chain started at the proved 13/8, with justified (step-up) table reads on a q-grid of step 0.01. The
    LP hulls at the non-binding twists are NUMERICAL.
  * **Q^{3/2+ε} for every ε > 0 (l3f/S52_LIMIT_PROOF.md).** The proof is analytic, given the framework:
    * one explicit two-chamber hull point;
    * an exact cone lemma (cone_exact.py);
    * closed-form S50 bookkeeping on all of q ∈ [⅔, 5/6], with no grids and no LP;
    * the table-offset recursion Δ ↦ (δ′+Δ)/(5.5 + 12(δ′+Δ)), started from S49's Δ₀ = 1/12.
  * The numerical chains (1.5055 at grid step 0.01, 1.5020 at 0.004, 1.5002 with linear reads) agree with it.
  * Q^{3/2} is the Cauchy floor of S50, and it matches what GRH gives (§10b).
* **Method.**
  * S50's Cauchy route needs first-moment bounds (T–theta correlation) for twists h ≤ H = Z^{1/3}.
  * These come from the double series over the 12 W(G₂) chambers, continued by Bochner's tube theorem from *mean-value* regions instead
    of absolute convergence. Each dyadic block takes the better of:
    * the large sieve over c;
    * the mean square over d of the c-sums, using c-FE reflection (fix (a)) and a recursive table S(q) of T*-family exponents.
  * DDHL's clean FE is for ζ(4s−3)ψ, so normalisers restrict evaluation to x, y ≥ 0.
  * An envelope-of-holomorphy sup argument gives polynomial growth, convex h₀-exponents, and vector-valued (RMS-in-h₀) bounds.
  * S50's j⁴-coprimality is handled exactly (effective twists h₀e²), with the true quartic large sieve throughout.
* **Numbers.**
  * The 1.5156 fixed point (s(1) = 2.25) was too optimistic: S52 §14, where the generic twist range at q = 1 forces s(1) ≥ 7/3.
  * After the theta polar term (S52 §16) and the audit fixes of S52 §17 (proved bootstrap start, no table clamping, step-up reads), the
    block regions give 1.5104 → 1.5055 (recur11).
  * The binding mechanism mixes the Region_d boundary point (⅓, 0) across chambers s_c and s_d s_c s_d s_c, with x_r = −1/12 on the (3,1)
    normaliser line and b = ¾η. It closes on the (L³)-shaped block D = C^{2/3}, so the recursion contracts to s(2/3) = 2, i.e. Q^{3/2}.
  * First-pass values (1.5417, 1.5365, 1.5234) were superseded by corrections documented in S52 §10–§12.
* **Why it stops at Q^{3/2}.**
  * The earlier claim that Q^{3/2} needs S(Y,Y) ≪ Y^{2+ε} (S52 §13, §15) only tested the residue-line point of Region_d. The cross-chamber
    combination avoids it (S52 §17).
  * Below Q^{3/2}, S50's Cauchy step is the floor: s(2/3) ≥ 2 from the generic large sieve at the top twist. The second-moment hull gives
    nothing (m-FE → 1/L). (L³) itself needs a non-Cauchy argument.
* **θ.** Unchanged: the verified θ chain is sextic, and the quartic route needs full (L³) (see §7).
