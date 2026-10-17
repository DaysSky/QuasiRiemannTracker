# Large improvements of the quasi-Riemann bound: what the OpenAI architecture needs

Status as of this session. The verified record on the tracker is θ = 0.874861. No new bound is proved here. This note records
(i) exact reductions, (ii) proved identities, (iii) barrier statements, and (iv) the precise open inputs that would give large gains.
Status labels: **PROVED**, **CITED**, **MODEL** (exact exponent algebra of the 7/8 paper's architecture, calibrated), **NUMERICAL**, **HEURISTIC**.

## 1. The architecture and its exponent algebra

The 7/8 proof (OpenAI, 30 Sep 2026) uses a probe I(X,Y,Z) built from the sextic family over K = Q(√−3), with two representations.

- **Low side** (reflection + additive Gram bound): |I| ≪ Z^{x/2 + b/(2k)}, with X = Z^x, Y = Z^y, b = y − x, and k = 6.
- **High side** (Poisson): a principal Mellin signal Z^{C(s)}/L(s,η) with C(s) = s − 1 + x/2 + h/k, plus non-principal rows u.
  The rows are controlled by a zero detector and moment-based row counts.

With compensation length ℓ, row range h = 1 − x + ℓ, and x + y + ℓ = 1 (MODEL, reproduces the paper's 11/12 − ℓ/4):

    σ0 = 1 − 1/(2k) − 3ℓ/(2k)            (k = 6: 11/12 − ℓ/4; k = 4: 7/8 − 3ℓ/8)

The high side at d = h (MODEL; matches the paper's E(h) = −1/48 + (2/3)δ + q/6 − h(1−R) at ℓ = 1/6, x = 17/48):

    E = C0 + (1/2 + ℓ)δ + qℓ − h(1 − R(δ,q)),   C0 = −x/3 − 1/12 + 13ℓ/12   (k = 6)

The generic-row constraint C0 < 0, together with b ≥ 0, forces ℓ < 1/5. This gives the **architecture ceiling**: sextic 13/15, quartic 4/5.

The calibrated numerical model (sextic/fastmodel.py) reproduces OpenAI's endpoint: θ_max ≈ 0.87504, against the exact 7/8.

## 2. Quartic route (family of quartic residue symbols over Q(i))

**Identities (PROVED; NUMERICAL checks, 0 mismatches).** For squarefree primary c ∈ Z[i]:

- g̃(χ_c)² = μ(c)α(c), with α(c) = c/|c|. This follows from J(χ,χ) = −χ(−1)π and twisted multiplicativity.
- μ(c)γ₋₁(c) = ε(c)ᾱ(c)γ₁(c), where ε(c) = ±1 is determined by c mod 4.
- μ(c)g(χ_c,h) = ε(c)α(c)(h/c)₂ g(χ̄_c,h).

So Möbius cancels after Poisson, exactly as in the sextic case.

**Missing lemma (L³).** The rows are Kubota coefficients at index m³ (PROVED). The quadratic twist (m/c)₂ = χ_c(m)² is a Kummer twist
for the 4-fold cover, so it folds into the index. The conductor is N(m)³: a thin family (size Q, conductor Q³). The needed estimate is

    Σ_{N(m)~Q} |ψ_α(m³, ½+it)|² ≪ Q^{1+ε}   (G at M ≥ 2/3).

Assuming it, θ = 5/6 (M = 2/3), or 13/16 with full compensation (ℓ = 1/6).

**Barriers (all bounds below for Σ|ψ(m³)|²):**

| Statement | Status |
|---|---|
| (L³) ⇒ \|ψ_α(m³,½)\| ≪ N(m)^{½+ε} = conductor^{1/6}: Weyl strength. The lossy version still gives conductor^{0.19}. | PROVED |
| Theorem Z: (L³) uniform in Q with loss Q^λ ⇒ ζ and every quartic Hecke L(s,(k/·)₄) over Q(i) zero-free for Re s > 7/8 + λ/2. | PROVED; re-checked |
| Any approximate functional equation + arbitrary-coefficient large sieve gives ≥ Q^{3/2−ε}. | PROVED |
| Gao–Zhao quartic large sieve gives Q^{15/8}. | PROVED |
| Reflection is the involution M ↔ M/(3M−1), with fixed point M = 2/3 = (L³). Reflection + optimal dual sieve reaches only M ≤ ½. | PROVED, given Kubota FE |
| Uniqueness of the sextic: for an order-k family the leftover quadratic twist is Kummer iff the converted character has even order. Only order 3 (k = 6) is usable with a known theta. | PROVED (Stickelberger + Jacobi sums) |
| Spectral (thin Kuznetsov) route stalls at Q²; complete sums are rank-6 hypergeometric (moments 1, 3, 15). | PROVED modulo Katz; NUMERICAL |
| Moment/MDS route: hyperbolic Coxeter group (2,6,6); continuation stalls at convexity; GRH gives only Q^{3/2}. | PROVED modulo tube-domain growth lemma |
| Large-sieve recursion for the Möbius dual: loss ≥ 1/12 (floor), proved 1/6. | PROVED |

**Evidence that (L³) is true (NUMERICAL).**

- α-twisted rows match random phases at Z = 10⁵ for M ∈ [0.54, 1.02] (ratio 0.87–0.93).
- The Möbius dual shows no off-diagonal excess up to D = 10⁵ with H = D^{4/3} (1.000 ± 0.004).
- Untwisted Gauss sums show the predicted bias, which validates the sensitivity of the test.

## 3. Sextic route (improving the existing proof)

**Exact conditional formula (MODEL).** Suppose the row count meets density-hypothesis strength with loss λ, i.e. #rows ≪ U^{1+λ−δ} for all δ. Then

    θ = 11/12 − (1 − 2λ)/(4(5 + 6λ))        (λ = 0: 13/15 = 0.8667; λ = 0.02: 0.8698; λ = 0.04: 0.8728)


### Conditional Theorem D (PROVED, relative to the 7/8 architecture)

**Hypothesis.** In the setting of the 7/8 paper, keep the low-side estimates, the principal-signal estimates, and the small/large-row
estimates. Replace the moderate-row count by a density-hypothesis-type count with loss λ ≥ 0: for every retained bin with
δ = 2a − 1 ∈ (0, κ], every row dyad U = Z^d with d ∈ [2/5, h + ζ], with the floor class lowered to a0 = 1/2 + λ/8 (§50, §53), and
every amplitude class. For d ≤ 2/5 the paper's own no-prime count R = 76/75 − 2δ/3 suffices; it fails only on (0.426, 1/2] near
δ = 2σ0 − 1 (sextic/DH_SEXTIC.md §51, §61). The required count is

    #B ≪ U^{1 + λ − δ + ε}.

**Conclusion.** Let ℓ < (1 − 2λ)/(5 + 6λ), and choose x with

    (13ℓ − 1 + 12(1 + ℓ)λ)/(4 + 12λ) < x ≤ (1 − ℓ)/2.

Then the contradiction argument closes, so β* ≤ 11/12 − ℓ/4. Letting ℓ increase to its supremum gives

    θ = 11/12 − (1 − 2λ)/(4(5 + 6λ)).

**Proof.** Use the paper's high exponent at d = h = 1 − x + ℓ:

    E(h) = C0 + (1/2 + ℓ)δ + qℓ − h(1 − R),   C0 = −x/3 − 1/12 + 13ℓ/12,   q ≤ δ/2.

Substitute R ≤ 1 + λ − δ and q ≤ δ/2:

    E(h) ≤ C0 + hλ + δ(1/2 + 3ℓ/2 − h) = C0 + hλ + δ(x + ℓ/2 − 1/2).

- The δ-coefficient is ≤ 0 iff x ≤ (1 − ℓ)/2.
- With h = 1 − x + ℓ, the condition C0 + hλ < 0 is equivalent to x(1/3 + λ) > 13ℓ/12 − 1/12 + (1 + ℓ)λ, i.e. to the stated lower bound on x.
- An x satisfying both exists iff (13ℓ − 1 + 12λ + 12ℓλ)·2 < (1 − ℓ)(4 + 12λ), i.e. iff ℓ(30 + 36λ) < 6 − 12λ.

So E(h) ≤ C0 + hλ < 0 uniformly in δ, q, with a target-independent margin.

The remaining conditions hold:

- **d < h.** The slope R + δ/2 − 17/50 ≥ 1 − δ/2 − 17/50 > 0, so these exponents are dominated by d = h.
- **Small rows.** The margin is h(17/50 − 1/6) − l_y/2 + 2d_min. It equals −0.0417 at the limit parameters (exact check in
  sextic/DH_SEXTIC.md §61, which re-evaluates every margin of paper §20 there; all stay positive).
- **Low side.**
  - The completed-row excess (5ℓ − 1 + d)/8 − d ≤ 0 for d ∈ [0, ℓ] whenever ℓ ≤ 1/5. This holds because (1 − 2λ)/(5 + 6λ) ≤ 1/5.
  - The Gram condition b = 1 − ℓ − 2x ≥ 0 holds.
- **No slots needed.** The counts above are used without slot amplitudes, so the prime-supply condition is not needed for the count.
  The slots are still used for the compensation. ∎

**Values.**

| λ | θ |
|---|---|
| 0 | 13/15 |
| 0.005 | 0.8675 |
| 0.02 | 0.8698 |
| 0.04 | 0.8728 |

These agree with the numerical model.

**Caveat.** This is the paper's exponent bookkeeping with the count replaced. A full proof would re-run the paper's order-of-choices
argument, which carries over verbatim because all margins are uniform.

**Single-input improvements are weak (MODEL).** Each input below was pushed to an extreme in isolation:

| Input changed | θ_max |
|---|---|
| Plain capacity 6κ → 1 | 0.8736 |
| Inverse amplification α: 5/6 → 1/2 | 0.8723 |
| Inverse capacity doubled | 0.8734 |
| All three together | 0.8702 |

The binding structure is the two-witness zero detector, whose worst case is max(r, 2m) ≥ 2/3.

**Mixed moment (MODEL + PROVED reduction).**

- Restricted to r + m ≤ 1, it gives **no gain** (0.875058). The detector allows r + m ∈ [t, t + ½], and the adversary escapes.
- The full version MM′ (all r ≤ 1, m ≤ ½) would give the density hypothesis and hence 13/15.
- After Poisson, length s = r + m has U^{2s−1} dual rows and U^s columns. Cubic-theta reflection shortens the dual only if
  r + 2m < 1. This is exactly the paper's marked constraint r + 2z ≤ 1, so the factor 2 is intrinsic.
- Detector pairs always have r + m ≥ 1, so the method cannot reach them.

**Required zero-density exponent (MODEL; A defined by #rows with a zero of Re ≥ σ ≪ U^{A(1−σ)}, density hypothesis A = 2):**

| Target θ | Required A, on σ-range |
|---|---|
| 0.8745 | ≤ 2.155–2.28, on [0.70, 0.80] |
| 0.874 | ≤ 2.14–2.25, on [0.68, 0.82] |
| 0.872 | ≤ 2.09–2.21, on [0.62, 0.86] |
| 0.870 | ≤ 2.065, on [0.58, 0.86] |

**Calibration (HEURISTIC analogy).** For ζ, the best known exponent around σ ∈ [0.70, 0.78] is Guth–Maynard's 30/13 ≈ 2.31. The density
hypothesis is known only for σ ≥ 25/32 (Bourgain). In the sextic family, the Halász off-diagonal involves χ_u χ̄_{u′} = χ_{uu′⁵}, whose
conductor is about U², so large-values tools are weaker than for ζ.


**Mollifier calibration (HEURISTIC, Selberg's method N(σ) ≪ U^{1−θ_m(2σ−1)}).**

- Target 0.874 needs θ_m ≥ 0.878 at σ = 0.70, θ_m ≥ 0.90 at σ = 0.75 and θ_m ≥ 0.916 at σ = 0.80, i.e. a mollifier of length ≈ U^{0.9}.
- Known mollified second moments reach θ < 4/7 for ζ (Conrey) and θ < 1/2 for Dirichlet characters mod p (Michel–VanderKam).
- MM′ is the family analogue of a full-length (θ = 1) mollifier.

**Proven-input limits (PROVED, by exponent algebra).** With the paper's inputs, no rearrangement beats max(r, 2m):

- Hölder between the inverse moment and the plain fourth moment gives count 1 − δ·max(r, 2m).
- Reflecting the plain witness keeps the lower bound U^{δm} and lengthens it.
- An inverse-only detector lets the adversary push r > 1, where e(r) − δr = 1/6 + (5/6 − δ)r grows.
- A marked inverse moment with integer marks k·k′⁵ needs r + 4m ≤ 1, while the detector forces r ≥ 1 − m.
- Cubic-theta reflection with the plain factor as a mark needs r + 3m < 1 (twist conductor N(h)N(k)), versus the paper's r + 2z ≤ 1 for pre-summed prime slots.

**Model confirmation.**

- MM′ gives 13/15 exactly.
- MM′ with loss U^{0.02} gives 0.8695.
- MM at critical length only gives 0.875058 (no gain).

## 4. Open inputs that would give large gains

1. **(L³)**, Lindelöf-on-average for quartic Kubota coefficients at cube indices. It gives θ ≤ 5/6, or 13/16 with full compensation.
   It implies Weyl-strength subconvexity and a quasi-RH (7/8) for quartic Hecke L-functions over Q(i).
2. **Density-hypothesis-strength zero density for the sextic family on σ ∈ [0.6, 0.85]** (equivalently MM′). It gives θ → 13/15.
3. **A reflection for products of Möbius and plain sums beyond r + 2m < 1** in the cubic-theta setting. This is the concrete mechanism for (2).

## 5. Files

- quartic/NOTES.md: full log.
- quartic/*.py: identities, dual mean squares. Subfolders condB, po, spectral, moments, barrier: subagent work.
- sextic/*.py: exponent model, mixed moments, density requirements.

## 6. Direct attempt at the density hypothesis for the sextic family

**Re-derivation of ρ(δ) (PROVED, exponent algebra).** For a zero in bin δ, take the detector tail block n ~ U^s with s ≥ t, in two representations.

- (i) l ≤ U^t: pair (r ≤ 1, m = s − r). Separate moments give 1 − δ·max(r, 2m) ≥ 1 − (2/3)δs.
- (ii) Complement, l > U^t: pure Möbius of length s. The amplified inverse moment e(s) = (1 + 5s)/6 gives 1/6 + (5/6 − δ)s.
- The adversary balances at s* = 5/(5 − 2δ), giving ρ(δ) = 1 − 10δ/(3(5 − 2δ)). This is exactly Cycle61's Lemma 2.1, and it is valid for all δ.
- At s = 1, branch (ii) already gives the density hypothesis. All loss comes from blocks s > 1.

**Mollified approximate-functional-equation detector (PROVED, exponent algebra).** From 0 = M_X A + ε M_X B there are three branches.

- The root-number branch is forced to plain length exactly ½ by the pointwise bounds; the plain fourth moment then gives the density hypothesis.
- The E1 branch has m ≥ ½.
- The E2 branch is the original adversarial pair, unchanged. No gain.

**Exact sufficient improvements.**

- (a) **Amplification slope.** If e(r) = 1 + α′(r − 1) with α′ ≤ 0.55, θ = 0.874 follows at δ = 0.4. Currently α′ = 5/6 = 1 − 1/6, from sixth powers.
  Cube lifts leave a quadratic twist (a/n)₂, so they do not amplify.
- (b) **MM′.**
- Neither follows from proven inputs.

**Why ζ differs (structural).** For ζ, products at length ≤ T are covered by the perfectly orthogonal mean value theorem. The sextic large sieve
gives U^{5/3} at length U. OpenAI's special-coefficient orthogonality covers μ alone and plain × plain only.

**Type I large-sieve attempt.**

- Reflecting the long plain factor through the root number, with |ε_u| = 1, gives a sextic large sieve over LK characters.
- That sieve is orthogonal only when LK ≤ U^{1/5}.

**Conclusion.** The density hypothesis here contains the s > 1 problem that keeps it open for ζ below σ = 25/32, with weaker tools.

## 7. Breakthrough search in the October 5 framework (11/12 via family mean square)

The October 5 paper reduces 11/12 to MS(H): Σ_{N(u)≤H} |A_u(D)|² ≪ D^{1+ε}H, where A_u(D) = Σ μ(n)ν(n)χ_n(u)W(N(n)/D) and H = D^{1+ϑ}.
It then extracts using the H^{1/6} copies u = p⁶: |A_1|² ≪ D·H^{5/6}.

**(a) Copy extraction is optimal (PROVED, elementary).** Take any linear amplifier a(u) supported on N(u) ≤ H, and let
K(n) = Σ_u a(u)χ_n(u) ≈ c on squarefree n.

- The mean of K over n is c = Σ_v a(v⁶), since the other sixth-power classes average to 0.
- Hence ‖a‖² ≥ Σ_v |a(v⁶)|² ≥ c²/H^{1/6}.
- So Cauchy–Schwarz against the box mean square never beats |A_1|² ≤ H^{-1/6}·Σ_u |A_u|².
- For an order-k family, θ = 1/2 + (1 − 1/k)·log H/(2 log D), and the bound is no better than this.

**(b) The threshold H ≥ D is the real wall (PROVED, exponent algebra).**

- Including n = n′, Poisson in u gives Σ_u |A_u|² w(u/H) = (H/D)·Σ_{h≠0} ŵ(hH/D²)·|B_h(D)|². The diagonal DH sits inside this sum.
- With absolute values, the right side is ≍ (H/D)·(D²/H)·D = D². So DH needs H ≥ D.
- For H < D one must use that ∫ŵ = w(0) = 0 (annular w). That is a mean-zero cancellation in Σ_h ŵ|B_h|².
- Poisson in h returns exactly the original u-family (an involution), so this route is circular.
- Theta reflection shortens the dual only when 𝓗 = D²/H ≤ D.
- So MS(D^c) with c < 1 (which would give θ = 1/2 + 5c/12) is a beyond-large-sieve, beyond-density-hypothesis statement.

**(c) Degree-of-reflection meta-formula (MODEL).**

- An order-k family whose Gauss-sum dual reflects through a degree-d theta has dual length 𝓗^d/X. This shortens only for H ≥ D^{2−2/d}.
- That gives θ = 1/2 + (1 − 1/d)(1 − 1/k).
- d = 2, k = 6: 11/12. d = 2, k = 4: 7/8. d = 2, k = 2: 3/4.
- d = 3, k = 4: θ = 1, trivial. So a GL(3) quartic theta would not help even if its coefficients were known.
- CITED: Kazhdan–Patterson. The 4-fold cover of GL(3) has a unique Whittaker model only for certain covers. On GL(2), uniqueness holds only for n ≤ 3.
- So d = 2 with k < 6 needs a GL(2) theta of order 4 or 5, whose coefficients are unknown. That is (L³).

**(d) Weighted cube family {u·v³} (NEGATIVE).**

- It has F^{1/6}·H2^{1/2} copies, which formally would give 7/8 at H2 = D^{1/4}, H1 = D^{3/4}.
- But for each fixed v, the u-family is a short family (H1 < D). By (b), its dual diagonal costs D², and the gain is lost exactly.
- The 3-variable family u1u2u3, whose Gauss-sum cube is a Hecke character, has copy density H^{-2}: trivial.

**(e) Part II generic-row ceiling.** Even if every non-target row satisfied GRH (δ = 0, full count), C0 < 0 forces ℓ < 1/5.
So the probe architecture cannot go below 13/15 without a new way to bound generic rows other than absolute values.

**(f) Landau–Siegel paper (OpenAI, 1 Oct 2026).** It uses an interpolation-determinant method with Frobenius at primes where χ(p) = −1.

- It needs a log-mass of these primes above 3/4.
- A zero with 1 − β fixed only shifts this mass by O(1/(1 − β)) out of log X.
- So the method does not extend to fixed-width half-planes.

**Net.** Every route to a large decrease of θ reduces to one of three inputs:

1. (L³) for quartic Kubota coefficients: Weyl-strength subconvexity.
2. Density-hypothesis-strength zero density for the sextic family, giving at most 13/15 in the Part II architecture.
3. A beyond-large-sieve mean square MS(D^c) with c < 1.

No proof of any of these was found.

## 8. Direct attack on the density hypothesis for the sextic family

### 8.1 Exact target (NUMERICAL, sextic/needcurve.py)

The bound needed is #rows with a zero of real part ≥ σ ≪ U^{A(1−σ)}. The table gives the largest admissible A at each σ.

| Target θ | σ = 0.65 | σ = 0.70 | σ = 0.75 | σ = 0.80 | σ = 0.85 |
|---|---|---|---|---|---|
| 0.8745 | — | 2.175 | 2.218 | 2.282 | — |
| 0.872 | 2.102 | 2.114 | 2.132 | 2.160 | 2.205 |
| 0.870 | 2.066 | 2.067 | 2.068 | 2.069 | 2.071 |

The current count ρ(δ) gives A = 2.275 at σ = 0.70, 2.333 at 0.75 and 2.368 at 0.80.

For comparison, Guth–Maynard for ζ gives 15/(3+5σ): 2.31 at 0.70, 2.22 at 0.75, 2.14 at 0.80. So even ζ-quality bounds miss the 0.8745 target at σ = 0.70.

### 8.2 The block calculus (PROVED, exponent algebra)

**Detector.** A zero in bin δ forces a block M_r·S_m to be large at level U^{δs/2}, where s = r + m. The constraints are r ≤ t, m ≤ 1/2 (after the functional equation) and s ∈ (t, t + 1/2].

**Counts from Hölder.** Using a p-th moment for the Möbius witness and a q-th moment for the plain witness (both Lindelöf-on-average) gives the count

    1 − δ·s·pq/(2(p+q)).

- The density hypothesis (count ≤ 1 − δ) needs s·pq/(p+q) ≥ 2 for every block the adversary can choose.
- A p-th Möbius moment needs (p/2)·r ≤ 1, i.e. length below the conductor. A q-th plain moment needs (q/2)·m ≤ 1.
- **Available inputs:**
  - Möbius: p = 2 for r < 1.
  - Plain: q = 4 at every length.
  - Marked inverse: still a second moment of M·Q, so p = 2.
- **Blocks with s < 1.** These cannot reach the density hypothesis even with (4,4). It would need p = 2/r and q = 2/m, i.e. fractional critical moments. This is the classical integer-power obstruction, which for ζ is bridged only by Halász–Montgomery in limited σ-ranges.
- **Hence t ≥ 1, so r ≥ 1/2.** A Möbius 4th moment at r ≤ 1/2 never enters. That moment is the Möbius×Möbius moment, whose dual is the twisted A2 cubic Weyl group multiple Dirichlet series. It follows from the canonical marked family via Cauchy–Schwarz when r1 + 2r2 ≤ 1.
- **The worst blocks** are r ∈ (1/2, 1], m = s − r, s ∈ [1, 3/2). There (2,4) gives 1 − (2/3)δs. Every other usable pair needs a Möbius moment of total length > U, or a mixed moment with s > 1. Both are beyond-conductor objects.

### 8.3 Large values (Halász–Montgomery) for this family (PROVED, exponent algebra)

- The Gram entries are plain sums of χ_u·conj(χ_u′), with conductor up to U².
- The pairs family {χ_u·conj(χ_u′)} is essentially all sextic Hecke characters with local exponents {1,5} and conductor ≤ U², of size about U².
- Even with a Lindelöf-on-average 4th moment over that family, the Huxley-type bound is R ≪ N²V^{-2} + N³U·V^{-4}. At N = U this is U^{2−2δ}, worse than the density hypothesis.
- Pointwise Lindelöf for the pairs would give R ≤ N^{2−2σ} for σ ≥ 3/4. That still fails for the supercritical blocks.

### 8.4 Where the density hypothesis stops

Every route needs a Lindelöf-on-average mean value for a Möbius-containing polynomial longer than the conductor:

- e(r) = 1 for some r > 1; or
- MM′ (the mixed Möbius×plain moment at total length s > 1); or
- a Möbius 4th moment for r ∈ (1/2, 1].

In the theta toolkit these sit outside the reflection's shortening range (r + 2m < 1, or 3r − 2 < r). In large-values terms, the needed range σ ∈ [0.63, 0.79] lies below Bourgain's DH range for ζ (σ ≥ 25/32 ≈ 0.781). At σ = 0.70 the target exponent 2.175 is also below Guth–Maynard's 2.31.

Not proved here.

## 9. New reduction: the Möbius fourth moment implies the density hypothesis (and 13/15)

**Reduction (PROVED in the model; follows from the paper's simultaneous-saturated-witness detector).**
Let M4(r₄) be the statement Σ_{N(u)≤U} |M_u(U^r)|⁴ ≪ U^{1+ε} for all r ≤ r₄.

- The detector produces blocks with r ≤ t, m ≤ 1/2, r + m ≥ t, and both witnesses are saturated simultaneously.
- With M4, the Möbius witness count is 1 − 2δr. The plain fourth-moment count is 1 − 2δm.
- So the block count is 1 − 2δ·max(r, m), and at t = 1 the adversary's best is r = m = 1/2, giving 1 − δ: the density hypothesis.

**Model values** (sextic/m4model.py; calibrated baseline 0.87501):

| r₄ | ≤ 0.667 | 0.70 | 0.75 | 0.80 | 0.90 | 1.00 | 1.00, loss U^{0.03} |
|---|---|---|---|---|---|---|---|
| θ | 0.87501 | 0.87501 | 0.87425 | 0.87212 | 0.86918 | **0.86667 = 13/15** | 0.86864 |

**Numerics (sextic/mmnum/m4.py; U up to 3·10⁵, about 10⁶ rows).**

- The ratio M4/(2·M2²) at r = 1 grows 1.054 → 1.062 → 1.083 → 1.091 (U = 10⁴ … 3·10⁵).
- The random multiplicative control behaves the same: 1.043 → 1.065 → 1.076 → 1.090. Independent random signs give 1.005.
- The largest single row carries at most 0.2% of the moment.
- So μ behaves like a random multiplicative function, the excess is the multiplicative (log-type) factor, and M4 holds numerically with no power loss.

**Why the existing toolkit stops at r = 1/2 (PROVED, exponent algebra).**

- **Dual structure.** Poisson in u turns M4 into Σ_h ŵ(hU/D²)·|B_h|² with D = U^{2r}. Here B_h = Σ_N ᾱγ₂(N)·τ_w(N)·χ̄_N(h): cubic-theta coefficients weighted by a smooth divisor function, i.e. the twisted cubic A2 Weyl-group multiple Dirichlet series (H(c₁,c₂) = g(c₁)g(c₂)(c₁/c₂)₃^{-1}).
- **Size.** With absolute values the dual is 𝓗·D = D³/U. The target is D², so r > 1/2 needs cancellation from ∫ŵ = w(0) = 0.
- **What that cancellation amounts to.** For r ≤ 3/4 it is enough to have square-root cancellation of the fluctuation |B_h|² − D against ŵ at scale 𝓗. Poisson in h is the involution back to the u-family, and theta reflection lengthens the sums when r > 1/2.
- **Cauchy–Schwarz over one factor** with the paper's canonical marked estimate gives r₁ + 2r₂ ≤ 1, i.e. r ≤ 1/3 in the symmetric case.
- **Sixth-power lift from the critical case** costs (5/6)(2r − 1).

**Next route (spectral, outside the Poisson/theta group).**

- Identify B_h with Fourier coefficients of a metaplectic Eisenstein series on the cubic cover of GL(3); A2 is its Whittaker coefficient (Bump–Hoffstein, Brubaker–Bump–Friedberg).
- Obtain Σ_h ŵ(h/𝓗)|B_h|² by Rankin–Selberg unfolding.
- **Needed analytic input:** R(s,s′;w) = Σ_h Z(s;h)·Z̄(s′;h)·N(h)^{-w} is holomorphic with polynomial growth on Re w ≥ α, apart from w = 1 + i(t − t′), with α ≤ 2r/(4r − 1).
  - That is α ≤ 3/4 for r = 3/4, 0.727 for r = 0.8, and 2/3 for r = 1.
  - The functional-equation (WMDS) group here is not of finite type, so this has to come from spectral theory.

## 10. Second round on the density hypothesis: both witness sides

### 10.1 Correction to the spectral route
B_h = Σ_N ᾱγ₂(N)·τ_w(N)·χ̄_N(h) splits as follows:

- the cubic part χ_N(h)² is a numerator twist g₃(h̄; d), i.e. a Whittaker index;
- the quadratic part (h/N)₂ is an extra quadratic Hecke twist of the inducing data.

So B_h is the (h̄, h̄) Whittaker coefficient of the cubic-cover GL(3) Eisenstein series twisted by (h/·)₂. The h-family is a family of quadratic twists of a GL(3) object, not a sequence of Fourier coefficients of one form. Rankin–Selberg unfolding does not apply as stated.

- A GL(3) reflection (conductor about N(h)³) shortens the dual only when D < U^{3/4}, i.e. r < 3/8.
- Cauchy–Schwarz plus the canonical estimate gives r ≤ 1/3.

### 10.2 The critical-moment form of the density hypothesis (PROVED, exponent algebra)

- At t = 1 the worst block is (r, 1 − r) with r ∈ (1/2, 1).
- The density hypothesis needs, for each such r, the critical Möbius moment Σ_u |M_u(U^r)|^{2/r} ≪ U^{1+ε} (length × exponent/2 = 1).
- Equivalently, plain moments of every order.
- Multiple mollifier lengths add no constraint: the adversary reuses witnesses.
- Saturation does not change the (2,4) count of 1 − (2/3)δs.

### 10.3 The plain side (MODEL, sextic/p6model.py)

| Plain input added | θ |
|---|---|
| Sixth moment for m ≤ 1/3 (critical) | 0.875005 (no gain; the adversary moves to m slightly above 1/3) |
| Sixth moment for m ≤ 0.4 | 0.874597 |
| Sixth moment for m ≤ 0.5 | 0.872743 |
| Eighth moment for m ≤ 0.5 | 0.871425 |

All of these are supercritical: the product of three plain sums has length 3m > 1. The paper's two-factor plain lemma is length-free only because each factor can be reflected to length ≤ 1/2 by the functional equation, so two factors stay critical.

### 10.4 Summary of conditional targets

| Input | θ |
|---|---|
| Möbius fourth moment, r ≤ 1 | 13/15 = 0.86667 |
| Möbius fourth moment, r ≤ 0.8 | 0.87212 |
| Möbius fourth moment, r ≤ 0.75 | 0.87425 |
| Plain eighth moment, m ≤ 1/2 | 0.87143 |
| Plain sixth moment, m ≤ 1/2 | 0.87274 |
| Plain sixth moment, m ≤ 0.4 | 0.874597 |

Each is a supercritical Lindelöf-on-average statement for a thin family (family size about equal to conductor). The classical analogues for characters mod q, namely the full-length mollifier fourth moment and the sixth moment of L(1/2, χ), are open as well.

## 11. Loss tolerances and where the paper's own plain lemma stops

### 11.1 Weakened higher moments (MODEL; sextic/lsmodel.py, liftmodel.py, p6loss.py)

**No gain at all (θ stays at 0.875005):**

- large-sieve-shape bounds Σ|S_m|^{2k} ≪ U^{max(1,km)} and Σ|M_r|⁴ ≪ U^{max(1,2r)};
- critical-length moments plus the sixth-power lift, i.e. a supercritical loss slope of 5/6 per unit of excess length. This holds for plain 2k-th moments with k = 3..8 and for the Möbius fourth moment.

**Plain sixth moment, m ≤ 1/2, with a constant supercritical loss λ:**

| λ | 0–0.02 | 0.04 | 0.06 | 0.08 | 0.10 | 0.12 | 0.15 |
|---|---|---|---|---|---|---|---|
| θ | 0.872743 | 0.872967 | 0.873736 | 0.874088 | 0.874411 | 0.874721 | 0.874983 |

- Every λ ≤ 0.12 beats the record 0.874861. Only m ≤ 0.45 is needed.

**Plain sixth moment with loss c·(3m − 1):**

| c | 0.6 | 0.5 | 0.4 | ≤ 0.3 | 5/6 (lift) |
|---|---|---|---|---|---|
| θ | 0.874421 | 0.873356 | 0.872938 | 0.872743 | no gain |

- With the plain eighth moment added at c = 0.3: 0.871526.

### 11.2 Why the paper's plain-lemma recursion cannot be pushed past A = M

- In the second transform the diagonal j = 0 exceeds the allowance by exactly A − M (paper eq. 2.12).
- This positive diagonal comes from Cauchy–Schwarz in the first dual variable h, which separates the two rectangles, followed by enlarging the row ball by positivity. That step discards the mean-zero oscillation of ŵ.
- For two plain factors this is harmless, because the functional equation reflects each factor to length ≤ M/2, so A ≤ M (the "padded zero-slot core").
- Three factors stay supercritical (A up to 3M/2). So the recursion gives at best the large-sieve shape, which by 11.1 gains nothing.

### 11.3 The sharpened target

The smallest input known to beat the record is:

> Σ_u |S_u(U^{m₁}) S_u(U^{m₂}) S_u(U^{m₃})|² ≪ U^{1+λ}, with m_i ≤ 0.45 and supercritical loss λ ≤ 0.12, or slope c ≤ 0.6.

It needs a step that keeps the ŵ-oscillation in the first dual (no Cauchy–Schwarz in h). The closest precedent is Li's unconditional second moment of quadratic twists of GL(2) L-functions, also a supercritical second moment in a thin family. There the dual collapses, through Hecke relations, to short L-functions. Here the plain dual is θ·μτ₃, a Möbius-type sum, so that collapse is not available in the same form.

## 12. The low side's exceptional Gram term is the main-term carrier

### 12.1 Numerics (lowside/bias.py, gram.py, gauss6.py, gauss6b.py)

- With c(s) = 1 over all squarefree s (split primes), Gram/diag ≈ 0.95–0.97 for P_a from 5 to 82. There is no Kubota bias in Σ_s W g̃(χ_s) χ̄_s(−m).
- This does **not** show that the paper's term is an artifact. The excess cancels here only because χ_s(−1) and the reciprocity phase χ_{s'}(s)χ̄_s(s') change with the class of s. The paper fixes a ray class σ, where both phases are constant.

### 12.2 The exceptional term is real for the paper's coefficients (PROVED, local computation)

- The paper's coefficient is c_σ(s) = 1_σ(s)·χ_s(b_*)/ξ(s).
- At an odd prime p of S it restricts to χ_{π_p}·ξ̄_p (sextic × quadratic), a character of order 3.
- The exceptional frequency k carries κ_k(m) ∋ χ_{π_p}(m)^{e_p}. The complete mean is nonzero exactly when e_p ≡ 2 (mod 6) at every such p.
- So the exceptional k form a fixed translate k₀·a⁶ of the sixth powers. The C = d = 1 term has fixed sign and gives (Q/Y′)·P_a^{1/6}. Under the paper's coefficients the P_a^{1/6} term is a genuine lower-order energy, not slack.

### 12.3 A canonical projector (answers the catalogue's "A^(6) undefined")

- Define A_θ(m) as the polar part of A_m(Y′) in the s-variable: the residue at w = 2/3 of Y′^w·Ŵ₁(w)·Σ_{s∈σ} c_σ(s) g̃_{χ_s}(s,−m) q_s^{−w}.
- This is the Kubota-type series of the twisted sextic metaplectic Eisenstein series. It has a pole because c_σ restricted to sixth powers is 1_σ(t⁶).
- So A_θ(m) = Y′^{−1/3}·r_σ(−m)·Ŵ₁(2/3+iν), where r_σ is the twisted sextic theta coefficient.
- The exceptional Gram energy corresponds to Σ_{N(m)~Q}|r_σ(m)|² ≈ Q^{5/6}, i.e. Y′^{−2/3}Q^{5/6} = (Q/Y′)P_a^{1/6}.
- The definition depends on no pair (s₁, s₂). It is a fixed linear map applied to A, as the catalogue (F04.P, A14) asks.

### 12.4 Why removing the term is at least as strong as the theorem (exponent algebra plus the projector)

- The theta part of the probe is I_θ = Q^{−1/2}·Y^{−1/3}·J(Q,Z), with J = Σ_m ξ(m) r_σ(−m) (N(m)/Q)^{−iν} Ω B_m(Z).
- It depends on X and Y only through Q = q_{b*}XY.
- The paper's binding inequality C(θ) > x/2 + b/12 is exactly θ > 5/6 − ℓ/6 + q/12 with q = log_Z Q, a function of q alone. It is the Cauchy–Schwarz bound |J| ≤ ‖r‖·‖B‖ = Q^{5/12}·Q^{1/2}, set against a main term of size Q^{5/6}·Z^{s−5/6+ℓ/6}. (The main term X^{1/3}Y^0 has exactly the scaling Y^{−1/3}Q^{1/3}.)
- The row-energy threshold ‖B‖² ≈ Q needs Q ≥ Z^{1−ℓ}. So q = 1 − ℓ, which gives θ = 11/12 − ℓ/4, the paper's value.
- Fix Q and let Y → ∞ (so X = Q/(q_{b*}Y)). In the high representation the Y^{−1/3} coefficient comes from singularities on the hyperplane w + z = 7/6, and the principal double pole (w,z) = (1,1/6) lies on it. So the main term X^{1/3}Z^{C(s)}/L(s,η) sits in the theta part I_θ.
- The non-principal rows contribute to J only through the continuation of the row sum Σ_u, which is the same object the paper controls by zero density. So J = κ·S(Z) + (row terms), with κ ≍ Q^{5/6} and S(Z) the η-twisted Möbius-type sum (Mellin transform H(s)/L(s,η)).
- **Conclusion.** Any bound for the theta-paired part better than Cauchy–Schwarz must in particular bound the main term κ·S(Z). It is therefore at least as strong as the zero-free statement it is meant to prove. The A14 / α-lever ("remove b/12") is not an independent input. My earlier model values 0.8696 (α = 1/2) and 0.8642 (α = 1) are not reachable through the low side.
- The probe works as a theta-dual large sieve: |S(Z)| ≤ ‖r‖‖B‖/κ = Q^{1/12}. The loss Q^{1/(2n)} is the alignment defect of the n = 6 theta. Linear amplification of the test vector is already excluded (§7, copy extraction is optimal).
- The remaining levers are q < 1 − ℓ (the row energy below its threshold, i.e. the Möbius/plain moments of §§9–11) and larger ℓ (the high-side rows).

## 13. Model of the Cycle61 record and its binding constraints (c61/c61model.py, opt.py, opt2.py)

### 13.1 The model reproduces the record

- **Low side** (Cycle61 eq. 5.5): θ ≥ 5/6 − ℓ/6 + M/12 + (1 + 3ℓ − 2M)₊/8, with M = x + y. It is valid under M + ℓ ≤ 1, ℓ < x ≤ y and y − 2x + ℓ < 0.
- **Rows** (eq. 3.14): θ ≥ max over d ≤ h, δ ≤ κ, p ≤ 1/2 of G_κ(d,δ,p), where
  G = (1+δ)/2 + 13h/75 − y/2 − ℓ/2 + ℓδp + d(1/4 + 3R_κ/4 − 17/50) − δpC_d/2.
- At the record geometry the low side is exactly tight (θ = 11/12 − ℓ/4 = 0.874861). The rows give 0.874850 at d = h, δ ≈ 0.475, p = 1/2.

### 13.2 The cap M + ℓ ≤ 1 is binding

- Rows want y large: dG/dy ≈ −1/2 − δp/Λ ≈ −0.56. The low side charges only 1/12 per unit of y. So the record sits at y = h − 2ℓ, i.e. M + ℓ = 1.
- **Origin of the cap** (paper, proof of Lemma probe-row-norm): T_d ≤ H + (M + ℓ − 1). Past the cap the reflected dual is longer than the row range, and the hybrid large sieve charges max(H, v + ℓ_b). The row energy then grows by (M + ℓ − 1), i.e. by (M + ℓ − 1)/2 in θ.
- **Optimised θ** if the excess cost were c·(M + ℓ − 1):

| c | 0 (free) | 0.05 | 0.10 | 0.15 | 0.25 | 0.5 (actual) |
|---|---|---|---|---|---|---|
| θ (grid model) | 0.87148 (ℓ → 1/5) | 0.87367 | 0.87491 (≈ cap) | cap | cap | cap |

  (Grid model with the cap gives 0.874915; the record's finer tuning gives 0.874861.)
- Break-even is c ≈ 0.1, against the natural 1/2.
- Removing the excess means bounding the mean square of the theta rows beyond the reflection range. By Poisson duality this is the Möbius fourth moment of §9 at r ≈ 0.7. So the cap is the same input again.

### 13.3 Strength of the remaining inputs

| Input | Rough strength |
|---|---|
| M4 at r ≤ 1/2 | Fourth moment of L-functions in the family (critical; follows from the large sieve) |
| M4 at r = 3/4 | Sixth moment |
| M4 at r = 1 | Eighth moment |
| Plain sixth moment, m_i ≤ 0.45, loss ≤ 0.12 | Below ζ's known sixth-moment loss of 1/4 |

- Even for the quadratic family (Heath-Brown's optimal large sieve), an unconditional Lindelöf-on-average sixth moment is not known.
- Single-modulus Möbius variance (the slope α′ < 5/6 of §6) is a power-saving averaged Chowla statement over shifts kq.
- So every lever found so far is at least sixth-moment-strength for a thin family.

### 13.4 Conditional values in the record's own framework (c61/c61m4.py, c61m4opt.py)

The count model (min over t of max of the short and long envelopes) reproduces Cycle61's closed form R_κ to 10⁻³. Adding the Möbius fourth-moment witness count 1 − q − (2δ − q)r for r ≤ r₄ and re-optimising the geometry (cap kept):

| Input | θ |
|---|---|
| none (record inputs, grid) | 0.87485 (record 0.874861) |
| M4, r ≤ 0.75 | 0.87464 |
| M4, r ≤ 0.85 | 0.87121 |
| M4, r ≤ 1 (= density hypothesis count) | 0.86675 ≈ 13/15 |
| Density hypothesis for every row (R = 1 − δ) | 0.86675; ℓ → 1/5, then generic rows plus the cap bind |

- Under the density hypothesis the rows have about 0.035 of slack. The ceiling then comes from ℓ < 1/5, which is the generic rows (δ = 0) together with y ≤ h − 2ℓ (the cap).
- So below 13/15 the architecture needs both the density hypothesis and a row-energy bound past the cap (§13.2). Even with every input, the low side gives θ ≥ 5/6 − ℓ/6 + M/12 ≥ 5/6.

## 14. Attempt on the Möbius fourth moment (M4)

M4(r) means Σ_{q_u≤U} |M_u(U^r)|⁴ ≪ U^{1+ε}, with M_u(D) = D^{−1/2} Σ_n μ(n)ν(n)χ_n(u)^{±1} W(q_n/D).

### 14.1 What is actually proved

- **r ≤ 1/3, rigorous.** Write the second factor M_u(D) as a finite sum of products of prime slots of total length r (μ factorises over prime dyads). The paper's marked inverse moment (Lemma marked: r + 2z ≤ m, 2r + 8z ≤ 3m, with z = r, m = 1) then gives r ≤ 1/3.
- **Not proved for r ∈ (1/3, 1/2].**
  - Writing M_u² as one polynomial of length U^{2r} meets the BGL large-sieve term (U·U^{2r})^{2/3}, which exceeds U once r > 1/4.
  - Avoiding that term needs the canonical marked estimate (Lemma canonical-moment) for columns ᾱγ₂(c)·t(c), with t(c) = Σ_{c₁c₂=c} W W (c₁/c₂)₃. That is the A2 cubic double series. The paper's lemma explicitly excludes any residual column coefficient.
  - Even granted, the canonical regime requires dual rows ≤ columns: D⁴/U ≤ D², i.e. r ≤ 1/2.

### 14.2 Why r = 1/2 is a wall for the Poisson/theta toolkit (PROVED, exponent algebra)

- **Self-duality.** Poisson in u sends rows U and length D² to dual rows 𝓗 = D⁴/U and columns D². The two sides exchange roles exactly at D² = U.
- **Reflection lengthens past r = 1/2.**
  - Cubic theta reflection of a column of length C at row conductor q_k has length q_k²/C. In the dual this is 𝓗²/D² = D⁶/U², which exceeds D² iff r > 1/2.
  - The A2 functional equations reflecting c₁ at fixed c₂ have conductor q_k·q_{c₂}, so they lengthen too.
- **Bilinear (Type II) pieces.** For general coefficients, Σ_u |A_u(M)B_u(N)|² ≪ U||α||²||β||² is false once MN > U. Take α, β conjugate-aligned with ψ_{u₀}. So any argument must use the arithmetic of μ inside the bilinear pieces.
- **Classical analogue.** For all characters mod a single q, M4 at r = 3/4 is the Lindelöf-on-average sixth moment for an individual modulus, which is open. M4 at r slightly above 1/2 is the individual-modulus variance of μ*μ in residue classes at x = q^{1+δ}, which is a power-saving averaged-Chowla statement over shifts kq.

### 14.3 Range needed (Cycle61 framework, §13.4)

| r₄ (M4 assumed for r ≤ r₄) | ≤ 0.72 | 0.75 | 0.85 | 1 |
|---|---|---|---|---|
| θ | no gain (0.87487, grid baseline) | 0.87464 | 0.87121 | 0.86675 |

- M4 has to reach r ≈ 0.73 before the record moves at all.
- The proved range (r ≤ 1/3), and even the conjectural canonical-family range (r ≤ 1/2), give nothing.
- The needed range r ∈ (0.73, 1] is beyond the r = 1/2 self-duality wall of §14.2.

### 14.4 Large values instead of moments
- Halász–Montgomery with the actual pair kernel K(u,u') = Σ_n W(n/D)χ_n(u/u′) has |K| ≲ U (Pólya–Vinogradov). This needs V > U, but the witnesses have V = U^{δr} < U.
- Burgess at conductor U² saves only when D > U^{3/4}, while the witness lengths are U^{0.62–0.75}.
- So large-values methods give nothing in the range that matters.

### 14.5 Three concrete proof attempts at r ≈ 0.75 and where each stops

**(a) Sixth-power lift of M4(1/2).** If M4 held for r ≤ 1/2 on every family, amplification would give Σ_{u~U}|M_u(U^r)|⁴ ≪ U^{1/6+5r/3}. That beats the second-moment count for r < 0.70 (at δ = 0.475).
- But the record's binding point is the crossing of the inverse and plain counts at t ≈ 1.167, r_* ≈ 0.748, where the lifted count (0.643) does not bind (R stays 0.5849).
- Checked in c61/c61lift.py: no change at any (δ,p). So M4 must hold directly at r ≈ 0.75.

**(b) Canonical recursion with the second Möbius factor as a coherent label.** After Poisson in u, the Gauss identity gives column coefficient αGγ̄₂(N)·τ_W(N). All cross phases cancel exactly, as checked via G(vw) = G(v)G(w)R(v,w) and γ₂(ab) = γ₂(a)γ₂(b)χ_b(a)⁴.
- Splitting N = fn gives Σ_k |Σ_f γ₂(f)χ_f(k)·C(k,f)|², where C(k,f) is exactly the paper's canonical column with label f.
- The paper's lemma only controls Σ_f |C(k,f)|², with the label outside the square. Cauchy–Schwarz over f loses the full √X.
- More decisively, the canonical estimate is the large-sieve shape (K + N)·‖a‖² in the regime of rows shorter than columns. M4 at r > 1/2 is exactly that regime, but needs K·‖a‖², i.e. beating the N-term. No step of the recursion ever beats the N-term.

**(c) Large values.** As in §14.4: Pólya–Vinogradov needs V > U, and Burgess needs D > U^{3/4}.

**What a proof must do.**
- In the t-aspect analogue, M4(3/4) is ∫₀^T |Σ_{n≤T^{3/4}} μ(n) n^{−1/2−it}|⁴ dt ≪ T^{1+ε}.
- Its off-diagonal is the shifted convolution Σ_{h≤X^{1/3}} Σ_{m≤X} (μ∗μ)(m)(μ∗μ)(m+h) with X = T^{3/2}. It needs a power saving. That is a power-saving averaged Chowla estimate with shifts up to X^{1/3}; Matomäki–Radziwiłł–Tao give only logarithmic savings.
- In zero-density terms, the binding count sits at σ = (1+δ)/2 ≈ 0.74. There the record's family bound already matches Guth–Maynard-quality exponents for ζ, and the improvement needed is beyond them.
- So a proof of M4(3/4) for the sextic family would in particular give a new result of this type for a thin family with weaker tools than ζ has.

### 14.6 The spectral route, and why M4 past r = 1/2 is a super-density-hypothesis statement

**Spectral large sieve.** Kuznetsov on the cubic cover (Livné–Patterson, Proskurin) gives the Deshouillers–Iwaniec-type spectral large sieve
Σ_j |Σ_n a_n ρ_j(n)|² ≪ (family size + length)·‖a‖²,
plus Eisenstein terms (Gauss-sum products) and residual (theta) terms. This is the "family + length" shape. At best it removes the BGL loss (UD²)^{2/3} and gives M4 for r ≤ 1/2. Beating the length term needs arithmetic input on μ.

**M4 at r > 1/2 implies bounds stronger than the density hypothesis (explicit-formula heuristic).**
- A simple zero ρ = σ + iγ of L(s,ψ_u) gives M_u(D) ≈ D^{ρ−1/2}Ŵ(ρ)/L′(ρ) + … .
- Rigorously, Plancherel over a range of lengths gives a lower bound for the mean of |M_u|² from the pole of 1/L at ρ.
- So M4 on lengths around U^r forces
  N(σ) ≪ U^{1 − 2r(2σ−1) + ε}.
- At r = 1/2 this is the density hypothesis U^{2(1−σ)}. For r > 1/2 it is strictly stronger.
- Equivalently, under the density hypothesis alone the zeros contribute U^{1 + (2σ−1)(2r−1)} > U, so the density hypothesis does not imply M4(r) for any r > 1/2.
- At r = 3/4 the requirement is N(σ) ≪ U^{2.5−3σ}, e.g. U^{0.25} at σ = 3/4. No family of L-functions has anything like this unconditionally, including ζ: Guth–Maynard gives T^{0.58} at σ = 3/4.

**Conclusion.**
- The targets M4(0.75) → 0.87464 and M4(1) → 13/15 are each stronger than the density hypothesis for the sextic family, which by itself already gives 13/15 (§13.4).
- A spectral (large-sieve-shape) method gives at most r ≤ 1/2, which changes nothing (§14.3, §14.5a).
- So the spectral route cannot reach the needed range unless it carries zero-density-strength arithmetic information.

## 15. A small improvement below the record: θ = 0.87485 (c = 0.12515), conditional on one extension

Full write-up: c61/IMPROVEMENT.md. Scripts: c61/fastR.py, cert.py, sidechecks.py, thopt2.py.

**The slack.**
- Cycle61's detector (eq. 3.10) counts plain witnesses with the old capacity (1−2m)/(6κ).
- At the binding class (δ ≈ 0.466, p = 1/2) the witness length is m ≈ 0.417. That lies inside the window where Cycle61's own interior affine moment (Thm 2.2: n₁+n₂+Λz ≤ B) applies, so the capacity rises to (B − 2m)/Λ ≈ 0.0384.
- The window is m ∈ [0.4137, 0.4235]. Its upper end keeps any shorter selected slot set paid by the old moment.

**Effect.**
- At the binding, R drops by about 1.1·10⁻⁴ and the row margin at the record geometry grows from 5·10⁻⁶ to 7·10⁻⁵.
- Re-optimising h (with the cap M + ℓ = 1 kept) gives a model optimum θ ≈ 0.8748446.

**Certified (exact rationals).**
- θ = 87485/100000, h = 79577/100000.
- Rows d = h, d₋, ½ certified with Cycle61's reserves 10⁻⁶, 10⁻⁵, 10⁻³.
- All side conditions (3.3), (3.5), (3.8), (4.5)–(4.7), (5.6), (5.14), (5.15) hold.
- κ-feedback: |∂G/∂κ| < 0.16 < 1/2.

**Status: conditional.** In Cycle61's Lean source, `affine_fourier` (B064) is stated only for the physical numerator pair (η = fixedSourcePrincipal). The detector uses the κ-marked moment for every dictionary character. So the improvement needs Theorem 2.2 for the detector's dictionary characters, an extension of Cycle61's proof whose inputs are stated for arbitrary fixed finite-order characters. Not formalized, not submitted.

## 16. Session 3 attempt at 13/15 (sextic/SESSION3_1315.md)

- **Lemma B (PROVED).** Inside the contradiction argument (no zeros right of β*), harmonic majorant + Jensen give
  |L(σ+it,ψ)|^{q₀/(1−λ)} ≤ q^{O(ε)}·∫P|L(½+iv,ψ)|^{q₀}, σ = ½ + λ(β*+ε−½). With the family fourth moment this yields
  Σ_u |L(σ+it, ηχ_u)|^{4/(1−λ)} ≪ U^{1+ε}: every moment is Lindelöf-on-average strictly right of ½.
- **Proposition B2 (PROVED, exponent algebra).** These moments never improve a detector count: for a zero at
  b = β−½ ≤ b₁ the savings 4m(b−λb₁)/(1−λ) (plain) and 4(b−λb₁)/(3−λ) (mixed) decrease in λ, and the
  Halász–Montgomery off-diagonal would need 2β−1 > β*+ε. The self-bootstrap route is closed.
- **Quartic alternative.** Reaching 13/15 needs M > ½ in (G), which is inside Theorem ZF's quasi-RH regime.
- **Status.** 13/15 OPEN; the missing input is still DH-strength zero density for the sextic family on σ ∈ [0.51, 0.867].
