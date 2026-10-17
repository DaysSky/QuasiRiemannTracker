# The density hypothesis for the sextic twist family and θ = 13/15: what is proved, what is needed

DaysSky, with Claude assistance. 16 October 2026.

Labels: **PROVED** (a complete argument inside the stated framework); **MODEL** (exponent algebra on the record's constraint
system, evaluated numerically); **NUMERICAL** (computer experiment); **CITED** (literature); **CONJECTURE**.
Details, scripts and intermediate reasoning: `sextic/DH_SEXTIC.md` §1–§63 and `sextic/DH_SUMMARY.md`.

## 0. Summary

- Within the 7/8 architecture, θ = 13/15 follows from a density-hypothesis count for the sextic twist family on
  σ ∈ (1/2, 0.875] at the row scales d ∈ [2/5, h+ζ], together with lowering the paper's bin floor a0 = 51/100 towards 1/2
  (DH_SEXTIC §50–§52). The floor change and the needed extension of the principal Euler region below 0.87 are margin
  changes, checked against the exact local Euler tables. With the paper's floor, DH on [0.51, 0.875] gives θ = 0.86979.
  Every margin of the paper's §20 has been re-evaluated exactly at the new geometry (DH_SEXTIC §61, `DH_THEOREM.md`). Below
  d = 2/5 the paper's own count suffices. At the binding class it needs a DH-level count. At the record's binding class (σ = 0.733) the current count is 0.5935 and DH asks for 0.534.
- **Theorem A:** with every available mean-value input, the count is exactly R(δ, p). An explicit adversary row realises it.
- **Theorems B and C:** the two classical routes beyond the toolkit, Guth–Maynard-type large values and Poisson for the long
  Möbius mean square, are exact involutions for this family. They cannot gain without a third ingredient.
- **Proposition C′:** the decisive input LM is the family's weak-Mertens problem. Its ζ-analogue (power-saving variance of
  μ in short intervals) is open. **Theorem D:** combinatorial identities for μ cannot replace it.
- **Numerics:** the long Möbius mean square shows no off-diagonal mass up to 76× past the effective conductor (U = 10³).
  Against exact coprime diagonals, the off-diagonal is within 1.5% at the effective conductor for U = 10⁴, 3·10⁴ (all rows)
  and 10⁵, consistent with LM (and with GRH-on-average, which no finite test can separate). The true marked and plain-slot moments show no excess beyond the
  proven ranges, for any slot coefficients including Gauss-correlated ones. This is consistent with natural slot costs.
  Both are supporting evidence, not asymptotic tests (§9 caveat).
- **MODEL:** LM at every scale gives exactly 13/15 (with the floor lowered). Natural slot costs alone give ≈ 0.8698.

## 1. Setting

Rows u: sixth-power-free elements of O = Z[ω] with N(u) ≍ U. Characters ψ_u(n) = ν(n)χ_n(u)^{±1}, ν in a fixed finite
group Θ. Polynomials, normalised by length^{−1/2}:
- Möbius witness M_u(U^r) = Σ_n μ(n)ψ_u(n)W(N(n)/U^r);
- plain witness S_u(U^m) = Σ_n ψ_u(n)W(N(n)/U^m);
- prime slots P_u(U^z) = Σ_p a(p)ψ_u(p)W(N(p)/U^z).

A row is in bin (δ, p) if its maximal zero (bounded height) has real part (1+δ)/2 and its prime slots carry the fraction
p ∈ [0, 1/2] of their in-bin cap. **Conditional Theorem D** (REPORT §3): a count U^{1+λ−δ} in every bin gives
θ = 11/12 − (1−2λ)/(4(5+6λ)), and λ → 0 gives 13/15. "Every bin" includes the floor. The paper's floor a0 = 51/100
has the trivial count R = 1, i.e. λ = δ0 = 1/50, which caps θ at 0.86979. So λ → 0 also needs a0 → 1/2. The floor
is the holomorphy boundary sr ≥ 51/100 of the Euler-correction region D1. The exact local Euler table allows a0 = 1/2 + ε:
the binding exponent 4 − 6sr − 6zr < −1 holds for sr > 0.4933 at zr = 17/50 (§52). Every occurrence of 51/100 in the paper
(Lemma 7.1, Props. 8.3, 9.2 and 16.1, Lemma 19.1, the floor bound (20.5)) was checked to be a margin, so any a0 > 1/2 works (§53). DH is needed at the top row scales. Below them the record's own
counts suffice: with Liu's cruder count R = 76/75 − 2δ/3 at d ≤ 1/2, θ would stop at 0.8701 (§51).

**Conditional Theorem D′ (corrected statement; PROVED relative to the OpenAI paper, DH_SEXTIC §49–§53, §61).**
*Hypothesis (DH_λ).* Take the floor a0 = 1/2 + λ/8. For every fixed bin width e > 0, every bin a ∈ (a0, β*] and every row
scale U = Z^d with d ∈ [2/5, h + ζ], the number of sixth-power-free rows in an (i, a) bin and prime-amplitude class is
≪ U^{1+λ−δ+ε}(1+T1)^A, with δ = 2a − 1.
*Conclusion.* Every finite-order Hecke L-function over Q(√−3), and every Dirichlet L-function, is zero-free on
Re s > 11/12 − (1 − 2λ)/(4(5 + 6λ)) + O(λ). In particular DH_λ for all λ > 0 gives the zero-free half-plane Re s > 13/15.
*Proof outline.* The OpenAI argument is unchanged except for the following points.
- Paper §20.2 (Lemma 10.4 exponent): E(d) = a − B + h(17/50 − 1/6) − a·l_y − ℓ/2 + qℓ + d(R + δ/2 − 17/50). At d = h, with
  R ≤ 1 + λ − δ and q ≤ δ/2, this is ≤ C0 + hλ + δ(x + ℓ/2 − 1/2) (REPORT §3). For d ≤ 2/5 the paper's own no-prime count
  76/75 − 2δ/3 (§20.5, t = 1) suffices: E(2/5) ≤ −0.0142, and it stays below (1+δ)/2 − σ0 for bins δ ≤ 5/6 (DH_SEXTIC §61,
  exact rational check). That count fails only on (0.426, 1/2] for δ near 2σ0 − 1 (§51). There the hypothesis applies, or
  alternatively the detected count of §2, whose λ ≤ 0.084 lies inside the tolerance (h/d − 1)(0.66 − δ/2) ≥ 0.17. The
  detected count would, however, need the slot mesh redone below d = 1/2.
- The floor bin a0 contributes C0 + (3/4)δ0 with δ0 = 2a0 − 1 ≤ λ (paper (20.5)). Lowering a0 from 51/100 is legitimate:
  every use of 51/100 (Lemma 7.1 region (7.14), Props. 8.3, 9.2 and 16.1, Lemma 19.1) is a margin (DH_SEXTIC §52–§53).
- The principal correction H_η is holomorphic and close to 1 for Re s > 2/3 at the residue point (w, z) = (1, 1/6), so the
  continuation criterion (OpenAI Prop. 2.1 "Continuation from a common signal"; Xu–Sang Lemma 8.1) applies with B < 0.87 (§52).
  Accordingly the D2 edge 7/8 moves to 0.86, and Lemma 10.5 is used with σ0 = B.
- The small-row lines (β* + e, 1/2, 17/50) lie in D1(1/3) instead of D1(3/8). This is a choice of ε0, which the paper
  allows for every ε0 > 0. The exponents in the proof of (16.15) stay
  ≤ 0 (selected primes) and < −1 (unselected), so the small rows keep the margin −0.0417 (§61).
- Every other margin of paper §20 was re-evaluated exactly at these parameters: principal m_w = l_y/20, m_z = h/600, prime
  supply, low-side excess and Gram, large rows. All remain positive (DH_SEXTIC §61, margins/part2_margins.py).
Optimising ℓ < (1 − 2λ)/(5 + 6λ) and x as in REPORT §3 gives the bound. Numerically this reproduces
θ(λ) = 0.866683, 0.869792, 0.872774 for λ ≈ 0, 0.02, 0.04 (floor/scales.py). ∎

**In-bin caps (PROVED; paper Lemma buffered-bins).** |L|, |1/L| ≪ U^ε on Re s ≥ a + 6e. Hence normalised excesses
(|X|² = U^{len+exc}) satisfy:
- exc ≤ δ·len for Möbius and prime polynomials of any length;
- exc ≤ δ·min(c, 1−c) for plain polynomials of length c ≤ 1 (Hecke reflection);
- plain polynomials longer than 1 are O(U^{−A}).

## 2. The count and the binding (MODEL; sextic/bindpt.py)

Witness costs at bin (δ, p), with q = δp:
- A_I(r) = 1 − q − (δ − q)r for r ≤ 1 (marked inverse moment, slot cost c_I = 2);
- A_P(m) = 1 − 2δm − 2q·z_P(m) (plain fourth moment, slot capacity z_P from cost 6κ or Cycle61 Thm 2.2);
- long(t) = (1+5t)/6 − δt for t > 1 (sixth-power copies).

The detector holds for every t ∈ [1, 3/2]. The count is R(δ, p) = min_t max(short(t), long(t)), where
short(t) = max_{r+m=t} min(A_I(r), A_P(m)). At the record (θ = 0.8748389) the binding class is δ = 0.466, p = 1/2, with
t* = 1.162, r* = 0.745, m* = 0.417. There A_I = A_P = long = 0.5935, against the DH value 0.534.

**Closed form and the DH gap (DH_SEXTIC §49).** The count above equals Xu–Sang's closed form (arXiv 2610.15291, (21)) to
machine precision: R = 1 − δ + λ with λ = (α − δ)δP(x)/(2J), where P = (2 − 2cx)(1 − x), J = (α−δ)(3 − (1+2c)x) + δP,
c = 1/(3κ) and α = 5/6 is the long (LM) slope. The DH gap is proportional to α − δ. A long slope α′ therefore gives DH for
the family on σ ≥ (1 + α′)/2, and nothing below that unless the plain slot cost is natural and x = 1/2. At the binding
class λ = 0.0596.

## 3. Theorem A (toolkit closure; PROVED within the exponent algebra of the record's constraint system, sextic/closure.py)

**Inputs considered:**
- marked inverse moment (r + 2z ≤ 1); plain lemma (n1 + n2 + 6κz ≤ 1, or Thm 2.2);
- sixth-power copies for any polynomial of length > 1; Möbius fourth moment via copies;
- marked Möbius×Möbius (r1 + 2r2 ≤ 1); theta-reflection MM′ (r + 2m < 1);
- pure prime-slot moments in the BGL range; the pointwise plain sixth moment;
- all Hölder combinations, and the detector at every t ∈ [1, 3/2].

**Statement.** The best bad-row count obtainable from these inputs is exactly R(δ, p).

*Proof.* Upper bound: by construction. Lower bound: exhibit a row type, the adversary, consistent with every input
whose cheapest witness costs R. The detector at all t is equivalent to the covering condition
[1, 3/2] ⊆ ∪_{r∈𝒟}[r, r + m_max], where 𝒟 is the set of saturated Möbius lengths and m_max = max 𝒩 (§12).
- A cover avoiding r = 1 (which costs 1 − δ) needs consecutive points p_a ≤ 1 < p_b with p_b − p_a ≤ m_max.
- The row then costs min(A_I(p_a), A_P(m_max), long(p_b)). Maximising gives max_s min(short(s), long(s)), which equals R
  because short decreases and long increases.
- The adversary 𝒟 = {r*, t*}, 𝒩 = [0, m*] attains it.
- Every other listed input evaluated on this row is ≥ R. At δ = .466: M4-copies 0.714, (r*, t*)-product copies 0.867,
  pure slots 0.767, plain-6th pointwise 0.792; MM′ and marked M×M are inapplicable. Checked at δ = .3, .4, .466, .55, .65.
- Hölder never goes below its vertices. ∎

**Robustness (PROVED, §24, §27).** The adversary extends consistently to every constraint coming from:
- the exact identities μ∗1 = δ and μ·log = −μ∗Λ: at each length the number of large pieces is never exactly one;
- Vaughan and Heath-Brown decompositions combined with the in-bin caps.
These only force descending chains of saturated Möbius witnesses bridged by saturated plain witnesses. The extension is
𝒟 = {0, 0.328, r*, t*, r″ ∈ (1.5, 1.579]}, and every added witness costs more than R (e.g. A_I(0.328) = 0.69,
long(1.55) = 0.74).

## 4. Theorem B (exact self-duality; PROVED, §9)

Let M_B = (ψ_u(n)W(N(n)/N))_{u∈B, n}. Poisson in n modulo uu′ and the factorisation
τ(ψ_{u1}ψ̄_{u2}) = ψ_{u1}(u2)ψ̄_{u2}(u1)τ(ψ_{u1})τ(ψ̄_{u2}) (coprime moduli) give, after splitting B into boundedly
many classes,
  M_B M_B^* = N·I + (N/U)·D(M′_B M′_B^* − H·I)D^*,  H = U²/N,  D = diag(normalised Gauss sums).
The cubic trace tr((M_BM_B^*)³) is likewise the same trace at dual length H. Guth–Maynard (ζ) and Chen–Gupta–Li (all
characters mod q) gain from a free archimedean or residue variable that yields a linear constraint on dual frequencies.
Here the moduli u1u2, u2u3, u3u1 differ and the Gauss sums cancel cyclically, so no such variable exists. Large-values
estimates at length N and at length U²/N are the same problem.

## 5. Theorem C (the long slope is a fixed point; PROVED, §22, §26) and Proposition C′ (weak Mertens; PROVED, §34)

**Theorem C.** For U < D (all elements u, smooth weight),
  MS(U;D) := Σ_u w(u/U)|M_u(D)|² = UD·ŵ(0) + (U/D)·Σ_{h≠0} ŵ(hU/D²)(|B_h(D)|² − D),
with B_h(D) = Σ_n μ(n)g̃(n)χ̄_n(h)W(n/D) and K = D²/U. Applying Poisson in h returns MS(U;D) with coefficient exactly 1.
By Hasse–Davenport (checked numerically at 2248 primes, DH_SEXTIC §56), μ(n)g̃(n) = ε(n)(n/|n|)conj(g̃3(n)) with ε of finite
order. So the dual coefficients are Grössencharacter-twisted cubic Gauss sums, and spectral tools would live on the
cubic cover, where Patterson's theory is explicit. The self-duality is unaffected.
- Copies (u ↦ ua⁶, P = (H/U)^{1/6}) give MS ≪ H/P at H = D^{1+o(1)}, i.e. slope 5/6. On the dual side this is a relative
  saving (K/D)^{−1/6}. LM (α′ = 0) needs relative saving D/K.
- Copies map a slope α′ to min(α′, 5/6), so 5/6 is the fixed point from slope 1.
- In OpenAI's canonical recursion, the condition r < 1 enters through the invariant F − M − Q − z0 > 0, which equals
  m − r − 2z (the puncture cancels identically). At r > m the child is rows-dominant and only the diagonal survives.

**Proposition C′.** Shift Perron's integral to Re s = −ε, where |1/L(−ε+it)| ≪ U^{−1/2}log U by the functional equation.
This gives
  M_u(D) = Σ_ρ W̃(ρ)D^ρ/L′(ρ, ψ_u) + O(D^{−ε}U^{−1/2}log U).
So LM requires (i) DH-type control of zeros off the line and (ii) the family average of Σ_ρ |W̃(ρ)|²/|L′(ρ)|²
(Gonek–Hejhal J_{−1}) to be ≪ U^{1+ε}. For ζ the analogue is the weak Mertens conjecture. In the t-aspect it is the
power-saving variance of μ in short intervals of relative length 1/T, where only (log)^{−A} (Ramachandra 1976) and
o(H²) (Matomäki–Radziwiłł 2016) are known [CITED; Gorodetsky, Oxford thesis §1]. The quadratic-character analogue has an
asymptotic only below the conductor, and only under GRH [Gao–Zhao, arXiv 2212.09241].

## 6. Theorem D (decompositions cannot prove LM; PROVED, §24, §30a)

Each of the standard combinatorial identities for μ·1_{[D,2D]} (Vaughan, Heath-Brown k-fold, Linnik) contains an explicit
Möbius-heavy piece with no smooth factor of length ≥ U^c, e.g. (μ_{≤y})^{∗k} with y^k ≍ D. Heuristically this is
forced for any such identity, since 1/ζ has no long ζ-factor. For such a piece, LM is the 2k-th moment of a short Möbius
polynomial at total length D > U, which is again beyond the conductor.
- The pieces that do have long plain factors are fine: plain factors of length > 1 vanish, and two plain factors summing
  to t > 1 cannot be saturated by the caps.
- The three-plain pieces give count 1 − 2tδ/3, which beats copies only for t > 5/(5−2δ) (1.229 at δ = .466), not at t*.

## 7. Proposition E (cost-1 slots from the Möbius fourth moment; PROVED, §33)

Let P_u be one prime slot of length Z = U^z with z < 1/4. Sum over all u ≤ U with a smooth weight (positivity). Then
  Σ_u |M_u(N)P_u|² ≤ Σ_p|a(p)|²·Σ_u|M_u(N)|² + Z·(U·Σ_u|M_u(N)|⁴)^{1/2}·U^ε.

*Proof.* The off-diagonal is Σ_{p1≠p2} a(p1)ā(p2) Σ_u χ_{p1}χ̄_{p2}(u)|M_u|². Apply Cauchy over the ≤ Z² pairs and
large-sieve duality for 𝒮 = {χ_{p1}χ̄_{p2}}. Distinct members of 𝒮 differ by a nonprincipal character of modulus
≤ Z⁴ ≤ U^{1−ε}, whose smooth lattice sum is O(U^{−A}). So the Gram matrix of 𝒮 is U·I + O(U^{−A}). ∎

A sharp Möbius fourth moment therefore gives slot cost 1 for z < 1/4. At r* this beats the record only if the
fourth-moment loss is < 0.114; copies give 0.408.

## 8. Where the slot costs and the slope come from (PROVED readings of the paper, §32, §37)

- **c_I = 2.** Active slot primes sit in the theta-reflection conductor. The maximal dual length T_d = 2H + 2A0 + 2z_a − N*
  grows by 2z_a, while the saving term u = min{v, z_a, (v+z_a)/3} is capped at z_a, which nets +z_a. A secondary BGL penalty
  (2z_a − v)/3 appears only on small-dual dyads. So c_I = 2 is intrinsic to detecting slots by theta reflection.
- **c_P = 6κ and α′ = 5/6.** Both come from sixth-power copies, whose efficiency is 1/6: the plain lemma finances
  pointwise slot removal by averaging h ↦ hp⁶, with row extension 6ℓ per allowance ℓ. Sixth powers are the only
  multipliers fixing every column Gauss sum. A replication mechanism with efficiency e would give α′ = 1 − e and c_P = 1/e.

## 9. Numerical evidence (NUMERICAL; sextic/mmnum/)

**LM** (lm.py; lm_big.py with the table-free kernel mmk2.c). Ratio of mean|M_u(U^s)|² to the diagonal. Rows are all
sixth-power-free N(u) ≤ U, or a random sample where the size is given in parentheses:

| U | s = 1.0 | 1.1 | 1.2 | 1.3 | 1.4 | 1.5 | cols_eff/rows at s_max |
|---|---|---|---|---|---|---|---|
| 3·10³ | 0.946 | 0.938 | 0.926 | 0.931 | 0.953 | 0.928 | 1.05 |
| 10⁴ (6000) | 0.919 | 0.952 | 0.950 | 0.930 | 0.924 | 0.930 | 1.92 |
| 3·10⁴ (4000) | 0.938 | 0.923 | 0.934 | 0.953 | 0.913 | – | 1.18 |
| 10⁵ (4000) | 0.934 | 0.937 | 0.914 | – | – | – | 0.19 |
| 3·10⁵ (4000) | 0.927 | 0.941 | – | – | – | – | 0.07 |
| 10⁶ (4000) | 0.921 | 0.944 | – | – | – | – | 0.08 |

Standard errors are ≈ 0.015. The random multiplicative control is identical, and max|M|² stays ≲ 19·mean.
**Caveat (DH_SEXTIC §44).** Rows number ≈ 3.62U, while the effective column count is cols_eff. Even a sharp large sieve
allows a ratio up to 1 + cols_eff/rows. So only the points with cols_eff/rows ≈ 1–2 genuinely test LM: U = 10⁴, s = 1.5;
U = 3·10⁴, s = 1.4; U = 3·10³, s = 1.5. Those allow ×2–3 and show ×0.91–0.93. The other points are rows-dominant. The
support for LM is real but modest. Likewise the slot tests below are in the rows-dominant regime for M·P. They exclude
hidden coherence (excess) but do not by themselves separate cost 1 from cost 2 asymptotically.
**Far beyond the conductor (DH_SEXTIC §45, lm_beyond.py).** At U = 10³ (all 3636 rows) and U = 3·10³ (3000 sampled rows),
going to cols_eff/rows = 76 and 57 respectively (sharp-sieve room ×77 and ×58), μ stays at 0.88–0.96 of the diagonal, as does
the random multiplicative control. The squarefree indicator, once the trivial row is excluded, falls to 0.14 because its
sums are mostly smooth. So the test separates coefficient types. The μ plateau is what Proposition C′ predicts under GRH
with Gonek–Hejhal-size J_{−1}. It supports LM, but no finite computation can separate LM from GRH-on-average.
- LM for 1/L² (coefficients μ∗μ, lm2.py) is equally flat: 0.84–0.91 up to s = 1.5.
- On the dual side, B_h(D) = Σμ(n)g̃(n)χ_n(h)W is uniform in h (mean 0.92–0.95 on every range of N(h)), with no
  small-h bias (dualB.py).

**Cross-scale tails** (xscale.py, U = 3·10⁴). Rows in the top 1% at both r* = 0.745 and t* = 1.162 are 7.7× more common than
independence predicts. 93% of them are above the median at r = 1. The adversary's "dip at r = 1" is atypical. By the
explicit formula plus Nazarov–Turán (DH_SEXTIC §40), it requires ≍ log U zeros clustered near the bin, modulo a
far-zero remainder whose control needs high moments.

**Slot costs** (slotnum.py, slotplain.py). ρ = mean|X·P|²/(mean|X|²·mean|P|²), with slot coefficients 1 / conj(g̃(p)) /
random, far outside the proven ranges:

| X, (lengths) | U = 10⁵ | U = 10⁶ | U = 3·10⁶ |
|---|---|---|---|
| M(U^0.5), z = 0.40 (r+2z = 1.3) | 1.001 / 0.998 / 1.000 | 0.999 / 1.000 / 1.000 | 1.0005 / 0.9993 / 0.9997 |
| M(U^0.5), z = 0.45 (r+2z = 1.4) | 0.999 / 1.004 / 1.001 | 1.000 / 1.000 / 1.001 | 1.0000 / 1.0000 / 1.0001 |
| S(U^0.35)², z = 0.30 (2m+4.5z = 2.05) | 1.000 / 0.998 / 0.998 | 0.999 / 1.002 / 1.000 | 0.9992 / 0.9994 / 0.9995 |

There is no drift from 10⁵ to 3·10⁶: a power loss U^c would drift by 30^c, so c ≲ 10⁻³.

**At the binding Möbius length** (slotbig.py, r* = 0.745, 20 000 sampled rows, U = 10⁷, 10⁸, 10⁹; slots with
r + 2z = 1.145 and 1.245):
- random-phase slots give ρ = 0.974–1.009 (±0.012);
- a ≡ 1 gives a stable 1.04–1.07 that does not grow over a factor 100 in U. This is exactly the Möbius pair term
  n1 = kp2, n2 = kp1 (μ(kp1)μ(kp2) = μ(k)²), of relative size ≈ 1/(2 log Z), i.e. diagonal mass of the decoupled sequence
  b(m) = −μ(m)ω_a(m). No power-type penalty appears.

**Other checks:**
- Sextic Möbius–Gauss sums (gbias.py, twisted multiplicativity verified by brute force): Σ_{N(n)≤X}μ(n)g̃(n) = −43.9,
  +20.9, +72.5 at X = 1875, 7500, 30000. These are square-root sized with sign changes; no Patterson-type bias.
- Möbius fourth moment, plain sixth moment and large-value tails all match the random multiplicative model
  (DH_SEXTIC §4, SESSION3 §6.1, §6.7).

**Formulation (DH_SEXTIC §46, PROVED).** LM must be stated for sixth-power-free rows. Over all elements, the rows εa⁶ are
V^{1/6} copies of the unit-type rows ψ_ε (conductor O(1)). An all-elements slope α′ would then give a zero-free region
Re s > 1/2 + (5 + 6α′(s−1))/(12s) for L(s, ψ_ε): 11/12 at α′ = 5/6, which is the architecture's base exponent, and 7/9 at
α′ = 0. The sixth-power-free LM never reaches D → ∞ for a fixed row, so it is strictly weaker, and it is what the bins need.

**The dip as a divided difference (DH_SEXTIC §47).** A cluster of zeros contributes F[ρ_1, …, ρ_n] with F = W̃e^{(s−a)ξ}/g
(exact). If W̃/g has nearly constant phase on the cluster, Hermite–Genocchi makes this ξ^{n−1}× a Laplace transform of a
positive measure. A dip of the adversary's depth then needs relevant zeros spread over heights ≳ 0.2, not a local cluster.
Equally spaced progressions fail at the saturated ends. Closing the route still needs a multi-zero density bound plus
control of L′/L, so it is a sharper form of §40, not a proof.

**Exact off-diagonal (DH_SEXTIC §57).** With each row's exact coprime diagonal (kernel mmk3.c), MS/truediag is the exact
relative size of the off-diagonal. Over all rows at U = 3·10³, it is 0.986–1.006 for s ∈ [0.8, 2.0], up to 57× past the effective
conductor. At the effective conductor it is 0.996 ± 0.006 (U = 10⁴, all rows) and 1.0013 ± 0.0034 (U = 3·10⁴, all 108 606 rows). This is the same as
the random multiplicative control, with no term of the diagonal's order.

**Function-field analogue (DH_SEXTIC §58).** Over F_q[t], M_u(k) = q^{k/2}h_k(Frobenius). With big monodromy, Deligne–Katz
equidistribution gives LM with constant 1 at every length as q → ∞, because Sym^k is irreducible (CITED mechanism, as in
Keating–Rudnick). At accessible degrees (d ≤ 6) the sextic Kummer family has finite or CM monodromy, and its mean square
is driven by rows with repeated Frobenius eigenvalues. So it is not a proxy for the number-field question.

## 10. θ under combined inputs (MODEL; c61/capc4.py, Cycle61 grid, baseline 0.87492 ≈ record 0.8748389)

| slot costs \ LM slope α′ | 5/6 (now) | 0.70 | 0.50 | 0.30 | 0 |
|---|---|---|---|---|---|
| c_Pκ = 4.5, c_I = 2 (now) | 0.87492 | 0.87372 | 0.87055 | 0.86844 | 0.86667 |
| c_Pκ = 2, c_I = 2 | 0.87339 | 0.87267 | 0.87040 | 0.86847 | 0.86667 |
| c_Pκ = 4.5, c_I = 1 | 0.87289 | 0.87150 | 0.86922 | 0.86793 | 0.86667 |
| c_Pκ = 2, c_I = 1 | 0.87164 | 0.87082 | 0.86904 | 0.86780 | 0.86667 |
| c_Pκ = 1, c_I = 1 (natural) | 0.86979 | | | | 0.86667 |

**Literature, 16 Oct 2026 (DH_SEXTIC §48).** Xu–Sang (2610.15291, 0.8748867) and Liu (2610.12234, 0.874957) introduce no
DH or LM input. Xu–Sang's averaged-numerator high side is the one already in our Cycle61-type model. Xu's large-sieve
lower bound (described there) shows the general-coefficient sextic large sieve is not sharp for D ∈ [U^{16/11}, U²]. So
even the existing slope 5/6 depends on the arithmetic of μ, and so must any improvement.

## 11. The frontier, ranked

| input | gives | status | classical analogue |
|---|---|---|---|
| LM: Σ_u\|M_u(U^s)\|² ≪ U^{1+ε+α′(s−1)}, α′ → 0, s ≤ 3/2 | 13/15 (floor → 1/2) | CONJECTURE, numerically supported (§9, §57) | weak Mertens; μ in short intervals with power saving (open for ζ) |
| LM only near the conductor, s ∈ [1, s0] (§55) | 0.8732 (s0 = 1.05), 0.8714 (1.1), 0.8695 (1.16); DH on σ > 0.781 for s0 = 1.16; slope 0.7 on [1, 1.05] alone gives 0.87467; a U^{0.003} saving over copies on s ∈ [1.10, 1.25] passes the record (MODEL) | CONJECTURE; smallest high-value form (§54) | power-saving variance of μ in intervals of length X^{small} (open) |
| natural slot costs c_I = 1, c_Pκ = 1 | 0.8698 | CONJECTURE, numerically supported (rows-dominant range) | critical-length mean values with prime factors; needs a non-reflection slot treatment |
| Möbius 4th moment at r ≈ 0.75 with loss < 0.11 | partial c_I | open | LM for 1/L² at length 1.5 |
| MLV / MM′ / plain sixth moment at supercritical length | 13/15 jointly | open | Lindelöf-on-average sixth moments (open for quadratic characters) |

Every route to 13/15 inside this architecture passes through a statement of LM type: Möbius randomness against the family
beyond its conductor. Within the theta-reflection/large-sieve toolkit, this statement is a fixed point (Theorems A–D). All
experiments so far are consistent with it (up to 76× past the effective conductor at U = 10³). The sextic family's current slope 5/6 already beats the ζ analogue (slope 1), so any α′ < 5/6
would be a new kind of power-saving Möbius-variance theorem.

## 12. Reproducibility
- `sextic/bindpt.py`, `closure.py`: binding point and closure check.
- `c61/capc3.py`, `capc4.py`: θ tables.
- `sextic/mmnum/lm.py`, `lm_big.py`, `lm_beyond.py`, `slotnum.py`, `slotplain.py`, `slotbig.py`, `gbias.py`,
  `gbias_check.py`, `dualB.py`: numerics. The table-free kernel is `mmk2.c`.
- `sextic/floor/floor.py`, `scales.py`: the floor and row-scale corrections (§50–§51). `sextic/dip/prog.py`: dip model (§47).
- `c61/lmrange.py`: partial-LM θ and DH-range tables (§55). `sextic/mmnum/mmk3.c`, `lm_exact.py`: exact coprime diagonals (§57).
  `hdcheck.py`, `hdcomp.py`: Hasse–Davenport and cubic Gauss sum identities (§56).
- `sextic/ff/ff_lm.py`, `brute.py`: function-field analogue (§58).
- Build the kernels with `gcc -O3 -march=native -fopenmp mmk3.c -o mmk3 -lm` (likewise `mmk.c`, `mmk2.c`).
- Paper references: OpenAI, "The Quasi-Riemann Hypothesis" (30 Sep 2026). Line numbers refer to paper.tex or, where stated,
  to the public PDF's text (scratchpad/oai_pdf/oai5.txt).
