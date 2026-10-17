# Pushing the Cycle61 boundary below 874861/1000000

**Claim (relative to Cycle61's pinned inputs).** Every finite-order Hecke L-function of trivial infinity type over Q(√−3), and every Dirichlet L-function, is nonzero for Re s > θ = 87485/100000 = 0.87485 (c = 0.12515), apart from principal poles. The record is c = 0.125139.

The proof is the Cycle61 proof with **one change**: the physical detector's plain-witness count uses Cycle61's own Theorem 2.2 wherever it applies, instead of the old paired moment (1.6).

## 1. The change

**Cycle61's count.** Cycle61 §3, eq. (3.10), counts rows with a saturated plain witness |P(U^m)|² ≥ U^{δm} through the paired plain moment |P P Q_I|², with the old slot capacity
z_P^{old}(m) = (1 − 2m)/(6κ),  giving  A_P = 1 − 2δm − 2q z_P(m),  q = δp.

**What Theorem 2.2 allows.** Theorem 2.2 (interior affine moment, same canonical positive family as (1.6)) gives the paired moment ≪ U^{1+ε} whenever
.41 ≤ n₁, n₂ ≤ .44,   .034 ≤ z_I ≤ .04,   n₁ + n₂ + Λ z_I ≤ B,
with Λ = 19312/4787 and B = 14200/14361 − 10⁻⁵.

**The new capacity.** Take n₁ = n₂ = m. On the window
m ∈ W = [m_lo, m_hi] = [(B − Λ/25)/2, min(.44, (B − .034Λ)/2, (1 − .034·9/2)/2)] ≈ [0.41370, 0.42350],
use z_P(m) = max(z_P^{old}(m), (B − 2m)/Λ). Off W, keep z_P^{old}.

**Every selected slot set is paid.** The whole-slot selector (3.4) returns I with z_I ≤ z_A = z_P(m) and gain ≥ q·z_A.
- If z_I ≥ .034, Theorem 2.2 applies, since 2m + Λ z_I ≤ B and z_I ≤ .04 on W.
- If z_I < .034, the old moment (1.6) applies, since 2m + 6κ z_I ≤ 2m_hi + (9/2)(.034) ≤ 1 for every κ ≤ 3/4. The bound m_hi ≤ (1 − .153)/2 is chosen exactly for this.
- If I is empty, the unrestricted no-slot moment applies.

All capacities stay below the supply ℓ/d₊ (checked).

**The new count.** R(δ,p) = min over t ∈ [1, 3/2] of max(short(t), long(t)), where
- short(t) = max over r ∈ [t − ½, min(t,1)] of min(A_I(r), A_P(t − r)),
- A_I = 1 − q − (δ − q)r,
- long(t) = 1 − δ + (5/6 − δ)(t − 1).

With z_P^{old} this reproduces Cycle61's closed form (3.9) to 10⁻⁹. With the upgrade it is smaller. At the binding class δ ≈ 0.466, p = ½ the witness sits at m ≈ 0.417 ∈ W with capacity 0.0384, and R drops by about 1.1·10⁻⁴.

## 2. Geometry and certificate

**Geometry.**
- θ = 87485/100000, h = 79577/100000.
- ℓ = 11/3 − 4θ (the low side is tight with M + ℓ = 1).
- x = 1 + ℓ − h, y = h − 2ℓ.
- κ_* = 2θ − 1.

**Row inequalities** (Cycle61 eq. 3.14, with the new R), all checked in exact rational arithmetic (cert.py). The certificate uses adaptive box subdivision of [0, 3/4] × [0, 1/2]. It relies on two facts: R is nonincreasing in δ and in p (every A_I, A_P and long(t) is), and the remaining part of G is bilinear in (δ, p).

| Row | Reserve | Result |
|---|---|---|
| new, d = h | 10⁻⁶ | certified, 112,253 boxes |
| old, d = d₋ = h − 10⁻³ | 10⁻⁵ | certified, 29,547 boxes |
| old, d = ½ | 10⁻³ | certified, 29 boxes |

The structural steps are as in Cycle61: affine in d on [d₋, h] and [½, d₋], junction (4.5), convexity below ½ with reserve (4.6), floors (4.7), and the 10⁻¹⁰ extension to d₊.

**Side conditions** (sidechecks.py, all exact):
- (3.3): n ∈ (0.41966, 0.42139), z_A ∈ (0.03619, 0.03705), 41/100 + N + (9/2)z_A,max < 1.
- (3.5) supply.
- (3.8) S_b(t₀) = 5.37·10⁻⁵ > 0. Here b = 0.16068, slightly above Cycle61's convenience bound 4/25, but S_b(t₀) is checked directly.
- (4.5) junction.
- (4.6), (4.7) reserves.
- (5.6) low side: y − 2x + ℓ < 0 and 1 + 3ℓ − 2M ≤ 0.
- (5.14), (5.15) source reserves.

**κ-feedback** (Cycle61 eqs. 4.8–4.10). The new R depends on κ only through the old capacity pieces. So |∂R/∂κ| ≤ max 2q(1 − 2m)/(6κ²) ≤ 0.255, and |∂G/∂κ| ≤ (3/4)·d₊·0.255 < 0.16 < 1/2. The old rows add Cycle61's < 1/8. Since κ − κ_* = 2μ, every strict row gap persists at the actual supremum.

## 3. What remains to be confirmed

1. **The one new assumption (now pinned down from the Lean source).**
   - Cycle61's Lean proof states the interior affine moment `affine_fourier` (B064.lean) only for the physical numerator pair. Its plain polynomials use `rowCharacter … u` with η = `fixedSourcePrincipal`, and its analytic premise is `AffineAt F` for the source F.
   - The detector count (B023/B024, `actual_baseline_source_moments`) uses the κ-marked moment (1.6) for every character η ("No affine interior estimate … is used").
   - So the change needs Theorem 2.2 for the detector's dictionary characters η, i.e. the affine premise `AffineAt` with a general fixed η in place of the source principal. This is an extension of Cycle61's proof, not a direct citation.
   - Earlier draft of this point, kept for context: Theorem 2.2 must apply to the detector's plain witness P(U^m), the saturated witness of a dictionary character ν·ψ_u^{∓1}. Cycle61 states Theorem 2.2 for "the actual fixed parent of §1", the same canonical positive family whose moment (1.6) its detector already uses in (3.10). §3 also remarks that the physical pair is "the finite_physical_pair_positiveAt contract, not a freely selected character family". If that contract is narrower than (1.6)'s family, the change needs Theorem 2.2's proof re-run for the witness family. Its proof uses only the family's Poisson/reflection structure, the slot mesh and Lemma 2.1.
2. **The Lean formalization.** Cycle61 encodes R_κ by its closed form (countingR in B022.lean). A Lean version of this improvement would need the piecewise count and the box certificate. This has not been done, and nothing is committed or submitted.
