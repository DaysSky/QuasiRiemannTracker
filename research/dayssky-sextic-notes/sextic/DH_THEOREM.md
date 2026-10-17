# A density hypothesis for the sextic twist family gives Re s > 13/15

DaysSky, with Claude assistance. 16 October 2026. This is a conditional result. It is proved relative to the OpenAI paper
"The Quasi-Riemann Hypothesis" (30 Sep 2026), whose argument is used unchanged apart from the points listed in §3. Full
derivations are in `DH_SEXTIC.md` §49–§61. The checks are `margins/part2_margins.py` (exact rational arithmetic) and
`floor/scales.py`.

## 1. Hypothesis

Notation as in the paper, Part II (§12, §19, §20):
- O = Z[ω];
- rows u are sixth-power-free with q_u ≍ U = Z^d;
- a row lies in the dynamic bin (i, a) when its selected zero (bounded height) has real part in the bin around a;
- δ = 2a − 1;
- β* is the supremum of real parts of zeros of primitive finite-order Hecke characters over Q(√−3), and β* ≤ 11/12 by
  Part I.

**DH_λ** (λ > 0). Put a0 = 1/2 + λ/8. For every admissible bin width e > 0, every bin a ∈ (a0, β*], every row scale
d ∈ [2/5, h + ζ] and every prime-amplitude class (main-slot mean q ∈ [0, δ/2]), the number of retained rows satisfies

  #C ≪_ε U^{1 + λ − δ + ε} (1 + T1)^{A}.

The trivial count U^{1+ε} already gives this for δ ≤ λ. So DH_λ has content only for σ = (1+δ)/2 > (1+λ)/2.

## 2. Theorem

**Theorem (conditional).** Assume DH_λ. Then every finite-order Hecke L-function over Q(√−3), in particular ζ and every
Dirichlet L-function, has no zeros with

  Re s > 11/12 − ℓ/4  for every ℓ < (1 − 2λ)/(5 + 6λ).

So θ(λ) = 11/12 − (1 − 2λ)/(4(5 + 6λ)): 0.866826 at λ = 0.001, 0.869792 at λ = 0.02 and 0.872774 at λ = 0.04. If DH_λ
holds for every λ > 0, the zero-free half-plane is Re s > 13/15.

## 3. Proof (relative to the paper)

Take the paper's Part II with the geometry σ0 = B = 11/12 − ℓ/4, x = l_x in the open interval

  (13ℓ − 1 + 12(1 + ℓ)λ)/(4 + 12λ) < x < (1 − ℓ)/2,

with h = 1 − x + ℓ and l_y = 1 − ℓ − x. Argue by contradiction with β* > σ0, as in paper (12.1). Every estimate of the
paper is kept, with these changes.

1. **Moderate rows, top scales (Lemma 20.1 with the DH_λ count).** Paper (20.4) at d = h with R ≤ 1 + λ − δ and q ≤ δ/2
   gives E(h) ≤ C0 + hλ + δ(x + ℓ/2 − 1/2), where C0 = −x/3 − 1/12 + 13ℓ/12. The left inequality on x is C0 + hλ < 0, and
   the right one makes the δ-coefficient negative. E is linear in d, and E(2/5) < 0 with the same count, so all scales
   d ∈ [2/5, h + ζ] follow, with ζ chosen after the margins as in paper Prop. 20.3.
2. **Moderate rows, d ∈ [1/100, 2/5].** Use the paper's own no-prime count R = 76/75 − 2δ/3 (§20.5, t = 1). At d = 2/5,
   E ≤ −0.0142 for δ ≤ 2σ0 − 1, and E < (1+δ)/2 − σ0 for 2σ0 − 1 < δ ≤ 5/6. This is the requirement of paper (20.12).
3. **Floor.** The floor bin a0 = 1/2 + λ/8 with R = 1 gives E ≤ −0.0006 (λ = 0.001). Its use requires the D1 edge
   Re s ≥ 51/100 of paper Lemma 7.1 to be lowered to 1/2 + ε.
   - At the central contours the local exponents are (−1−10e, −2.54, −2.54, −1.04−96e, −2.04) at a = 1/2. All are < −1,
     and the binding one holds for a > 0.4933. Row primes give O(1) per prime.
   - Every occurrence of 51/100 in the paper was checked to be a margin: Lemma 7.1, Props. 8.3, 9.2 and 16.1, Lemma 19.1,
     and (20.5) (DH_SEXTIC §52–§53).
4. **Principal row.** The D2 edge Re s ≥ 7/8, and Lemma 10.5's σ0 ≥ 7/8, are lowered to 0.86.
   - At (w, z) = (1, 1/6) the local exponents are (−1−s, −1−s, −2, 3−6s, −1−s), all < −1 for Re s > 2/3.
   - At the D2 corner the largest is −1.807 at s = 13/15.
   - B_p = −1 + O(q_p^{−0.867}).
   - The margins m_w = l_y/20 and m_z = h/600 of (20.3) stay positive.
5. **Small rows.** The line (β*+e, 1/2, 17/50) lies in D1(1/3) rather than D1(3/8). The paper's Definition 10.1 already
   gives holomorphy and bounds on every D1(ε0), ε0 > 0, so this is a choice of ε0, not a region change.
   - The exponents in the proof of (16.15) are −0.867, −2.04, −2.37 and −1.54 for selected primes, all ≤ 0. For unselected
     primes they are all < −1.
   - The margin h(17/50 − 1/6) − l_y/2 + 2d_min is −0.0417.
6. **Low side (Prop. 15.3).** The bound Z^{l_x/2 + b/12 + ε}, with b = l_y − l_x > 0, sits exactly β* − B below C(β*), since
   C(s) − l_x/2 − b/12 = s − 11/12 + ℓ/4. The proof's length inequalities hold with room ≥ 0.11, including the
   completed-row excess (5ℓ − 1 + d)/8 − d ≤ 0 on [0, ℓ].
7. **Unchanged.** Large rows, the prime normaliser A_T, κ_P = 7ℓ/(16K), the prime supply ℓ/(h+ζ) > 7/37, and the order of
   choices (Prop. 20.3, Lemma 11.1).
8. **Not used.** The selected-prime row moments (Lemma 18.1 mesh, Prop. 19.2 with selected primes, the endpoint (20.7)) are
   replaced by DH_λ. This matters because κ = 2β* − 1 can now lie in (11/15, 3/4), below the range [3/4, 1] where Lemma 18.1
   is stated. The Cycle61 paper works on κ ∈ [7/10, 3/4] after importing β* ≤ 7/8 (its (1.3), (4.8)), so the
   detected-count route also stays available. DH_λ simply makes it unnecessary. The slot system stays only in the compensation (high correction and normaliser), which does not depend on κ.
   The no-prime count of item 2 (t = 1: max_r min(1 − δr, 1 − 2δ(1 − r)) = 1 − 2δ/3) uses only the in-bin caps and Hecke
   reflection, and is independent of κ.

All margins in items 1–7 are evaluated in exact rational arithmetic for λ = 0.001, 0.01, 0.02 and 0.04 (`part2_margins.py`,
no failures), at ℓ = ℓ_sup − λ/100 with x the midpoint of its interval. An independent floating-point optimiser
(`floor/scales25.py`) with the same split reproduces θ = 0.866986 at floor 0.501 and 0.866699 at 0.5001. Only items 1 and 3 depend on λ at leading order. For
smaller ℓ (x fixed) both C0 and h decrease, so they only improve. The other margins stay at least 0.0013. The same code reproduces the paper's own numbers at its geometry: (20.3), (20.5), (20.6) and (20.11). The
contradiction then closes exactly as in paper §20.6, which gives β* ≤ σ0. ∎

**Level of the region checks.** Items 3 and 4 change numerical region boundaries; item 5 only changes ε0. The claim is that the paper's proofs of
Lemma 7.1, (16.13) and (16.15) go through with the changed boundary, because every exponent they use stays on the
required side. This was checked from the exact local Euler formulas, not by rewriting those proofs line by line.

## 3a. Corollary: LM gives DH_λ without selected primes

**Corollary.** Suppose that at every row scale U = Z^d, d ∈ [2/5, h+ζ], LM_α holds:
Σ*_u |M_u(U^s)|² ≪ U^{1+α(s−1)+ε}·U^s for s ∈ [1, 3/2]. Then DH_α holds, so θ ≤ 11/12 − (1−2α)/(4(5+6α)), which tends to
13/15 as α → 0.

*Proof.*
- Use the paper's detector (Prop. 8.3) with t = 3/2, exactly as in its endpoint step (20.7). The Möbius witness then has
  length U^r with r ∈ [1 − O(ε), 3/2].
- The no-slot amplified estimate (19.3) with slope α counts the bin-δ rows by U^{1 − α + (α − δ)r + ε}. For δ ≥ α this is
  at most U^{1−δ+ε}, the maximum being at r = 1.
- For δ < α, the trivial count U^{1+ε} is at most U^{1 − δ + α + ε}.
- Hence DH_α holds, and the Theorem applies with λ = α.

No selected-prime moment (Lemma 18.1) enters, so the restriction κ ∈ [3/4, 1] of that lemma is irrelevant. This route is
crude for fixed α > 0, since the uniform λ = α is much worse than the detected count of DH_SEXTIC §49. It gives 0.86825 at
α = 0.01, 0.87130 at α = 0.03 and 0.87421 at α = 0.05, and it passes our local record 0.8748389 only for α < 0.05441 (and the current upstream frontier 0.874836, Gao, only for
α < 0.05439). The detected
count (MODEL, §55) needs only a U^{0.003} saving over copies. ∎

## 4. Status of the hypothesis

DH_λ is open. It is a density hypothesis for a thin family (family size ≍ conductor); no such result is known for any
fixed-order family, quadratic included. In this architecture it follows from the long Möbius mean square (LM):

  Σ*_u |M_u(U^s)|² ≪ U^{1+ε} U^s  for s ∈ (1, 3/2],

via the closed form R = 1 − δ + (α − δ)δP/(2J) of the detected count (Xu–Sang). See `LM_PROBLEM.md`. The current slope is
α = 5/6, from sixth-power copies. Within the theta-reflection and large-sieve toolkit this is a fixed point. Numerically,
LM holds to within the row-fluctuation level up to 57 times past the effective conductor (DH_SEXTIC §57).
