# Session 16 Oct 2026, 14:37–22:37Z: DH for the sextic twist family (option 1)

DaysSky, with Claude assistance. Labels as in DH_NOTE.md. No commits, no pushes. Record unchanged: 0.8748389 (9acccf3).

## Established this session (DH_SEXTIC §24–§63)

**The conditional chain, made precise.** §50–§52 correct REPORT §3 and earlier notes. The precise statement is:

> DH counts U^{1+o(1)−δ} for the sextic family, at every row scale U = Z^d with d ∈ [d_min, h+ζ] and every bin
> σ ∈ (1/2, 0.875] ⟹ θ = 13/15.

Stated with a proof outline as Conditional Theorem D′ in DH_NOTE §1. It uses two margin changes, both checked against the exact local Euler tables:
- the bin floor a0 = 51/100 moves to 1/2 + ε (valid for sr > 0.4933);
- the principal Euler region extends below 0.87 (holomorphic for Re s > 2/3).

Without these corrections the chain stalls:
- with the paper's floor, DH on [0.51, 0.875] gives 0.86979;
- the record's own counts are needed below the top row scales; with Liu's cruder count there, θ would stop at 0.8701.

With the corrections, the formula θ = 11/12 − (1−2λ)/(4(5+6λ)) is reproduced exactly (λ = 0, 0.02, 0.04).

**The DH gap in closed form (§49).** Xu–Sang's closed form for the OpenAI count equals our model to 10⁻¹⁶. It reads
λ = (α − δ)δP/(2J), with α the long Möbius slope. So DH at bin δ holds iff δ ≥ α (for c_P > 1), and a slope α′ gives DH on
σ ≥ (1+α′)/2. The whole problem is the single slope α: LM.

**Structure of LM (PROVED reductions and barriers).**
- Toolkit closure (Theorem A) and robustness under decompositions (§24, §27).
- Self-duality and fixed point (Theorems B, C).
- Weak Mertens (C′). Decompositions fail (Theorem D).
- Cost-1 slots ⟸ Möbius fourth moment (Prop. E).
- Copy efficiency 1/6 behind both α′ = 5/6 and c_P = 6κ (§37).
- Formulation check: all-elements LM would give a zero-free region 7/9 directly, so LM must be sixth-power-free (§46).
- The general-coefficient sextic large sieve is provably not sharp (de Faveri–Dunn–Hoffstein, Xu), so any improvement
  must use the arithmetic of μ (§48).

**Where LM is needed (§54–§55).**
- LM only on s ∈ [1, s0] already gives θ = 0.8732 (s0 = 1.05), 0.8714 (1.1), 0.8695 (1.16), 0.8674 (1.3).
- A power-saving asymptotic near the conductor would also improve the slope beyond s0, via copies plus subtraction of
  scales.
- So the smallest high-value target is a near-conductor asymptotic for the Möbius mean square, equivalently the dual
  Gauss-sum family just above its self-dual point.
- The paper's canonical recursion (Lemma 17.2) is an ε-lossy upper bound and does not supply it.

**The dip route (§40, §47).**
- A cluster of zeros contributes a divided difference. With phase control, Hermite–Genocchi gives a positive-measure
  representation.
- A dip needs relevant zeros spread over heights ≳ 0.2, and equally spaced progressions fail at the saturated ends.
- Closing the route needs multi-zero density and L′/L control; it is about a factor 10 short in moment order.

**Numerics (§28, §31, §38, §41–§45).**
- mean|M|²/diag stays 0.88–0.96 from rows-dominant sizes (U ≤ 10⁶) to 76× past the effective conductor (U = 10³).
- The random multiplicative control matches; the squarefree control collapses as expected.
- Slot moments show no power-type excess at r* = 0.745 up to U = 10⁹.
- This is consistent with LM and with GRH-on-average. No finite computation separates them (§44–§45).

**Literature (§48).**
- Xu–Sang 2610.15291 (0.8748867) and Liu 2610.12234 (0.874957) add no DH or LM input.
- Xu–Sang's averaged-numerator high side is already in our Cycle61-type model.

## Added after 19:00Z (DH_SEXTIC §55–§63)

- **D′ fully margin-checked (§61, PROVED relative to the paper).** Every margin of the paper's §20 was re-evaluated exactly
  at the 13/15 geometry. All remain positive. The hypothesis also shrinks: DH is needed only at row scales
  d ∈ [2/5, h+ζ], because below 2/5 the paper's own no-prime count already works (E(2/5) ≤ −0.0142). Two region edges
  move: D1 to Re s ≥ 1/2 + ε and D2 to Re s ≥ 0.86. The small rows only need ε0 = 1/3 in D1(ε0), which the paper
  already allows. Each is checked against the exact local
  Euler exponents. REPORT §3's open "re-verify small rows" item is closed (−0.0417).
- **Slot-free bridge (§63, PROVED relative to the paper).** LM_α gives DH on δ ≥ α through the paper's own no-slot
  detector at t = 3/2, which is its endpoint step (20.7). So LM_α gives θ ≤ 11/12 − (1−2α)/(4(5+6α)) → 13/15 without
  selected primes. This avoids Lemma 18.1, whose stated range is κ ≥ 3/4; Cycle61 extends it to [7/10, 3/4]. On this route,
  α < 0.0544 already beats the record.
- **Partial slopes (§55, MODEL).** Slope α′ on [1, s0] then 5/6: even α′ = 0.7 on [1, 1.05] gives 0.87467, i.e. 0.00025
  below the model baseline (three times the baseline–record gap).
- **DH windows (§59, MODEL).** DH above σ = 0.75 alone is worth nothing. A window must reach down past the binding point
  0.733 and up to about 0.8. Partial LM is worth more than the DH window it implies, so the LM slope is the right
  intermediate target.
- **LM is stronger than DH (§62, HEURISTIC).** LM needs sub-DH zero counts above its slope, or cancellation in 1/L′(ρ).
  Single rows cannot block it for s ≤ 1.334, because of the proved strip. That is unlike the ζ-analogue.
- **Record-moving frontier (§55, MODEL).** A U^{0.003} saving over copies on s ∈ [1.10, 1.25] would pass the record.
- **Twisted copies (§60, PROVED).** Cubes and squares are not copies, so sixth-power copies are the family's only exact
  symmetry and 5/6 is the copy density.
- **Function-field analogue (§58).** With big monodromy, LM holds with constant 1 as q → ∞ (Sym^k irreducible;
  Keating–Rudnick mechanism). At accessible degrees the sextic Kummer family has finite or CM monodromy, so it is no proxy.
- **Numerics (§57).** Exact over all rows at U = 3·10³: off-diagonal −1.4% to +0.6% for s ∈ [0.8, 2], up to 57× past the
  effective conductor. At the effective conductor: 0.996 ± 0.006 (U = 10⁴, all rows), 1.0013 ± 0.0034 (3·10⁴, all rows), 0.986 ± 0.006
  (10⁵, 30 000 rows); 1.004 ± 0.007 at U = 10⁵, 3.4 times past; 0.964 ± 0.025 at U = 10⁴, 121 times past.

## Open input

LM (LM_PROBLEM.md): Σ*_u |M_u(U^s)|² ≪ U^{1+ε}U^s for s ∈ (1, 3/2]. Any slope α′ < 5/6 at s ≈ 1.16 improves the model baseline; α′ ≤ 0.814 on [1, 1.16] (a U^{0.003} saving) passes the record (MODEL, §55). α′ → 0 gives
13/15. Its ζ-analogue (power-saving variance of μ in short intervals) and its quadratic-character analogue are open.

## Files

- `sextic/DH_THEOREM.md`: the conditional theorem (statement, proof relative to the paper, corollary).
- `sextic/README.md`: index.
- `sextic/DH_NOTE.md`: consolidated statements.
- `sextic/DH_SEXTIC.md`: full log.
- `sextic/DH_SUMMARY.md`: status.
- `sextic/LM_PROBLEM.md`: one-page problem statement.
- Scripts: `sextic/floor/`, `sextic/dip/`, `sextic/mmnum/lm_beyond.py`, `lm_exact.py`, `mmk3.c`, `sextic/ff/`, `sextic/margins/`,
  `c61/lmrange.py`, `c61/dhrange.py`, `c61/lmfront.py`, `c61/lmeta.py`, `c61/lmwin.py`, `sextic/floor/scales25.py`.
