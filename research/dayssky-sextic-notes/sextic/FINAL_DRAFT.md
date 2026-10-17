The 8-hour session ran 14:37–22:37Z. The main result of this block: θ = 13/15 is now a checked conditional theorem relative to
the OpenAI paper, with one precisely stated hypothesis, and that hypothesis is smaller than before.

**Proved, relative to the paper** (`DH_THEOREM.md`)
- **Hypothesis DH_λ:** the sextic twist family has at most U^{1+λ−δ+ε} rows in each bin a > 1/2 + λ/8, at row scales
  d ∈ [2/5, h+ζ].
- **Conclusion:** ζ, every Dirichlet L-function and every finite-order Hecke L-function over Q(√−3) are zero-free on
  Re s > 11/12 − (1−2λ)/(4(5+6λ)). As λ → 0 this gives 13/15.
- **Checking:** I re-evaluated every margin of the paper's §20 in exact rational arithmetic at the new geometry, and all
  are positive. The same code reproduces the paper's own numbers (20.3), (20.5), (20.6) and (20.11) exactly. An
  independent optimiser gives 0.866986 → 0.866699 → 13/15 as the floor goes to 1/2.
- **Changes to the paper:**
  - the D1 edge moves from 51/100 to 1/2;
  - the D2 edge moves from 7/8 to 0.86;
  - both moves are checked against the exact local Euler exponents;
  - small rows take ε0 = 1/3, which the paper already allows;
  - below d = 2/5 the paper's own count suffices.
- **Corollary (§3a):** LM at slope α at those scales gives DH_α directly, through the paper's own no-slot endpoint
  detector (t = 3/2). No selected-prime machinery is used, and θ → 13/15 as α → 0.

**The remaining input:** DH_λ, which follows from LM (the long Möbius mean square past the conductor). New this block:
- **Smallest record-moving form (MODEL):** a saving of U^{0.003} over the sixth-power-copies bound at lengths U^s,
  s ∈ [1.10, 1.25].
- **DH windows:** DH above σ = 0.75 alone is worth nothing. Partial LM is worth more than the DH window it implies.
- **Twisted copies (PROVED):** they give nothing, so 5/6 is exactly the copy density.
- **LM is stronger than DH (HEURISTIC):** LM needs sub-DH zero counts above α, or cancellation in 1/L′(ρ). Single rows
  cannot block it for s ≤ 1.334, because of the proved strip. That sets it apart from the ζ-analogue.
- **Function fields:** with big monodromy, LM holds with constant 1 as q → ∞ (Katz equidistribution, the
  Keating–Rudnick mechanism). At accessible degrees the sextic Kummer families are CM, so they are not a proxy.

**Numerics (exact coprime diagonals)**
- U = 3·10³, all rows: the off-diagonal is between −1.4% and +0.6% up to 57× past the effective conductor.
- At the effective conductor: 0.996 ± 0.006 (U = 10⁴, all rows), 1.0013 ± 0.0034 (U = 3·10⁴, all 108 606 rows) and
  0.986 ± 0.006 (U = 10⁵).
- At U = 10⁵, 3.4× past the conductor: 1.004 ± 0.007.
- At U = 10⁴, 121× past the conductor: 0.964 ± 0.025.

**Record:** 0.8748389 (commit 9acccf3, local, not pushed). Xu–Sang (0.8748867) and Liu (0.874957) are both above it. No
commits or pushes this session.
