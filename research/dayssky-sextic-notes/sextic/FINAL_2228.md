**Superseded (05:40Z): the current end-of-session summary is FINAL_MESSAGE_DRAFT.md. This file is kept as the 03:40Z draft.**

DRAFT of the end-of-session message (22:28Z–06:28Z). Updated through 03:40Z; final numbers filled in at the end.

**Frontier.** The proved boundary is Gao's **0.874836** (GTimG/QuasiRiemannTracker PR #62, merged 16 Oct). It contains our detector-affine and
timbalud's radical count. The latest upstream PR is still #62. The local 0.8748389 (9acccf3) is superseded and was never pushed.

**What moves it, with values (MODEL on the record pipeline; all geometry, LOW-validity, floor, slot-supply and affine-window conditions enforced):**

| statement (open) | θ |
|---|---|
| none | 0.874836 |
| inverse slot cost c_I = 1.8 / 1.5 / 1.2 (saturates at the slot supply ℓ/h below ≈ 1.2) | 0.874353 / 0.873204 / 0.872762 |
| LM slope 0.80 / 0.75 / 0.70 instead of the copies' 5/6 | 0.874589 / 0.874186 / 0.87366 |
| plain slot weight 3 / 2 / 1 instead of 6κ ≈ 4.5 | 0.874313 / 0.873382 / 0.870524 |
| sixth moments of the AFE blocks m ∈ (1/6, 1/2] with small loss (§91) | roughly −2·10⁻³ |
| full LM, DH for the family, or the Möbius fourth moment at r ≈ 0.75 | 13/15 = 0.866667 |

Certified in exact rational arithmetic, relative to the radical record's proof (§73, plus this session):
- U^{0.001} over copies on [11/10, 5/4] ⟹ θ ≤ 0.87480;
- c_I = 39/20 ⟹ θ ≤ 0.87474;
- **LM slope 3/4 (instead of the copies' 5/6) ⟹ θ ≤ 0.87425** (§90; all three rows certified, side conditions pass, the top row on the paper's old capacity);
- **LM slope 7/10 ⟹ θ ≤ 0.87368**;
- **inverse slot cost c_I ≤ 6/5 ⟹ θ ≤ 0.8728** (the slot supply saturates there);
- **inverse slot cost c_I = 3/2 ⟹ θ ≤ 0.87325**;
- **all three critical-length inputs together (c_I ≤ 6/5, LM slope 7/10, plain weight 2) ⟹ θ ≤ 0.8708**;
- **plain slot weight 2 ⟹ θ ≤ 0.8734**.

**Proved this session (beyond the earlier half).**
- §87, the LOW side. It has the steepest slope of any lever (dθ/dΔ ≈ −0.77), and its whole b/12 term comes from the sixth-power exceptional
  frequencies of Prop. 15.2. Beating Cauchy–Schwarz there is the same as bounding the η-Möbius sum itself (REPORT §12.4), so it is no independent lever.
  The HIGH constants 13/75 and 17/50 cancel at the binding row d = h.
- §88, a large sieve over the slot twists. Σ_u|M_uQ_u|² ≤ Σ_u|M_u|² + (U·Σ_u|M_u|⁴)^{1/2} for slots up to U^{1/4}. Given the Möbius fourth moment,
  slots therefore cost nothing. Unconditionally the square root loses (U^{0.708} against Lemma 17.1's U^{0.592} at the binding point). For a single plain block, the
  corollary is unconditional: slots up to U^{1/4} are free, against Lemma 18.1's 0.11.
- §89, joint prices.

**Literature checked (CITED; abstracts read as data).**
- Diaconu–Ion–Paşol–Popa (2607.27131): twisted second moments for r-th order families over fields containing μ_{2r}. That covers cubic over Q(ω), not
  sextic, and second moments only.
- de Faveri–Dunn–Hoffstein (2607.07911): the general cubic and quartic large sieves are not orthogonal, and they conjecture the same for all orders. So the inputs must use the
  arithmetic of μ.
- Sang (2610.15291, θ = 0.8748867) and McGoldrick (2610.17144) supply none of the inputs.
- Analogue status: the sixth moment needed (λ < 0.20) is beyond ζ's best (loss 1/4) and is open in every family whose size equals its conductor (ζ, a
  single modulus, quadratic characters). It is known only in thick families (Huxley; Stucky for Γ₁(q)).

**Correction.** This session's quick price scripts (c61/cp/price*.py) lacked the geometry and LOW-validity checks, the slot-supply cap and the affine-window
check. All values above come from the corrected, window-aware runs (c61/cp/combo_cert2.py, combo_wa.py). Of the earlier numbers only c_I = 1 changed
materially (0.870184 → 0.872762, because of the supply cap). The certificates and the §71 table were never affected.

Files: sextic/DH_SEXTIC.md §83–§89, FOUR_INPUTS.md, PLAN_INPUT2.md, POST_DRAFT.md, SESSION_2228.md, README.md, c61/cp/*.py and *.log, c61/certG.py. No
commits or pushes.
