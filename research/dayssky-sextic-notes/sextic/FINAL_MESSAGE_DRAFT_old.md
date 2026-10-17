**Session 22:28Z–06:28Z: final summary**

**Frontier.** The proved boundary is Gao's 0.874836 (GTimG/QuasiRiemannTracker PR #62). Upstream still ends at #62. Everything below is either a proved
reduction or a conditional bound certified in exact rationals. No commits or pushes.

**Conditional bounds certified in exact rational arithmetic** (relative to the radical record's proof; all three row inequalities plus side conditions):

| hypothesis (open) | θ ≤ |
|---|---|
| LM: U^{0.001} below sixth-power copies on s ∈ [1.10, 1.25] | 0.87480 |
| inverse slot cost c_I = 39/20 (Lemma 17.1 has 2) | 0.87474 |
| LM slope 4/5 (copies give 5/6) | 0.87462 |
| LM slope 3/4 | 0.87425 |
| LM slope 7/10 | 0.87368 |
| c_I = 3/2 | 0.87325 |
| plain slot weight 2 (Lemma 18.1 has 6κ ≈ 4.5) | 0.8734 |
| c_I ≤ 6/5 (the slot supply saturates) | 0.8728 |
| c_I ≤ 6/5 and LM slope 7/10 | 0.8716 |
| all three critical-length inputs together (c_I ≤ 6/5, slope 7/10, plain weight 2) | 0.8708 |
| the same at natural strength (c_I ≤ 6/5, slope 1/2, plain weight 1) | 0.8692 |
| DH for the family on d ∈ [2/5, h+ζ] (Conditional Theorem D) | 13/15 |

**Proved this session.**
- **One statement behind inputs 1 and 3 (FOUR_INPUTS.md).** It is a rows-dominant asymptotic for the paper's canonical family (RD_{c₀}): c₀ = 1/9 gives LM slope 3/4. With marks it gives the inverse
  cost.
- **§88, a new inequality:** Σ|M_uQ_u|² ≤ Σ|M_u|² + (U·Σ|M_u|⁴)^{1/2} for slots up to U^{1/4}. Given the Möbius fourth moment (numerically flat), slots are free. For a single plain
  block it is unconditional: slots up to U^{1/4} are free, against 0.11 in Lemma 18.1.
- **§87, the LOW side.** It has the steepest slope (−0.77), but beating its Cauchy–Schwarz step is equivalent to bounding the η-Möbius sum itself, so it is not a separate lever.
- **§89, joint prices** with every condition of the proof enforced. Three critical-length improvements at once reach about 0.8708. Only supercritical inputs (full LM, DH, M4) reach 13/15. The long-block sixth
  moments are worth about −2·10⁻³ once the slot term is counted (§91, which corrects §81's −7·10⁻³).

**Checks.** Each certificate fails just below its model optimum (c_I = 6/5 at 0.8726, slope 3/4 at 0.8740), so the bounds are tight. An independent float re-implementation
agrees at all seven points (worst margin −7.5·10⁻⁵).

**Corrected along the way.** My quick price scripts lacked the LOW-validity, slot-supply and affine-window conditions. The corrected numbers are in §89. The certificates were
never affected.

**Literature (Oct 2026).** Diaconu–Ion–Paşol–Popa cover only second moments, and not the sextic family over Q(ω). de Faveri–Dunn–Hoffstein show the general
sextic sieve cannot be orthogonal. The sixth-moment input is beyond ζ's best (loss 1/4). The inverse-cost input has the shape of Bombieri–Friedlander–Iwaniec
for well-factorable moduli pp′.

Files: …/scratchpad/sextic/FOUR_INPUTS.md, DH_SEXTIC.md §87–§90, SESSION_2228.md, …/scratchpad/c61/certG*.py and the logs.
