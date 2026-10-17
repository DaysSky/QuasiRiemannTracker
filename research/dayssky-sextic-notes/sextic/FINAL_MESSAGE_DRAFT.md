**Session 22:28Z–06:28Z: final summary**

**Frontier.** The proved boundary is still Gao's 0.874836 (GTimG/QuasiRiemannTracker PR #62), and upstream still ends at #62. I made no commits and no pushes.

**What moves it, certified in exact rationals.** Each row replaces one or more lemma inputs of the radical record's proof by an open hypothesis. For each row:
- all three row inequalities are certified by exact-rational box subdivision;
- every side condition is re-checked;
- a run just below the value fails, so each value is tight;
- an independent float re-implementation agrees.

**Classical-strength targets (new this segment).** These rows are analogies, not proofs. Each says what the sextic version of a proved classical level-of-distribution theorem would give.
The barrier for the twisted Möbius mean square in progressions to the slot moduli pp′ is exactly c_I = 2. A level gain by factor f gives c_I = 2/f. Caveat: the classical theorems concern fixed residue classes, while our slots are a character twist (one sextic character per prime). So these rows price a barrier improvement by the classical factors; they are not direct transcriptions.

| open hypothesis | θ ≤ |
|---|---|
| c_I = 7/4: Bombieri–Friedlander–Iwaniec strength (level exponent 8/7 of the large-sieve barrier, slot moduli pp′) | **0.87425** |
| c_I = 7/4 and LM slope 4/5 (rows-dominant asymptotic RD at c₀ ≈ 1/24, marked and unmarked) | **0.87405** |
| c_I = 5/3: Maynard strength (level x^{3/5} against x^{1/2}) | **0.8740** |
| c_I = 5/3 and LM slope 4/5 | **0.8738** |
| c_I = 107/66: Lichtman strength (level x^{66/107}) | **0.87382** |
| c_I = 107/66 and LM slope 4/5 | **0.87362** |

**The full ladder** (CONDITIONAL_THEOREMS.md):

| open hypothesis | θ ≤ |
|---|---|
| LM: U^{0.001} below sixth-power copies on s ∈ [1.10, 1.25] | 0.87480 |
| inverse slot cost c_I = 39/20 (Lemma 17.1 has 2) | 0.87474 |
| LM slope 4/5 (copies give 5/6) | 0.87462 |
| LM slope 3/4 / 7/10 / 3/5 / 1/2 | 0.87425 / 0.87368 / 0.87255 / 0.87148 |
| plain slot weight 2 (Lemma 18.1 has 6κ ≈ 4.5) | 0.8734 |
| c_I = 3/2 | 0.87325 |
| c_I = 3/2 and slope 3/4 | 0.8728 |
| c_I ≤ 6/5 (the slot supply saturates) | 0.8728 |
| c_I ≤ 6/5, slope 7/10 and plain weight 2 | 0.8708 |
| the same at natural strength (slope 1/2, plain weight 1) | 0.8692 |
| full LM (model; the bin floor 51/100 nearly binds) | 0.8691 |
| DH for the twist family (Conditional Theorem D) | → 13/15 |

**Results this session** (each labelled PROVED, MODEL, NUMERICAL or HEURISTIC in the files).
- **RD.** One statement sits behind LM and c_I: a rows-dominant asymptotic for the paper's canonical family (FOUR_INPUTS.md, PLAN_RD.md). Poisson maps its off-diagonal exactly onto
  the original LM sum, so it is self-dual. At the binding r = 3/4, c_I = 2 is exactly the self-dual point, and c_I = 7/4 needs the marked family only U^{(1−r)/7} past it (§93).
- **Plain slot cost 5 without any zero-free input (§92).** Average over the slot primes before positivity: the sublattices m⁴O overlap boundedly, which recovers one factor Z.
  - This strictly improves Lemma 18.1's κ = 1 case, from 6 to 5.
  - It does not beat 6κ₀ ≈ 4.5 at the record.
  - It cannot be combined with deletion or with live slots, because the gain needs the zero-slot lemma's supercritical range.
  - Even 5κ₀ would be worth only 0.874661, because the inverse leg then binds.
- **Why the classical arguments do not transfer directly (§94).** BFI-type dispersion needs a smooth variable. In the sextic setting it carries cubic Gauss sums, whose smooth sums have the Patterson main term X^{5/6}, and that term cancels only in the full Möbius combination. So the first missing step is a "Möbius orthogonality to theta" statement. Numerically (Z[ω], X up to 2·10⁵), the smooth Type I pieces carry c(D)·X^{5/6} with c(D) fitting Σ_{N(d)≤D}μ(d)Π_{p|d}(N(p)+1)^{−1}. The theta measure has the local densities of squarefree numbers (multiples of p carry 1/(N(p)+1) of it: 0.2004 against 0.2000 at N(p) = 4, and likewise up to N(p) = 25), so Möbius kills the main term through the pole of ζ_{Q(ω)} at 1. As with Σμ = o(x), that gives no rate. The concrete missing input is a power-saving error term, uniform in d, for the twisted Kubota series Σ_m g̃(dm)N(m)^{−s}.
- **The wall, located precisely (§95, §96).**
  - A power saving below X^{5/6} for Σμg̃χ needs the self-dual Type III sums that Dunn–Radziwiłł mark as open.
  - Their corrected cubic sieve reduces the critical sieve for Möbius sequences to zero-density estimates: fewer than H^{2/3} cubic characters with a zero at σ = 5/6 suffices. Heath-Brown-strength densities hit this exactly, and the loop would self-start if the zero-free region went below 5/6.
  - Numerically, balanced Type II sums with Möbius structure cancel to square-root size (|Σμg̃w| ≈ 0.2–0.6·X^{1/2} at X ≤ 10⁵), while the theta-shaped sum sits at the X^{5/6} barrier (§97). The truth is far below the barrier, but no proof exists.
  - The c_I = 2 barrier is the dual diagonal (§98). The paper's copies bound already gives the orthogonal large-sieve upper bound for the Möbius dual up to K ≤ N^{4/3}, so what is missing for c_I < 2 is orthogonality between distinct sublattice shifts (the marked RD). Input 2, by contrast, is a pure upper-bound problem: a large sieve on thin subfamilies, missing by U^{0.13} against the allowed U^{0.02}. Numerically it holds with no loss at all: the mean of |S₁S₂|² on u ≡ 0 mod m² is 0.7–1.1 times the full mean for N(m) up to U^{0.335} at the record's plain length (U = 10⁶), and 0.77–1.25 up to U^{0.355} at U = 3·10⁶ (fitted loss exponent 0.05–0.07, falling with U), where the paper's rule gives a loss of N(m). Its closest analogue, the quadratic family's fourth moment, has only a log-saving asymptotic.
  - The proved region, fed back as pointwise bounds (§99, §100), improves input 4's known loss from 3/8 to 0.325 and gives the Möbius fourth moment loss 0.40 (copies: 0.408). Both are far from the thresholds 0.205 and 0.114. The bootstrap would self-start only below θ₀ ≈ 0.755 and 0.575 respectively. The dispersion route (§32's marks outside the conductor) needs a shifted Möbius fourth moment, the same input as §88.
  - Density estimates cannot reach RD itself. Zero counts give LM slope at least 2θ − 1 + E/(s − 1), and the detector's excess E = 0.04–0.08 over DH (peaking near σ = 0.73) puts this at the copies' 5/6 fixed point (§75).
- **§88 inequality.** Σ_u|M_uQ_u|² ≤ Σ_u|M_u|² + (U·Σ_u|M_u|⁴)^{1/2} for slot polynomials up to U^{1/4}. So slots are free given the Möbius fourth moment (numerically flat).
- **LOW side (§87).** It has the steepest slope (−0.77), but improving it is equivalent to bounding the η-Möbius sum itself.
- **Input 4 corrected (§91).** The HIGH-side sixth-moment lever is worth about −2·10⁻³, not −7·10⁻³.

**Corrected along the way.**
- My quick price scripts lacked the LOW-validity, slot-supply and affine-window conditions; §89 has the corrected prices.
- The first c_I = 7/4 trial point violated (3.8). The admissible strip at 0.87425 is only h ∈ [0.806, 0.807].
- No certified value was affected.

**Literature.**
- Diaconu–Ion–Paşol–Popa: second moments only.
- de Faveri–Dunn–Hoffstein: general sextic sieves are not orthogonal.
- Conrey–Iwaniec–Soundararajan: the asymptotic large sieve for Dirichlet characters, which for Möbius coefficients needs Chowla-type input.
- Classical models for the inverse slot cost:
  - Bombieri–Friedlander–Iwaniec, level x^{4/7}, gives c_I = 7/4.
  - Maynard, level x^{3/5} ([arXiv 2006.07088](https://arxiv.org/pdf/2006.07088)), gives c_I = 5/3.
  - Lichtman, level x^{66/107} for triply well-factorable weights ([arXiv 2309.08522](https://arxiv.org/pdf/2309.08522)), gives c_I = 107/66.
- None of these supplies the sextic input.

Files:
- /tmp/claude-0/-home-user-QuasiRiemannTracker/726b9ed3-631a-533e-b4f5-2a9b6d3bc0e0/scratchpad/sextic/CONDITIONAL_THEOREMS.md (the table with reproduce commands)
- /tmp/claude-0/-home-user-QuasiRiemannTracker/726b9ed3-631a-533e-b4f5-2a9b6d3bc0e0/scratchpad/sextic/FOUR_INPUTS.md
- /tmp/claude-0/-home-user-QuasiRiemannTracker/726b9ed3-631a-533e-b4f5-2a9b6d3bc0e0/scratchpad/sextic/PLAN_RD.md
- /tmp/claude-0/-home-user-QuasiRiemannTracker/726b9ed3-631a-533e-b4f5-2a9b6d3bc0e0/scratchpad/sextic/DH_SEXTIC.md (§87–§100)
- /tmp/claude-0/-home-user-QuasiRiemannTracker/726b9ed3-631a-533e-b4f5-2a9b6d3bc0e0/scratchpad/c61/ (certG.py, certG2.py, sidechecks_rad.py and every certG_*.log)
