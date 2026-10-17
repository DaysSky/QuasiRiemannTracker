# Conditional zero-free half-planes below the frontier, certified in exact rationals

DaysSky, with Claude assistance. 17 Oct 2026. Frontier (unconditional): θ = 0.874836 (Gao, GTimG/QuasiRiemannTracker PR #62).

Every statement below has the form "hypothesis H ⟹ L(s, χ) ≠ 0 for Re s > θ_H, for every Dirichlet character χ (and the finite-order Hecke
characters of Q(ω) in the tracker's catalogue)". It is proved relative to the radical record's proof (detector-affine + timbalud's
radical count, local commit 9acccf3, θ = 0.8748389), with H replacing exactly one or more lemma inputs. "Certified" means:
- all three row inequalities G_d(δ, p) ≤ θ − v_d (d = h, d₋ = h − 10⁻³, d = 1/2; v = 10⁻⁵, 10⁻⁵, 10⁻³, with v = 10⁻⁶ on the new-capacity top row of 2″ as in the record and §73) hold on [0, 3/4] × [0, 1/2]. This is checked by exact-rational
  box subdivision, with R bounded through its monotonicity in (δ, p);
- the geometry and LOW-validity conditions hold (ℓ < x < y < 2x, y < 1/2 < h < 2y, 1/6 < ℓ < 1/5, h > 4ℓ, 4ℓ > y, M + ℓ ≤ 1, y − 2x + ℓ < 0,
  1 + 3ℓ − 2M ≤ 0);
- Cycle61's reserves (4.6) and (4.7) are positive. Cycle61 §4 (4.8)–(4.10) shows only positivity is needed; the κ-feedback is absorbed by monotonicity;
- the inverse (and, for plain-weight hypotheses, plain) slot lengths are capped by the slot supply ℓ/d₊;
- at all points except 2″ the Theorem 2.2 window (3.3) fails at d = h, so the top row uses the paper's old capacity min(d, 2y − d)/c_P (OLDTOP = 1). The junction (4.5)
  is then not needed. At 2″ the window (3.3) holds, so the top row uses the detector-affine capacity (2y − (2 − B)h)/Λ, exactly as in the record and §73. The junction (4.5),
  (4.6) > 1/10 and the affine-leg conditions (3.7) and (3.8) all pass there; (4.7) = 0.0091 > 0.

| # | hypothesis H | θ_H (certified) | h | boxes (d = h / d₋ / ½), min slack |
|---|---|---|---|---|
| 1 | LM_η: Σ*_{N(u)≍U}|M_u(U^s)|² ≪ U^{1+(5/6)(s−1)−η+ε}U^s on s ∈ [11/10, 5/4], η = 1/1000 | 0.87480 | (§73) | (§73) |
| 2 | Lemma 17.1 with r + (39/20)z ≤ m − c₁ | 0.87474 | 0.7971 | (§73) |
| 2′ | LM slope 4/5 (RD with c₀ = 1/24) | 0.87462 | 0.789 | 113 055 / 20 597 / 31, 4.7·10⁻¹⁰ |
| 2″ | Lemma 17.1 with r + (7/4)z ≤ m − c₁ (BFI strength: level 8/7 of the large-sieve barrier for the twisted Möbius mean square; marked RD with c₀ ≈ 0.048) | 0.87425 | 0.806 | 63 519 / 35 621 / 29, 3.3·10⁻⁹ (Theorem 2.2 window holds: new top-row capacity) |
| 2‴ | 2″ + 2′ together (c_I = 7/4 and LM slope 4/5: RD at c₀ ≈ 1/24 for the marked and unmarked families) | 0.87405 | 0.795 | 61 909 / 24 501 / 31, 7.0·10⁻⁹ |
| 2ᵐ | Lemma 17.1 with r + (5/3)z ≤ m − c₁ (Maynard strength: level 6/5 of the large-sieve barrier, as x^{3/5} against x^{1/2}; second condition unchanged) | 0.8740 | 0.805 | 67 837 / 33 529 / 29, 1.5·10⁻⁸ |
| 2⁗ | Lemma 17.1 with r + (5/3)z ≤ m − c₁ (Maynard strength: level 6/5 of the barrier, as x^{3/5} against x^{1/2}; marked RD with c₀ ≈ 0.067) and LM slope 4/5 | 0.8738 | 0.80 | 58 957 / 29 005 / 29, 1.5·10⁻⁸ (second condition rescaled, as in rows 5 and 7) |
| 2ˡ | Lemma 17.1 with r + (107/66)z ≤ m − c₁ (Lichtman strength: level x^{66/107}; second condition rescaled) | 0.87382 | 0.805 | 81 247 / 38 591 / 29, 1.2·10⁻⁹ |
| 2ˡ′ | Lemma 17.1 with r + (107/66)z ≤ m − c₁ (Lichtman strength: level 132/107 of the barrier, as x^{66/107} against x^{1/2}) and LM slope 4/5 | 0.87362 | 0.80 | 78 781 / 33 635 / 29, 1.5·10⁻⁸ (second condition rescaled) |
| 3 | LM slope 3/4: Σ*|M_u(U^s)|² ≪ U^{1+(3/4)(s−1)+ε}U^s on s ∈ [1, 3/2] | 0.87425 | 0.773 | 18 691 / 11 485 / 35, 9.8·10⁻⁸ |
| 4 | LM slope 7/10 | 0.87368 | 0.758 | 25 107 / 11 371 / 41, 2.1·10⁻⁸ |
| 4″ | LM slope 3/5 (RD with c₀ = 7/18) | 0.87255 | 0.7648 | 28 491 / 11 641 / 37, 2.0·10⁻⁹ |
| 4′ | LM slope 1/2 (RD with c₀ = 2/3) | 0.87148 | 0.7713 | 27 675 / 11 669 / 35, 8.1·10⁻⁹ |
| 5 | Lemma 17.1 with r + (3/2)z ≤ m − c₁ (and 2r + 5z ≤ 3m − c₂) | 0.87325 | 0.8095 | 226 197 / 86 521 / 29, 1.9·10⁻¹⁰ |
| 5′ | 5 + 3 together (c_I = 3/2 and LM slope 3/4: the RD_{1/9} representative, PLAN_RD.md) | 0.8728 | 0.7957 | 95 979 / 39 447 / 31, 5.1·10⁻⁹ |
| 6 | Lemma 18.1 with n₁ + n₂ + 2z ≤ 1 (plain slot weight 2) | 0.8734 | 0.8097 | 148 965 / 45 905 / 19, 8.6·10⁻¹⁰ |
| 7 | Lemma 17.1 with r + (6/5)z ≤ m − c₁ | 0.8728 | 0.7707 | 39 699 / 21 091 / 45, 1.1·10⁻⁸ |
| 7′ | 7 + 4 together | 0.8716 | 0.7705 | 48 759 / 21 157 / 43, 6.0·10⁻⁹ |
| 8 | 7 + 4 + 6 together | 0.8708 | 0.7753 | 63 659 / 28 417 / 39, 9.0·10⁻⁹ |
| 8′ | c_I ≤ 6/5, LM slope 1/2 and plain weight 1 together | 0.8692 | 0.7853 | 14 995 / 12 267 / 19, 1.6·10⁻⁷ (floor reserve 6·10⁻⁴: the bin floor 51/100 starts to bind) |
| 9 | DH_λ for the family on bins a ∈ (1/2 + λ/8, β*], row scales d ∈ [2/5, h + ζ] (Conditional Theorem D, DH_THEOREM.md) | 11/12 − (1−2λ)/(4(5+6λ)) → 13/15 | | margins in exact arithmetic (§61) |

**Classical-strength reading of the c_I rows (analogy; §93).** For the twisted Möbius mean square in progressions to the slot moduli pp′ (FOUR_INPUTS input 3), the large-sieve barrier is exactly c_I = 2.
A level gain by factor f gives c_I = 2/f. Caveat (DH_SEXTIC §103): the classical theorems are for fixed residue classes, while the slot structure is a character twist (one sextic character per slot prime). So the rows say what a barrier improvement by the same factor would give, not that a direct transcription of those theorems exists.
- Bombieri–Friedlander–Iwaniec (x^{4/7}, f = 8/7): c_I = 7/4 ⟹ 0.87425 (2″), or 0.87405 with LM slope 4/5 (2‴).
- Maynard (x^{3/5}, arXiv 2006.07088): c_I = 5/3 ⟹ 0.8740 (2ᵐ), or 0.8738 with slope 4/5 (2⁗).
- Lichtman (x^{66/107}, arXiv 2309.08522): c_I = 107/66 ⟹ 0.87382 (2ˡ), or 0.87362 with slope 4/5 (2ˡ′).
- Continuing the curve:
  - f = 5/4 (Lichtman's 5/8 under Selberg's eigenvalue conjecture): 0.87370 (model);
  - f = 4/3: c_I = 3/2, which is row 5, 0.87325;
  - f ≥ 5/3: c_I ≤ 6/5, where the slot supply saturates at 0.8728 (row 7).
  So the inverse-slot route alone is capped at 0.8728 by the supply ℓ/h, and further progress needs LM or DH as well.

The c_I-type hypotheses saturate at c_I = 6/5 because the slot supply ℓ/h binds. In the model (combo_cert3.py), full LM (slope 0) gives 0.86911 and slope 2/5 gives 0.87042. At slope 0 the bin floor 51/100 is within 3·10⁻⁴ of binding, so going lower needs the floor lowered (§52–§53). Hypotheses 1–8 are partial power savings beyond the proved ranges. Reaching 13/15 needs full statements: 9, the full LM (any slope at every s), or the
Möbius fourth moment at r ≈ 3/4. The single statement behind 1–5 and 7 is the rows-dominant asymptotic (RD_{c₀}) for the canonical family
(FOUR_INPUTS.md).

**Reproduce** (scratchpad/c61):
```
python3 certG.py 87425/100000 806/1000 7/4 5/6                 # 2″ (window holds, no OLDTOP)
OLDTOP=1 python3 certG.py 87405/100000 795/1000 7/4 4/5        # 2‴
OLDTOP=1 python3 certG.py 87462/100000 789/1000 2 4/5          # 2′
OLDTOP=1 python3 certG.py 874/1000 805/1000 5/3 5/6            # 2ᵐ
OLDTOP=1 python3 certG.py 8738/10000 8/10 5/3 4/5              # 2⁗
OLDTOP=1 python3 certG.py 87382/100000 805/1000 107/66 5/6      # 2ˡ
OLDTOP=1 python3 certG.py 87362/100000 8/10 107/66 4/5          # 2ˡ′
OLDTOP=1 python3 certG.py 87425/100000 773/1000 2 3/4          # 3
OLDTOP=1 python3 certG.py 87368/100000 758/1000 2 7/10         # 4
OLDTOP=1 python3 certG.py 87255/100000 7648/10000 2 3/5          # 4″
OLDTOP=1 python3 certG.py 87148/100000 7713/10000 2 1/2          # 4′
OLDTOP=1 python3 certG.py 87325/100000 8095/10000 3/2 5/6      # 5
OLDTOP=1 python3 certG.py 8728/10000 7957/10000 3/2 3/4        # 5′
CPW=2 OLDTOP=1 python3 certG2.py 8734/10000 8097/10000 2 5/6   # 6
OLDTOP=1 python3 certG.py 8728/10000 7707/10000 6/5 5/6        # 7
OLDTOP=1 python3 certG.py 8716/10000 7705/10000 6/5 7/10       # 7′
CPW=2 OLDTOP=1 python3 certG2.py 8708/10000 7753/10000 6/5 7/10  # 8
CPW=1 OLDTOP=1 python3 certG2.py 8692/10000 7853/10000 6/5 1/2   # 8′
python3 sidechecks_rad.py θ h                                  # side conditions ((3.3), (4.5) are N/A with OLDTOP; (4.6), (4.7) need > 0)
```
Logs: certG_cI175b.log, certG_bfi_a080.log, certG_maynard_a080.log, certG_maynard.log, certG_lichtman_a080.log, certG_lichtman.log, certG_a080.log, certG_a075.log, certG_a070.log, certG_a060.log, certG_a050.log, certG_cI15b.log, certG_rd19.log, certG2_w2.log, certG_cI12.log, certG_cI12_a070.log, certG2_combo2.log, certG2_deep.log (and certG2_combo.log at 0.8709).

**Reproducibility.** Re-running the row 2′ command (04:52Z, certG_a080_repro.log) gives identical box counts and minimum slacks on all three rows.

**Side conditions, re-run at every certified point (sidechecks_rad.py).** The only items flagged are:
- (3.3) window and (4.5) junction: not used with OLDTOP;
- the printed thresholds 1/10 of (4.6) and 1/100 of (4.7): only positivity is needed, and the values are 0.085–0.112 and 0.0048–0.0233.
The affine-leg conditions (3.7) and (3.8) (S_b > 0) and the case conditions hold at every certified point. A trial point c_I = 7/4 at h = 0.8108 failed (3.8) and is not
reported.

**Independent float cross-check (c61/cp/xcheck.py, separate code built on fastRs).** On a 301 × 51 grid in (δ, p), max(G_d − (θ − v_d)) at d = h is:
−1.5·10⁻⁴ (#7), −7.5·10⁻⁵ (#6), −2.8·10⁻⁴ (#3), −1.7·10⁻⁴ (#5), −1.1·10⁻⁴ (#4), −1.1·10⁻⁴ (#8), −8.4·10⁻⁵ (#7′), −1.3·10⁻⁴ (#5′; cp/xcheck_rd19.log), −8.8·10⁻⁵ (#2‴), −1.7·10⁻⁴ (#2⁗), −1.8·10⁻⁴ (#2ᵐ), −1.3·10⁻⁴ (#2ˡ′), −1.8·10⁻⁴ (#2ˡ). The other rows are smaller still, and the maxima sit at
p = 1/2, δ ≈ 0.38–0.49.

**Sharpness.** The same scripts fail just below the model optima, at the binding class p = 1/2: c_I = 6/5 at θ = 0.8726, slope 3/4 at θ = 0.8740, #5′ at θ = 0.8726 (fails at δ = 0.517, p = 1/2), #2″ at θ = 0.8742 (fails at δ = 0.4755, p = 1/2), #2‴ at θ = 0.8740 (fails at δ = 0.488, p = 1/2), #2⁗ at θ = 0.8737 (fails at δ = 0.495, p = 1/2), #2ᵐ at θ = 0.8739 (fails at δ = 0.501, p = 1/2), #2ˡ′ at θ = 0.8735 (fails at δ = 0.508, p = 1/2), and #2ˡ at θ = 0.8737 (fails at δ = 0.511, p = 1/2)
(certG_sanity_*.log). The certified values are therefore the true values of this proof under each hypothesis, to within about 10⁻⁴.

**Status of the hypotheses.** All are open, and each is beyond the corresponding known result in its closest analogue (FOUR_INPUTS.md). Numerically, every one of them
holds with no visible power loss up to U = 10⁶–3·10⁶ (DH_SEXTIC.md §31, §41, §57, §76, §82, §88).
