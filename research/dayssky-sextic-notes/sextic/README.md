# Sextic route to θ = 13/15: index

DaysSky, with Claude assistance. Read in this order:

1. `DH_THEOREM.md`: the conditional theorem. DH_λ for the sextic twist family, at row scales d ∈ [2/5, h+ζ], gives the
   zero-free half-plane Re s > 11/12 − (1−2λ)/(4(5+6λ)), hence 13/15. It lists every change to the OpenAI paper, and every
   margin is checked in exact arithmetic (`margins/part2_margins.py`). Corollary §3a: LM at slope α gives DH_α through the
   paper's own no-slot detector, with no selected primes.
2. `LM_PROBLEM.md`: the open input, one page. It covers the long Möbius mean square beyond the conductor, its smallest
   record-moving form (a U^{0.003} saving over copies for s ∈ [1.10, 1.25]), proved barriers, analogues and evidence.
3. `DH_NOTE.md`: consolidated statements and proofs. Theorems A–D cover closure, self-duality, the fixed point and
   decompositions; Proposition E (cost-1 slots from the Möbius fourth moment), numerics and θ tables are also there.
4. `DH_SUMMARY.md`: short status.
5. `SESSION_1437.md`: what the 16 October 14:37–22:37Z session added. `SESSION_2228.md`: the 22:28Z–06:28Z session (exact
   prices of every open input at the record point, conditional certificates, the LM dual as a sharp sieve, Möbius–Gauss numerics).
6. `DH_SEXTIC.md`: the full working log, §1–§103, with labels PROVED / MODEL / NUMERICAL / CITED / CONJECTURE.
7. `FOUR_INPUTS.md`: the four open inputs that would move the frontier substantially, with values, mechanisms, evidence, literature status
   and their closest analogues.
8. `CONDITIONAL_THEOREMS.md`: conditional half-planes below the frontier, certified in exact rationals (classical strength: c_I = 7/4 (BFI-type) ⟹ 0.87425, with LM slope 4/5 ⟹ 0.87405; c_I = 5/3 (Maynard-type) ⟹ 0.8740, with slope 4/5 ⟹ 0.8738; c_I = 107/66 (Lichtman-type) ⟹ 0.87382, with slope 4/5 ⟹ 0.87362; LM slope 4/5 ⟹ 0.87462; slope 3/4 ⟹ 0.87425;
   c_I = 3/2 ⟹ 0.87325; c_I = 3/2 with slope 3/4 ⟹ 0.8728; c_I ≤ 6/5 ⟹ 0.8728; plain weight 2 ⟹ 0.8734; three inputs together ⟹ 0.8708, at natural strength ⟹ 0.8692; DH ⟹ 13/15),
   with the commands to reproduce them, sanity failures just below each value, and an independent float cross-check.
9. `PLAN_RD.md`: research plan for the single statement (RD) behind inputs 1 and 3, ending with the exact inequality whose proof gives θ ≤ 0.87425 (c_I = 7/4), checked numerically at the binding point. Its certified values, why the paper's tools stop, and a BFI-type dispersion with metaplectic Kuznetsov route.
9a. `OPEN_PROBLEM_DECORRELATION.md`: the single estimate (★) behind both remaining slot-cost inputs, with equivalent forms, numerics, obstructions and literature (17 Oct).
10. `PLAN_INPUT2.md`: research plan for the plain slot cost (twisted fourth moments on the cubic/sextic cover). `POST_DRAFT.md`: outreach draft
   of the four statements (not posted).

Scripts:
- `margins/`: exact margin check against paper §20.
- `floor/`: floor and row-scale models.
- `mmnum/`: Möbius mean-square numerics. The kernels are `mmk*.c`; `lm_exact.py` gives exact coprime diagonals.
- `ff/`: the function-field analogue.
- `dip/`: the zero-dip model.
- `../mugauss/`: Möbius-weighted cubic Gauss sums over Z[ω] (Patterson bias test). `../mugauss_t1/typeI.py`: theta bias of Type I pieces (§94 addenda).
- `../c61/`: θ models (`capc4.py`, `lmrange.py`, `lmfront.py`, `lmeta.py`, `lmwin.py`, `dhrange.py`); exact record pipeline
  variants `thrad.py`, `thLM.py` (partial LM), `thcI.py` (slot cost), `thq0.py`/`thBL.py` (Theorem 2.2 threshold); conditional
  certificates `certLM.py`, `certcI.py`, `sidechecks_rad.py`; `cp/combo_wa.py` prices the four levers jointly on the record pipeline
  with the LOW-side validity constraints enforced (`cp/combo_wa.log`); `cp/combo_cert2.py` adds the slot supply, the affine window and the floor
  (`cp/wa_*.log`); `certG.py` is the exact-rational certificate with inverse cost c_I and LM slope α as parameters (`OLDTOP=1` puts the
  paper's old capacity on the top row when the Theorem 2.2 window fails there), logs `certG_*.log`.

Frontier: 0.874836 (Gao, GTimG/QuasiRiemannTracker PR #62, merged 16 Oct). It supersedes the local 0.8748389 (commit 9acccf3, not pushed).
