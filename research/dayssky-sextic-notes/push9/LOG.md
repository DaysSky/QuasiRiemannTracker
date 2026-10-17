# 9-hour push log (started 13:15Z)

## Side lever (sextic), recorded before switching to (L³) at the user's request (~14:05Z)
- Upstream merged PR #60 (Gao, same PR50+PR56 combination, 0.874840). Our local 0.8748389 is still 1.1e-6 lower; the PR body needs a
  predecessor update before submission.
- B-bootstrap idea: feed the affine moment back into the high-zero count of Theorem 2.2 to raise B. Float model (combo/gen/fmodel2.py):
  dB = 2e-4, 5e-4, 1e-3 give θ = 0.8748373, 0.8748350, 0.8748345 (from 0.8748389); dB ≥ 2e-3 becomes infeasible in the model.
  The q0 trade-off (T2 q0 = 0.66) gives 0.8748364. So this lever is worth at most ~5e-6. Parked.

## (L³) — priority from 14:05Z
- (L³) round 3 logged as S42–S48 in quartic/L3_ATTEMPT.md; report §10e updated. No proof; obstruction sharpened (S46, S48).
- Submission maintenance (~14:50Z): local commit rebased onto upstream cb1bd454 (PR #60 merged) -> 9acccf34, author/committer DaysSky,
  backup branch backup-pre-rebase-0267c799. prepare.mjs at 9acccf34: 3819 modules, 3826 source hashes identical to the verified build,
  Solution.lean identical up to the head sha. PR body updated (credits PR #60, upstream cb1bd454, previous best 0.874840); admission checks
  pass (decimal vs gao-affine-asymmetric-20261015, metadata, proof-only). Not pushed.
- 22:12Z S49 logged as CANDIDATE (Poisson+BGL unbalanced quartic LS; (L3) Q^{13/8} conditional on verification items a-c). Window end.
- 22:14Z S49: Poisson step = BGL Lemma 3.3; improvement = Lemma 3.4 with l=3 + exact j-box multiplicity. Still CANDIDATE (items i-iii).
- 22:30Z S49 PROVED (Thm 1 unbalanced quartic LS; Thm 2 (L3) sum << Q^{13/8}); proof l3c/S49_PROOF.md; no theta consequence.
- 22:40Z θ levers after S49.
  - The θ chain uses only the cubic (Heath-Brown, sharp in all ranges by dFDH Thm 1.1) and quadratic sieves, so S49's trick has nothing to improve there.
  - Theorem 2.2 is a one-parameter family in q0. With Ka = q0/(9/2−6q0), Λ = 1+Ka X0·(7/2)/(q0+Ka X0) and B0 = (1−ρ0+Ka X0)/(q0+Ka X0), where X0 = 5−6ρ(q0),
    this reproduces (2.7) exactly at q0 = 17/25. Lemma 2.1 allows q0 ≥ 0.65 (real part ≥ 33/40).
  - Plain family (no bootstrap), B = B0(q0)+7/10000: q0 = .65 gives 0.87483979, worse than the 0.87483886 base. Scan continuing (combo/gen/q0plain.out).
  - Bootstrap T2 (feed (Λ,B) into the low bins): q0 = .655 / .66 / .665 / .67 give 0.87483757 / 0.87483641 / 0.87483504 / 0.87483745.
    The best is ≈3.8e-6 below the submission. It needs new math: low bins after greedy deletion can leave Theorem 2.2's window z ≥ .034,
    so the argument needs two regimes. q0 is hard-coded in 33 Lean modules, including timbalud's marked-count files and the 7/10000 B add-on.
- 22:45Z q0 lever details.
  - timbalud's saving rate is 519/4550 = q0 − (3/2)ρ(q0), tight at 17/25 (CoupledRadicalSaving.first_common_uniform). The 37/10000 is 0.033 × that rate,
    so lowering q0 shrinks the B add-on. The model must use dB = 0.0037·κ(q0)/κ(.68)/(q0+Ka X0).
  - Plain family with fixed dB = 7e-4: q0 = .65/.66/.67/.68/.69/.70 give 0.87483979/0.87483804/0.87483772/0.87483886/0.87484140/0.87484529.
  - Bootstrap window: after greedy deletion z can fall to ≈ .021, below Theorem 2.2's .034. For z < (1−.88)/4.5 ≈ .0267 the old moment applies directly.
    Cells need g_j > 0, i.e. z > (1−B)/(9/2−Λ) ≈ .0225. So the bootstrap needs Theorem 2.2 re-proved with cells from ≈ .0267 to .04, plus a new low-bin ledger.
- 22:55Z Plain q0 lever, with the B add-on scaled by timbalud's rate (combo/gen/q0plain3.out): q0 = .675 → 0.87483840, .67 → 0.87483841.
  The base is 0.87483886, so the gain is ≤ 4.6e-7. Lever closed: not worth re-proving the ~33 q0-dependent modules.
  Remaining small levers, ranked:
  - bootstrap T2 (≈3.8e-6): needs Theorem 2.2 window extension plus a new low-bin ledger;
  - detector marked count (dθ/dR ≈ 0.135; size unknown): needs new math in the detector fiber count.
  Both are large Lean projects.
- 23:12Z S50 logged: (L3) beyond 13/8 <=> T-theta correlation (P4); conditional Q^{1.543}; S44 orientation flagged.
- 23:41Z S51 (G2 WMDS, Chinta–Gunnells match), classification r=1,2,3 -> finite/affine/hyperbolic, theta extraction numerics inconclusive; report status block added.
- 23:46Z gen_cauchy: ||A*R||^2 = diag*(1-0.04..0.07); Cauchy loss = 2.2*Z/Q at Z=1e5..1e6 -> mechanism confirmed; conditional route = optimal LS shape for M<=0.9.
- 23:49Z twist growth: |X(h0)|/sqrt(S) flat (rms 0.74-1.09) for N(h0)<=100=H -> B'~0.
- 23:50Z L3 hour summary: S50 (13/8 = T-theta correlation), S51 (obstruction series = Chinta-Gunnells G2 WMDS, n=4, symbol eps^3; polar set = theta divisor; residue centre 3/8; (L3) << Q^{3/2} conditional on CG extended to l=1 twist & Q(i)); classification r=1,2,3 finite/affine/hyperbolic; numerics: Cauchy mechanism confirmed (||A*R||^2 = diag*(1-5%)), B'~0, first moment real & square-root size; theta extraction inconclusive.
- 23:57Z S51 §17: DDHL (4.28) checked against the twist rule (d: (lc,ld)->(lc+ld,-ld); c: (-lc, ld+3lc)). Gamma normaliser N = prod over positive root functionals of G4(1/2+f, |f(l0)|); |twists| y:0, x:1, x+y:1, 2x+y:2, 3x+y:3, 3x+2y:3 (l3d/g2_gamma.py, well-defined, one flip per simple reflection). 2-adic class group Q2(i)^x/4th powers has order 256 (matches the §16 conductor scale). A_ji poles on Re=centre/theta line: clear with 2-adic factors; exponents unchanged.
- 23:58Z S51 §18: mu8 ⊄ Q(i) costs only the quadratic sign (c,-1)_S = chi_c(-1) = (-1)^{(N c-1)/4} (class character mod 4; Omega-invariant, so it preserves M_t(Omega)). Open: one finite check that CG's sigma_c, sigma_d with the sign kept still satisfy (sigma_c sigma_d)^6 = 1 on the order-256 2-adic class space.
- 00:03Z S51 §19 (quartic/l3e/cg_local.py): exact check of CG Def 3.1 for G2, n=4. With s=chi_pi(-1)=-1, raw Gauss sums FAIL (28/48, sigma_i^2 != 1); gamma*s^i FAILS; gamma*s^{binom(i,2)} (Brylinski-Deligne (pi,pi)^{C(k,2)}) PASSES all for 4 twisting params x 2 seeds; s=+1 reproduces CG Thm 3.2. Global effect: chi_c(-1) on the c-node, absorbed into Psi. Task #26 remaining: global FE write-up, p-part comparison at non-coprime (c,d), Bochner+PL.
- 00:08Z S51 §20: CG p-part for G2 n=4 computed exactly (14 monomials; reproduces (3.26)). Naive T-series = CG / L(4s-3, alpha^4) on the c-axis (PROVED). The normalising L-factors for roots (1,0) [alpha^4] and (3,1) [alpha^12] sit at argument 1/2 at the residue centre X1=-1/8. Fix (a): convolution T = T* * lambda removes (1,0). The (3,1) factor remains -> X1 >= -1/12 -> C << Z^{5/12}. reduction_g2.py: (L3) << Q^{1.510} (fix a), Q^{1.531} (no fix), Q^{1.5} needs (3,1) handled. S51 §3's Q^{3/2} downgraded accordingly.
- 00:09Z S51 §21: leftover p-part E'=(1+Ax1+Bx2)/N has its (5,2) monomial (e=6) converging only for X1 > -1/10 at X2=1/4. So with (a) plus GRH(L(alpha^12)): C << Z^{2/5}, (L3) << Q^{1.506}. Route outputs: 1.531 / 1.510 / 1.506 / 1.5 (needs (5,2) matching). All < 13/8.
- 01:00Z S52 mean-value Bochner (l3f/): Bochner + envelope sup argument on mean-value regions (LoA in d, S49 mean squares in c, blockwise convexity), all 12 G2 chambers, DDHL normalisers (evaluate only at x,y>=0). LP bookkeeping through S50 with bootstrap: main box (L3) << Q^{1.5417} (1.5365 with fix (a)). Second-moment hull: no gain (m-FE gives 1/L). Caveat: near-main row boxes need the twisted-kappa version (theta<=1) - next.
- 01:07Z S52 §8-9: crude box reduction fails (kappa-exponent ~1.26); reflecting boxes through the c-FE (fix (a)) puts them strictly below the main box (d/dx2 ~ -2.8). S52 continuous in Y/Q^{3/2} (general_q.py). Full sum (L3) << Q^{1.5365+eps} modulo STANDARD bookkeeping. Large twists Z^{0.24}..Z^{1/3} remain the only loss before Q^{3/2}; Q^{3/2} is the Cauchy floor.
- 01:17Z S52 §10 correction: S50's X(h0,j) (coprime to j, N(j)<=(H/h0)^{1/4}) has effective twist h0*e^2 (e|j) with c-twist kappa=e; the S52 bounds are twist-sensitive, so e_X^eff(eta)=max_eps e_X(eta+2eps). Revised: (L3) << Q^{1.5625} (no fix a), Q^{1.5573} (fix a). Still < 13/8.
- 01:28Z S52 §11-12: block-optimal regions (B1 large sieve over c / B2 mean square over d with reflection + S(q) table), RMS over h0 (true quartic Xi), joint (tK,tE) exponents, recursion on S(q). Fixed point (L3) << Q^{1.5156}; s(1)=2.25 limits (j-pieces at q=1). Fine-grid verification running.
- 01:33Z S52 fine-grid verification: fixed point unchanged (L3) << Q^{1.5156}.
- 01:34Z e-tail cut: s(1) 2.2496 -> 2.2188 (one pass); continued recursion recur5 running.
- 01:36Z S52 §14 correction: generic twist range costs J*Xi4(K,Q) (>= (KQ)^{2/3}); at q=1 forces s(1)>=7/3. 1.5156 superseded; safe value 1.5253; recur6 running.
- 01:41Z S52 §15: corrected recursion iteration 0: s(2/3)=2.007 -> (L3) << Q^{1.5104} (direct check of X at binding twist OK). Q^{3/2} via this route needs s(1)<=13/6, blocked by balanced quartic large sieve (>= Y^{5/4} for general coefficients).
- 01:47Z S52 §16: block Region_d omitted the theta polar term in l_d=0 chambers (x>3/8+theta(1/4-y)); all block-region numbers (1.5234..1.5016) invalid pending recur7. Valid: 1.5573 (S52 §10).
- 02:18Z S52 §17: audit found 3 more pipeline errors in recur6/7 (bootstrap started at unproved lam=1.5253; np.interp clamp gave S(q>1.6)=2.933 instead of s49; reflected R_d points in polar chambers lacked polar exponent). None binds. recur9 (proved start 13/8): 1.5104, 1.5053, 1.5027, 1.5016 (floor = 1e-3 sampling offset); recur10 (offset 1e-4): 1.5008, 1.5004, 1.5002. §10 re-run from 13/8: 1.5573 (fix a) stands. Mechanism: Region_d boundary point (1/3,0) in chambers s_c and dcdc gives x_r=-1/12, b=(3x'-1/4)eta; X(1/3)=5/6+x'+b_loc; b_loc set by block D=C^{2/3} ((L3) shape) = 3eta((s-2)/2-delta'); contraction s_{n+1}-2=max(delta',(s_n-2)/2). Limit table S(q)=q+4/3 (q>=2/3) verified by fresh code (cone_check.py): b_loc=0 on whole cone, strict interior. Cocycle check of h0 factors: 12 chambers, 0 inconsistencies. Result: (L3) << Q^{3/2+eps} within the S52 framework (Cauchy floor); open: Lemma 4 constants, boxes/T<->T* transfer, rigorous-numerics pass. Interpolation: only step-up is justified; recur11 (step 0.01, step-up) running.
- 02:34Z S52 §18 / l3f/S52_LIMIT_PROOF.md: Q^{3/2+eps} made analytic (no LP, no grids, no interpolation). Explicit hull point: Region_d boundary (1/3+dp,0) pulled back from chambers s_c and dcdc (x_r=-1/12 on 3x+y=0, b=(3x'-1/4)(tK+tE)+b_loc). Cone lemma b_loc <= 3*Delta*(tK+tE)+tE/2: exact rational vertex enumeration at dp=Delta=0 (7323 vertices, max 0; cone_exact.py) + monotone in dp + dB2/dDelta <= 3(tK+tE) from sector geometry (tC<=9 on [5/9,2/3], tC<=3 on [2/3,5/6]); spot checks exact (max = max(-3Delta,-9dp)). Bookkeeping closed form at every q in [2/3,5/6]: Delta' = kappa/(5.5+12 kappa), kappa=dp+Delta; from S49's Delta0=1/12 -> Delta* <= dp/4.5 -> lambda <= 3/2+dp. Grid cross-check (limit_map.py): 1/12 -> 0.012917 (formula 0.012834) -> 0.002346 (0.002302). recur11 (step-up, h=0.01): 1.5055; recur12 (h=0.004): 1.5020 (cut by worker restart). Open: Lemma 4 constants, boxes / T<->T* transfer. (L3) at Q^{1+eps} still OPEN (needs non-Cauchy step); theta unchanged.
- 03:09Z Full sum (l3f/S52_FULLSUM.md): every row m=m1 m2^2 m3^3 m4^4 reduces (S49 corollary, Moebius over d|m4) to a twisted main-box problem T*theta with index-conductor F=N(m2)^2N(m3)N(d)^2. Cauchy step drops theta; F enters only via the c-reflection cocycle (path-independent, 7/12 at the two-chamber point; fcocycle.py) and the cone lemma (b^theta <= 3D(tK+tE)+tE/2+(5/6+kappa/2+2D)tF, via tD'=tD+tF/3 since D^3F=(DF^{1/3})^3; grid check twist_cone3.py gives 0.666). Twisted bookkeeping: excess (c t-2.75u)_+/5.5 above 2/3, Cauchy floor 2+(c t-2.75v)_+/5.5 below; c*<=1.57. Box sum exponents X2^{-1+.364c} X3^{-3+.18c} X4^{-5+.364c}: converges for c<2.75. => (L3) full sum << Q^{3/2+eps} given framework. Also l3f/BELOW_3_2.md: the (Z, Z^{4/3}] k-range is self-dual; S50/S52/GRH/large sieve all invariant; (L3) needs a duality-breaking input. Script bugs fixed on the way: NaN->max(0,nan)=0 in cone_g.py; tC'=0 edge in twist_cone3.py.
- 03:38Z CORRECTION (l3f/S52_CORRECTION_NAIVE_VS_CG.md): S52's non-identity chambers need the clean CG series. Exact CG p-part exponents (coef_5.txt): the binding point v2=(1/3,-1) (chamber dcdc) has worst E' exponent +1 and |Bx2|=q^{1/2}>1, while every chamber-dcdc Region_d point pulls back to y<-1/2, beyond any continuation of the naive series. So Lemma 4(f) ('verbatim') fails for the naive series. The argument holds for Z_CG with transfer naive=Z_CG*E'/norms at (-1/12,1/4) (E' abs. conv., worst -1.083) and dual points (x,-3x), -1/12<x<0. => Q^{3/2+eps} (and S52's 1.5573) are CONDITIONAL on task #26. Unconditional state: Q^{13/8} (S49). Next: #26 slice-property verification (CG p-part x1-polynomials vs Kubota local factors).
- 03:42Z #26 analysis: exact CG p-part shows the c-slices of CG = the naive fix-(a) series at squarefree d (H_{x2^1} = q^{1/2}), while the d-slices of CG include the d∋p^2 terms (H_{x1^1} = (A+q^2 x2^2)/(1-q^3x2^4)) that the naive coprime series lacks, giving its y=-1/4 natural boundary. Literature: Chen arXiv 2411.13143 proves the BBF Eisenstein conjecture (WMDS = Whittaker coefficients of metaplectic Eisenstein series) for adjoint metaplectic dual data (G2 included) but needs k ⊇ mu_{2n}; Q(i), n=4 lacks mu_8 ((pi,pi)_4=-1 at q≡5 mod 8 — exactly where S51 §19 needed the BD sign). #26 = Chen over k ⊇ mu_n with BD sign, plus the alpha and h0 twists.
- 04:35Z #26 local theory (quartic/l3g/S53_QI_LOCAL.md):
  * Derived CG's slice expansion over Q(i) keeping self-symbols. Lemma 5.4 is verbatim (signs s^{(i+1)(4-i)} = 1). Lemma 5.9 gains
    S(k') = s^{(K'+1)(k' + sum_{j>i} c(j,i) beta_j)}.
  * Rank-1 exact check (rank1sym.py): true Gauss-sum local factors satisfy the local FE for l = 0..12 only with S.
  * The matching local action is CG Def 3.1 with raw Gauss sums and Q multiplied by s^{mu(B_i+1)}. FLINT exact algebra (relflint.py,
    thm41sym.py; sqrt q and gamma(1) indeterminate):
    * W-relations hold identically for l in {0..7}^2 (all l by x^{-4 omega_k} conjugation);
    * N(x;l) is polynomial;
    * Thm 4.1' holds in 4191/4191 classes.
  * General l via the CG Thm 3.6 identities (I1), (I2), whose signs cancel since t is odd.
  * S51 §19's BD-signed action: a W-action, but it fails 343/554 slice classes. Its p-part differs from the correct one by monomial
    signs only, so the E' exponents are unchanged.
  * alpha-twist consistent (k-independent alpha(m_i A)^l).
  * Remaining for #26: CG Prop 5.7 Psi-bookkeeping over Q(i) (self-symbols are (y,-1)_S, Omega-trivial). Labelled STANDARD.
  * (L3) status unchanged: Q^{13/8} unconditional; Q^{3/2+eps} conditional on that bookkeeping; Q^{1+eps} OPEN.
- 04:55Z #26 global step (S53 §4):
  * CG Prop 5.6/5.7 over Q(i) proved, given Kubota's rank-1 FE in Brubaker–Bump form (T_b phi = sum_eta P_{b eta} phi~_eta):
    * equivariance in a_i: T_{b y^-2}(phi_y) = (T_b phi)_y, the only new symbol being (y,-1)_S^{2t} = 1;
    * Omega-covariance in a_j: new symbol (u,-1)_S^t = 1 for u in Omega.
  * So every self-symbol from mu8 ⊄ Q(i) is accounted for.
  * #26 is now CONDITIONAL on one rank-1 input: [BB06 Thm 1] for Q(i), n=4, with infinity-type twist. BB06 assume mu_{2n}; DIPP
    (arXiv 2607.27131, Rmk 1.4.1/4.1.8) call removing it routine; Dia04 (Invent. Math. 157, F ⊇ mu_r) / DDHL (4.28) give the FE in
    cusp-vector form. Converting cusp-vector form to T_b form is a finite 2-adic step, not done here.
- 05:30Z #26 global step closed (S53 §4 rewritten). (E1) no longer assumes the Brubaker–Bump form. It follows from consistency of
  DDHL's FE with the ϖ⁴ slice identity:
  * D(s,u0ϖ^4;φ) = (1+q^4x^4)D(s,u0;φ) − g(u0,ϖ)x D(s,u0ϖ^2;θ1φ), checked coefficient by coefficient by hand. The ϖ^4 coefficient
    q^4 − q^3 needs the self-symbol.
  * Injectivity via g(1,π)^4 = π^3 π̄, which makes the Gauss-sum ratio a Grössencharacter, plus Hecke equidistribution.
  * (R1) class dependence read from DDHL Remark 4.2 / (4.19).
  * #26 status: PROVED modulo DDHL (4.28) + Remark 4.2 (published) and STANDARD Bochner. So (L3) ≪ Q^{3/2+eps} modulo S52 Lemma 4
    constants. Q^{1+eps} still OPEN.
- 05:50Z S53 §6: S52 Lemma 4 checked for the CG series at the points the proof uses. Generic-prime corrections converge (worst
  -1.333 at (1/3,0), from (1,2)). Twisted local factors: e(k) <= k1+k2 for all l in {0..7}^2 and s = ±1, so they are O_l(1) at
  x',y' >= 0 (h0 is fourth-power free). Bochner needs no connectivity because the CG function is entire. Q^{3/2+eps} chain now open only
  for an end-to-end audit.
- 06:05Z BELOW_3_2 §7: the hyperbolic (2,6,6) second-moment series. The Tits cone contains the (L3) point, but Bochner from absolute
  convergence reaches the central line only at Re w >= 5/2 (LP over 5634 group elements, stable). That is convexity, so this route is
  closed quantitatively.
- 09:40Z (second 5-hour L3 block, l3j/). THEOREM_ZF: any (L3)-type bound Q^{3/2-delta}, uniform in index twists, implies that all
  quartic Hecke L(s,(k0/.)_4 lambda^e) over Q(i) are zero-free for Re s > (7+10θ')/(8+8θ'). (L3) itself gives 7/8.
  * Mechanism: repeated dual frequencies k0·ϖ^4, positivity from a narrow-annulus expansion, optimised row scale.
  * Checked numerically on the mu*T control family (excess 0.06/0.16/0.19 vs predicted 0.10/0.17/0.22 at Z = 500/2000/4000).
  * Corrections: the NOTES §8G corollary (positivity) and the BELOW_3_2 §7 model (wrong FE group). gen_rows' "untwisted" control is
    alpha-bar*T.
  * (L3) at Q^{1+eps}: still OPEN. A proof must contain a new quasi-RH for the quartic family plus beyond-GRH control of the 4th-power-free
    frequencies.
