# (L³) attempt log — session started 2026-10-15T08:51Z

Target. (G) at M = 2/3: for Q = Z^{2/3},
  S(Q,Z) = Σ_{N(m)≤Q, m primary} |Σ_{c sqfree primary, N(c)~Z} T(c) χ̄_c(m) W(N(c)/Z)|² ≪ Q·Z^{1+ε},
  T(c) = α(c) g̃(χ̄_c),  α(c) = c/|c|,  χ_c = (·/c)_4.

## Verified identities (NUMERICAL, exact to 1e-14)
1. g̃(χ̄_c)² = μ(c)·ᾱ(c) for all 248 squarefree primary c, N(c) ≤ 700. Hence T(c)² = μ(c)α(c).
2. T(ab) = T(a)T(b)·(a/b)₂·s(a,b), s(a,b) = −1 iff N(a) ≡ N(b) ≡ 5 (mod 8); 118/118 random coprime prime pairs.

## Exact reductions (derived; algebra checked)
- Poisson in m maps (G) at M = 2/3 to the Möbius mean square
    Σ_{N(k)≤H} |Σ_{N(c)~Z} μ(c) χ_k(c) W/√N(c)|² = H·C + O(Z),  H = Z^{4/3},
  using identity 1 (α cancels against g̃²). Poisson in k maps back: fixed point.
- AFE for the inner sum works for Hecke-character coefficients (short dual of length H^{1/4}, exact
  orthogonality over moduli ≤ H^{1/2}) but the Möbius coefficients give 1/L(s, χ_k): no functional equation;
  contour shifts cross zeros of L(s, χ_k).
- Bilinear (Type-II) route via identity 2: pointwise in m capped by convexity (trivial at Z = Q^{3/2});
  Cauchy–Schwarz to use the m-average loses a factor B. Needs joint square-root cancellation.
- Cube-index reformulation: ψ_α(m³,s) = twist of K_α(s) = Σ α(c) g̃(χ̄_c) N(c)^{-s} by χ̄_m (reciprocity);
  all-index Lindelöf on average (de Faveri–Dunn–Hoffstein) restricted to cubes gives Q³; large-sieve cube
  detection gives the Q^{3/2} barrier; Gallagher's larger sieve needs Σ_cubes |ψ|⁴ (circular).

## Session conclusions (2026-10-15, 08:51–11:51Z) — read this first
Status of (L³): NOT proved, NOT disproved.
- Evidence it is TRUE (NUMERICAL, S1/S11/S36): at Q = Z^{2/3}, Z = 1e5, 3e5, 1e6, 3e6, 6e6: S/diag = 1.005, 1.023, 0.996, 1.042, 1.010
  with Exp(1)-distributed rows; no growth beyond sampling noise.
- What it is (PROVED reformulations, S6/S9/S10): (G) ⟺ square-root cancellation, on average over N(m) ≤ Q, of the Möbius-weighted
  Kubota sums Σ_{c~Z} μ(c)g̃(χ̄_c, m), Z = N(m)^{3/2}; equivalently a Lindelöf-on-average bound for Q twists of a degree-3 metaplectic
  object by the quartic characters (m/·)_4 (conductor Q³), i.e. the analogue of Σ_{χ mod q}|L(½,π⊗χ)|² ≪ q^{1+ε} for GL3 π, which is open
  even for π = 1⊞1⊞1 (sixth moment, single modulus). It is a parity-sensitive statement.
- Where it sits (S4/S12): implied by the ratios conjecture for the quartic Hecke family; not known to follow from GRH; it implies quasi-RH 7/8 for
  that family and Weyl-strength subconvexity for each ψ_α(m³,·).
- Why every route closes (S3/S7/S8): Poisson duality is Fourier inversion at M = 2/3 (fixed point); Kubota reflection is self-dual
  there; no Voronoi formula in the modulus exists for quartic Gauss sums (no theta identification); Type-II splitting creates
  index-shifted rows of larger conductor; large-sieve shape is not even known (Gao–Zhao loss, dFDH non-orthogonality).
- Even proved, (L³) alone would not give 5/6: the quartic analogues of the marked inverse moment etc. are also needed (S8).
- Conditional proof written (S18, L3_REPORT.md §10): a ratios-type hypothesis (R) for the quartic Hecke family implies (G) for all
  M ∈ (0,1] with a power saving. Weakest sufficient form (S20): error H^{3/4 − 3σ0/2} at Re s = ½ + σ0.
- Pipeline validated against brute force (S17, S19); identity T² = μα checked for all 347 518 squarefree c ≤ 1e6.
- Bias test over 36 windows up to Z = 2.5e6: S/diag = 1.008 ± 0.007 (S15/S31); single run at Z = 6e6: 1.010 ± 0.013 (S36);
  M from 0.55 to 1 ≈ 1 (S14/S22); t ≤ 3000 uniform (S16); weights (sharp, wide) ≈ 1 (S30); 4th/6th moments Gaussian (S23).
- New obstruction for the quartic route (S32): with PRIME columns the α-twisted mean square is biased and growing, 1.23 → 1.29 → 1.44 → 1.67
  for Z = 1e5 → 3e6, ∝ Z^{1/2}/log Z (theta terms at primes, T(π) = −conj g̃), concentrated on N(m) ≡ 1 mod 8; quartic prime-slot moments
  would need explicit quartic theta coefficients.
- Side observation (S29): parity-blind theta residues at prime index vanish iff π ≡ 3+2i mod 4, with a size split by p = x² + 64y²
  (2-adic local behaviour, probably known; not needed for (L³)).
- Theorem Z re-derived in outline (S27).
- Side results this session: B headroom in the combined proof is ~4e-6 (negligible θ gain); local commit rebased onto upstream
  24fc1ccd (PR #53 merged, 0.87486648, not lower), all tracker checks re-run and passing.

## Session work
### S1. Scaling of (G) at M = 2/3 (l3/gen_rows.c, l3/rowstats.py) — NUMERICAL
| Z | Q=Z^{2/3} | rows | S_tw/diag | S_untw/diag | per-row var (Exp=1) | max row |
|---|---|---|---|---|---|---|
| 1e5 | 2154 | 847 | 1.005 | 0.996 | 1.09 | 6.3 |
| 3e5 | 4481 | 1759 | 1.023 | 1.088 | 0.94 | 7.3 |
| 1e6 | 9999 | 3928 | 0.996 | 1.008 | 1.01 | 7.6 |
| 3e6 | 20800 | 8172 | 1.042 | 1.064 | 1.13 | 13.1 |
| 6e6 | 33019 | 12961 | 1.010 | 1.045 | 1.06 | 9.1 |
Twisted rows |R(m)|²/rowdiag are Exp(1)-distributed (P(>3)=0.0502 vs 0.0498, P(>5)=0.0061 vs 0.0067 at Z=1e6):
consistent with (L³) AND with the much stronger pointwise statement max_m |ψ_α(m³,½)|² ≪ log Q.
Untwisted control: bias rows at N(m) = 17³, 29², 13² (theta residue at 4th-power-related indices), absent when twisted by α.
No growth in Z: no numerical disproof; (G) appears true with constant ≈ 1 (asymptotic S ~ diag).

### S2. Prime signs σ(π) (l3/primes_g.c, l3/sigtest.py, l3/sigtest2.py) — NUMERICAL + PROVED identity
- Computed g(χ̄_π) for all 78350 split primary primes, N(π) ≤ 1e6. Checks: T(π)² = −α(π) (err 6e-12); G(π̄) = conj G(π) (exact, proved
  by x ↦ x̄), so only Im π > 0 carries information.
- PROVED: at primes α(π)g̃(χ̄_π) = −conj(g̃(χ̄_π)) exactly (T² = −α and |g̃| = 1 give g̃ = ±i·ᾱ^{1/2}; checked on all 78350 primes,
  sign −1 in every class). Hence the α-twisted prime sum is minus the conjugate untwisted one (observed: k=0 and k=1 angular sums are exact negatives, 1078.2 vs −1078.2).
- The prime sums carry the quartic-theta (Patterson-type) bias: class π ≡ 1 (mod 4): Σ g̃ = 1838 (z = 9.3), class 3+2i: −760 (z = 3.8);
  scale X^{3/4}/log X = 2289 at X = 1e6. So σ(π) is NOT equidistributed in residue classes mod 4 at scale X^{3/4}; no Hecke-character
  correlation beyond that (891 tests, all other |z| < 4.3, the largest being sub-classes of the mod-4 bias).
- Consequence for (L³): the "double cover" decomposition T = σ·√(μα) carries no exploitable structure beyond the theta bias,
  and the bias disappears for the Möbius-weighted composite sums (μ kills the pole as for 1/ζ), consistent with S1.
  The sign-correlation route is CLOSED (negative).

### S3. The duality is a tautology at M = 2/3 (derivation)
Poisson in m, then CRT g(χ̄_cχ_{c'},k) = s(c,c')χ̄_c(c')χ_{c'}(c)g(χ̄_c,k)g(χ_{c'},k), then
  T(c)g(χ̄_c,k) = μ(c)χ_c(k)√N(c),   conj T(c')·g(χ_{c'},k) = μ(c')χ_{c'}(−1)χ̄_{c'}(k)√N(c')   (both PROVED from T² = μα, conj g(χ̄) = χ(−1)g(χ))
turns Off_G(Q,Z) into (Q/Z)·Σ_{classes mod 8} ±Off_Möb(H = Z²/Q, D = Z). The second Poisson (in k) returns Off_G(Q', Z) with Q' = Z²/H = Q.
The composite map is Fourier inversion: X = X. No contraction exists; any proof must inject information not visible to Poisson + CRT + (I1).
Also: the w-functional equation of D(s,w) = Σ_c T(c)L(w,χ̄_c)N(c)^{-s} has root number ε(χ̄_c) ∝ g̃(χ̄_c), so T(c)ε(χ̄_c) = μ(c):
the first-moment double series is mapped to Σ_m N(m)^{w-1}/L(s',(·/m)_4) (no FE). Only (Z/2)² acts — same obstruction as F's (2,6,6).

### S4. Is (L³) implied by GRH? (derivation) — not by any known argument
Under GRH for the quartic Hecke family, pointwise |Σ_{c~D} μ(c)χ_k(c)/√c| ≪ H^ε, giving Σ_k |M(χ_k)|² ≪ H^{1+ε}. But (G) needs the
ASYMPTOTIC H·C + O(D), D = H^{3/4}. Using GRH on each off-diagonal character sum Σ_{k≤H} ψ_{cc'}(k) ≪ H^{1/2+ε} gives D·H^{1/2} = H^{5/4},
short by H^{1/2}: square-root cancellation over the pair (c,c') is still required. ψ_α(m³,s) has no Euler product (metaplectic
Eisenstein coefficient), so GRH says nothing about it. (L³) is orthogonal to GRH; it is the metaplectic analogue of the sixth
moment Σ_{d≤Q}|L(½,χ_d)|⁶ ≪ Q^{1+ε} (GL3 object, family size = conductor^{1/3}), which is open even under GRH-free methods;
Young's third moment Σ_d L(½,χ_d)³ (first moment) is the solved analogue. The first-moment analogue here, Σ_m ψ_α(m³,½), is also
blocked: the w-functional equation sends T(c) to μ(c) (S3).

### S5. Literature check (web search, 2026-10-15)
- de Faveri–Dunn–Hoffstein arXiv 2607.07911: Lindelöf on average over ALL indices (already used; gives Q³ on cubes).
- Diaconu–Ion–Paşol–Popa arXiv 2607.27131: first/second moments of r-th order Hecke L-functions over fields with μ_{2r}
  (r ≥ 3) via MDS — standard families (size = conductor), not thin; for r = 4 needs Q(ζ8), not Q(i).
- Hamdar arXiv 2509.06152: notes the Gauss-sum Dirichlet series is understood only for cubic sums, not order ≥ 4.
- Goel–Ray arXiv 2505.12015: second moment of cubic Dirichlet L-functions over F_q(t).
None addresses a thin (cube-index) family or a Voronoi formula in the modulus c for quartic Gauss sums.
- GL3 comparison check (web search): the sixth moment of Dirichlet L-functions has an asymptotic only on average over moduli
  q ≤ Q (Huxley's upper bound; Chandee–Li–Matomäki–Radziwiłł arXiv 2409.01457 asymptotic). No single-modulus q^{1+ε} bound was
  found, consistent with labelling the GL3 analogue of (L³) OPEN.

### S6. Elliptic-curve form of the rows (PROVED identity) and the GL3 comparison
Grössencharacter of E_D : y² = x³ − Dx is ψ_D(π) = conj((D/π)_4)·π for primary π (Ireland–Rosen ch. 18). Hence for c coprime to 2m
  α(c) χ̄_c(m) = ψ_{E_m}(c)/|c| =: ψ̃_{E_m}(c),   so   R(m) = Σ_c g̃(χ̄_c) ψ̃_{E_m}(c) W(N(c)/Z),
and (T(c)χ̄_c(m))² = μ(c)·ψ̃_{E^{(m)}}(c) with E^{(m)} : y² = x³ − m²x (quadratic twist of the congruent-number curve).
So ψ_α(m³,s) = "Kubota series at index 1 twisted by the Hecke character of the quartic twist E_m". A's degree count
(conductor N(m)³ for a twist of conductor N(m)) says the Kubota series behaves like a DEGREE-3 object over Q(i), so
  (L³)  ≈  Σ_{N(m)≤Q} |L(½, Φ ⊗ (m/·)_4)|² ≪ Q^{1+ε},  Φ of degree 3, Q twists, conductor Q³.
Genuine-GL3 analogue: Σ_{χ mod q} |L(½, π⊗χ)|² ≪ q^{1+ε} (q members, conductor q³). This is OPEN even for π = 1⊞1⊞1
(the sixth moment of Dirichlet L-functions for a single modulus; known q^{4/3+ε} via Weyl/Petrow–Young + fourth moment;
Conrey–Iwaniec–Soundararajan need the extra average over q ≤ Q, family Q² ≥ length Q^{3/2}). In the metaplectic case even the
Weyl input (pointwise subconvexity for ψ(m³,½)) is unknown: no Voronoi formula in the modulus c exists for quartic Gauss sums.

### S7. Minimal breakthrough that would lower θ via the quartic route (derived from NOTES §3/§6 bookkeeping)
With (G) at M = 1 − ℓ holding up to loss Z^λ: σ0 ≈ 7/8 − 3ℓ/8 + λ/2 (HEURISTIC: linear interpolation between the two NOTES
endpoints λ = 0 → 7/8 − 3ℓ/8 and λ = ℓ (large-sieve shape) → 7/8 + ℓ/8, at g = 0).
Writing λ = ℓ − δ: gain = δ/2 − ℓ/8. To beat 0.8748389 need δ/2 − ℓ/8 > 1.61e-4, e.g. ℓ = 0.001, δ > 5.7e-4:
  Σ_{N(m) ≤ Z^{0.999}} |R(m)|² ≪ Z^{2 − 0.00057}   (beat (Q+Z)Z by a tiny power when rows ≈ columns).
Status of the ingredients:
  - (Q+Z)Z itself is NOT known: Gao–Zhao give (Q+Z+(QZ)^{3/4}) (loss Z^{1/2} at Q≈Z); de Faveri–Dunn–Hoffstein show the
    quartic large sieve is NOT perfectly orthogonal for general coefficients (Gauss-sum bias), so the optimal shape is false in
    general and could only hold for special coefficients such as T(c).
  - Pointwise, no bound beyond trivial |R(m)| ≤ Z exists at M ≈ 1 (dual length Q³/Z = Z^{2−3ℓ} > Z).
(In PROBLEM.md's parametrization σ0 = 1 − M/4 the same requirement sits at M slightly above 1/2; the two are exchanged by the
reflection M ↔ M/(3M−1).) So even the smallest θ gain from the quartic route needs (i) an optimal-shape quartic large sieve for T-coefficients and (ii) a power
saving beyond it. Both are open; (L³) is the strongest form.

### S8. (L³) is necessary but not sufficient for 5/6
The OpenAI sextic proof (openai-baseline Lean tree) rests on cubic-theta machinery: Eisenstein/ (Kubota continuation via
meromorphic resolvent/Fredholm), Reflection/ (266 modules), Descent/ (419), CubicSieve/ (91), Cusp/ResiduesModThree, RowCompletion/.
Their inputs: plain moments, marked inverse moments with prime slots (r + 2z ≤ m, 2r + 8z ≤ 3m), the canonical moment, the Gram low
side, completed-row energy. Each uses Patterson's identification "cubic Gauss sums = cubic theta coefficients" to reflect μ-weighted
or Gauss-weighted sums. The quartic route needs quartic analogues of ALL of these (NOTES §4). (L³) is the canonical-moment analogue;
the marked inverse moment analogue needs reflection of μ-weighted prime-slot products, which by S3 is again an (L³)-type fixed-point
statement. Kubota's continuation and the Gao–Zhao quartic large sieve do exist for the 4-fold cover; the theta identification does not.
So a proof of (L³) would be the main step, but turning it into θ = 5/6 means a new formalization on the scale of the sextic one.

### S9. Cleanest form: Möbius-weighted Kubota sums at index m (PROVED reformulation)
From T(c)² = μ(c)α(c) and |T| = 1:  T(c) = μ(c)·conj(g̃(χ̄_c)) for every squarefree primary c (at primes: T(π) = −conj g̃(χ̄_π), S2).
With g(χ̄_c, h) = χ_c(h)g(χ̄_c):   R(m) = conj( Σ_c μ(c) g̃(χ̄_c, m) W(N(c)/Z) ).
So (G) ⟺ Σ_{N(m)≤Q} |Σ_{c~Z} μ(c) g̃(χ̄_c, m) W|² ≪ Q Z^{1+ε},  Z = Q^{3/2}.
Contrast: WITHOUT μ, Σ_c g̃(χ̄_c, m) W(N(c)/Z) is a Kubota sum of conductor N(m) ≤ Q ≪ Z, i.e. far past its AFE length; it equals
the quartic-theta residue term plus a negligible dual. The Möbius weight raises the conductor from N(m) to N(m)³ (via I1),
which is exactly what makes the problem thin. Known pointwise technology for μ-weighted Gauss-sum sums (Heath-Brown–Patterson
Vaughan-identity + cubic/quartic large sieve type, X^{5/6}-strength) gives S ≪ Q Z^{5/3}: short by Z^{2/3}.

### S10. Parity formulation (derivation)
Let K_z(m) = Σ_{c sqfree, N(c)~Z} z^{ω(c)} g̃(χ̄_c, m) W.  z = −1: by S9 this is conj R(m), the (L³) row.  z = +1: the plain
Kubota sum at index m (parity-blind). CAUTION (numerics, l3/gen_k.c, S11): K_1(m) is NOT "theta residue + negligible dual" in the
squarefree-restricted form used here; the restriction to squarefree c and coprimality to m introduce 1/ζ_K-type factors, and
observed |K_1(m)| for small m grows roughly like √Z to Z^{3/4} (bias rows), not as a clean residue term.
With K_± = sums over ω(c) even/odd:  R̄(m) = K_+ − K_− = 2K_+(m) − K_1(m).  Hence
  (G) ⟺ Σ_m |K_+(m) − ½K_1(m)|² ≪ Q Z^{1+ε}:
the even-ω part must carry exactly half of the theta main term, up to square-root error, on average over m.
This is a parity-sensitive statement (Selberg's parity barrier): sieve or Type-I information sees only K_1, and the content of (L³)
is Type-II (bilinear) information with square-root strength. The twisted multiplicativity
  g̃(χ̄_{ab}, m) = s(a,b)·g̃(χ̄_a, m b²)·g̃(χ̄_b, m)
turns every Type-II piece into a sum over a of index-shifted rows (index m b²), i.e. rows of conductor N(m)³N(b)², longer still.

### S11. Z = 3e6 and arithmetic breakdown (l3/rows_Z3000000.bin) — NUMERICAL
- Ratio 1.042 (rows come in conjugate pairs with |R(m̄)| = |R(m)|, PROVED via T(c̄) = conj T(c); effective n = 4086, z ≈ 2.6).
  Tails: P(>3) = 0.057 (Exp 0.050), P(>5) = 0.0097 (Exp 0.0067); max 13.1 at N(m) = 5777.
- Breakdown (mean ± s.e.): N(m) prime 1.011±0.029; composite 1.057±0.025; non-squarefree norm 1.051±0.032;
  N(m) ≡ 1 mod 16: 1.072±0.036 (1.090±0.051 at Z = 1e6) — at most a 2σ effect, no clear structure.
- Small-m band (N(m) ≤ 0.1Q): 1.091 (n = 815). For these rows Z ≫ N(m)^{3/2}, so reflection makes them short dual sums
  (length N(m)³/Z), with slightly heavier tails.
- Index-m rows (gen_k.c, Q = 60): the Möbius rows K_{−1}(m) stay at diagonal size (ratio 0.73, 0.96, 0.99 for Z = 3e4, 1e5, 3e5);
  the parity-blind rows K_1(m) grow (ratio 2.5, 2.7, 4.9), led by N(m) = 17 at ≈ 3× diagonal size. The parity-blind sum carries bias;
  the Möbius sum does not.
- Z = 6e6 run started 09:35Z.

### S12. Where (L³) sits among standard conjectures (derivation)
Möbius mean square by Mellin:  Σ_{N(k)≤H} |Σ_c μ(c)χ_k(c)N(c)^{-1/2}W(N(c)/D)|²
   = (2πi)^{-2} ∫∫ W̃(u)W̃(v) D^{u+v} Σ_k [L(½+u, χ_k) L(½+v, χ̄_k)]^{-1} du dv   (up to local factors at 2 and coprimality).
- The Conrey–Farmer–Zirnbauer RATIOS CONJECTURE for {L(s,(k/·)_4) : N(k) ≤ H} (denominators only, shifts with Re u, Re v ≥ ε)
  gives Σ_k [·]^{-1} = H·A(u,v) + O(H^{1/2+ε}). Integrating at Re u = Re v = ε: Möbius MS = H·C + O(H^{1/2+ε}D^{2ε}).
  (G) at M needs error O(D) with H = D^{2−M}, i.e. H^{1/(2−M)} ≥ H^{1/2}: true for every M ∈ (0,1].
  So RATIOS CONJECTURE ⟹ (G) for all M ∈ (0,1] ⟹ (L³).
- GRH alone does not suffice by the natural argument (S4): pointwise 1/L bounds lose H^{1/2} on the off-diagonal (no implication known either way).
- Conversely (G)/(L³) ⟹ the zero-free half-plane Re s > 7/8 + θ'/2 for the whole quartic Hecke family over Q(i)
  (Theorem Z, subagent G), and Weyl-strength subconvexity for each ψ_α(m³,·).
Chain:  Ratios(quartic family)  ⟹  (L³)  ⟹  quasi-RH 7/8 for the quartic family + Weyl subconvexity for ψ_α(m³,·).

### S13. Small-index rows and the normalizing L-factor (l3/gen_T.c dump of all T(c), N(c) ≤ 1e6; l3/rz.py) — NUMERICAL + HEURISTIC
- Rows with N(m)³ ≪ Z (reflected length < 1, e.g. m = 1 at Z = 3e6: |R|²/diag = 3.8) are NOT small. Over 26 log-spaced Z in [3e3, 1e6]:
  m = 1: mean (|R|/Z^{3/8})² is 0.20 on [3e3, 3e4] and 0.24 on [1.6e5, 1e6] (stable), while |R|²/diag drops from 0.74 to 0.36.
  m = π5: (|R|/Z^{3/8})² goes 0.21 → 0.64. Noisy, single realizations.
- HEURISTIC explanation: as in Σ_c c_c(n)c^{-2s} = σ_{1−2s}(n)/ζ(2s), the Kubota series is (nice Whittaker coefficient)/(normalizing
  Hecke L-function), here of the shape L(4s'−1, α⁴). Its zeros give terms Res·Z^{(1+ρ)/4}: Z^{3/8} under GRH, o(Z^{1/2}) unconditionally.
  For small m these dominate (size Z^{3/8}N(m)^{3/16}); at the top rows N(m) ≈ Q they are √Z-sized and indistinguishable from the
  random part. Bookkeeping under the random model reproduces S ≈ QZ. Consistent with the conductor-N(m)³ picture; no new leverage.

### S14. Scan over the row exponent M at Z = 3e5 (l3/mscan3e5.txt) — NUMERICAL
| M | Q | rows | S_tw/diag | S_untw/diag |
|---|---|---|---|---|
| 0.55 | 1029 | 406 | 1.016 | 1.136 |
| 0.60 | 1933 | 758 | 1.050 | 1.152 |
| 2/3 | 4481 | 1759 | 1.023 | 1.088 |
| 0.75 | 12818 | 5037 | 1.014 | 1.033 |
| 0.85 | 45243 | 17758 | 1.020 | 0.979 |
| 1 (square case) | 300000 | 117802 | 1.006 | 1.003 |
(G) holds numerically across the whole open middle range (1/2, 1), including the reflection pair 0.6 ↔ 0.75.
Twisted ratios sit slightly above 1 (+1.4% to +5%; standard errors 0.011–0.05 with paired rows): possibly a small positive lower-order term,
far too small to threaten (G).

### S15. Bias test over 30 dyadic windows (l3/gen_multi.c, l3/multi1e6.txt) — NUMERICAL
Three runs (Zmax = 1e6, 7.94e5, 6.30e5), 10 disjoint dyadic windows each, Q = Z^{2/3} per window, Z from 1.2e3 to 1e6.
Row-weighted mean of S_tw/diag = 1.0059 ± 0.0088 (all windows); 1.0081 ± 0.0099 (Z ≥ 1e5). The windows from different runs overlap,
so the true s.e. is up to √3 larger. Hence NO systematic bias: the +1–5% values in S1/S14 were sampling noise. (G) holds with
constant 1.00 ± 0.01 asymptotically, i.e. S ~ diag.

### S16. t-aspect uniformity (l3/tasp.py on the T(c) dump, Z = 1.5e5, Q = Z^{2/3}, 1107 rows) — NUMERICAL
Rows twisted by N(c)^{-it}: S/diag = 0.995 (t=0), 0.995 (3), 0.964 (30), 1.000 (300), 1.058 (3000). Uniform in t over this range,
as (L³) with the (1+|t|)^A factor predicts. (Once |t| is comparable to Z the archimedean conductor dominates and this would change.)

### Side result: where the next PROVEN θ gain would have to come from (combo/B_headroom.txt)
Combined endpoint model: dθ/dR ≈ 0.135 for the detector count at the binding point (δ,p) = (0.467, 1/2). Per unit saving δ·η in
each detector branch: inverse −0.043, long −0.020, affine −0.011. The inverse branch rests on the marked inverse moment, whose
capacity r + 2z ≤ 1 is the theta-reflection limit (§8 of NOTES), so the most sensitive lever needs a beyond-theta input.
B is capped twice over (high reserve and the 37/1000 slot bound).

### S17. Independent validation of the C pipeline — NUMERICAL (exact to 4e-13)
- T(c) from gen_T (twisted multiplicativity over prime Gauss sums) vs brute-force gauss4.gauss(c, j=−1) for ALL 210 squarefree
  primary c with N(c) ≤ 600 (91 norms, composites included): 0 mismatches, max error 3.7e-13.
- Quartic symbol exponents e_k = log_i (m_k/c)_4 used to build rows vs brute-force gauss4.chi: 1656 checks for N(c) ≤ 400,
  0 mismatches. So all row data in S1/S11/S14/S15/S16 use the correct T(c) and χ̄_c(m).

### S18. Conditional theorem written out (L3_REPORT.md §10)
Hypothesis (R): ratios-type asymptotic for Σ_k V(N(k)/H)[L(½+u,ξ_{k,ν})L(½+v,ξ̄_{k,ν})]^{−1} with error H^{1/2+ε}(1+|u|+|v|)^B,
uniformly for ε ≤ Re u, Re v ≤ 1/4 (finitely many ν mod 16).
Theorem: (R) ⟹ S_w(Q,Z) = D(Q,Z) + O(Q·Z^{(2−M)/2+ε}) for 0 < M ≤ 1; at M = 2/3 error QZ^{2/3+ε}, so (G) and (L³) hold with a power
saving. Proof: open the square → Poisson in m (k = 0 vanishes off-diagonal) → CRT + T·g = μχ√N identities → Möbius off-diagonal with
ν-twists → Mellin in c, c′ → (R) → the main term is exactly the diagonal → error H^{1/2}. Standard steps (Poisson on the primary
coset, imprimitive pairs, separation of variables) are noted, not written out. (R) is UNPROVED; it is the single missing input.

### S19. Numerical checks of the conditional proof's algebra
- T(c)·g(χ̄_c,k) = μ(c)χ_c(k)√N(c) and conj T(c)·g(χ_c,k) = μ(c)χ_c(−1)χ̄_c(k)√N(c): 160 random (c,k), N(c) ≤ 500, error < 7e-12.
- CRT: g(χ̄_cχ_{c′},k) = χ̄_c(c′)χ_{c′}(c)·g(χ̄_c,k)·g(χ_{c′},k) exactly (41 coprime pairs; the factor s in S3 is 1).

### S20. Weakest sufficient ratios hypothesis (trade-off curve; derivation)
If the family asymptotic in (R) holds only for Re u, Re v ≥ σ0 with error E(H), integrate at Re = σ0. The diagonal main term is
analytic, so shifting it is free. The error contributes H·(E/H)·Z^{2σ0}·(normalisation), i.e. Off_M ≪ E(H)·Z^{2σ0}. (G) at M = 2/3 needs
Off_M ≪ Z = H^{3/4}, so it suffices that
   E(H) ≪ H^{3/4 − 3σ0/2}   at   Re s = ½ + σ0.
σ0 = 0: error H^{3/4} (exactly the Möbius-mean-square form of (G)); σ0 = 1/6: H^{1/2} (the CFZ-strength error at Re s = 2/3);
σ0 → 1/2: O(1) (impossible). The curve shows there is no "cheap corner": moving the shifts right buys exactly what it costs in the
needed error. Under GRH, negative moments of L near the line are themselves open (cf. Gonek's conjecture for 1/ζ), so even
GRH-conditional routes to (R_σ0) are not available.
- Coset Poisson factorization: Σ_{x mod q(2+2i), x≡1 (2+2i)} ψ(x)e(Im(kx/(q(2+2i)))) = ψ̄(u_q)·g(ψ,k)·e(Im(k/(2+2i))),
  u_q = (2+2i)^{-1} mod q: 48 brute-force cases (8 moduli × 6 k), error 7e-14.

### S22. Square case M = 1 (Q = Z = 3e5, 117 802 rows, l3/mscan_M1.txt) — NUMERICAL
S_tw/diag = 1.006 (s.e. ≈ 0.004 with paired rows), S_untw/diag = 1.003. So the optimal large-sieve SHAPE, and even the diagonal
asymptotic, holds numerically for T-coefficients at Q = Z. This is the input Theorem Z needs at η = 1. Gao–Zhao's proven bound is
(Q + Z + (QZ)^{3/4}) ≈ Z^{3/2}, versus the observed ≈ Z: a Z^{1/2} gap between what is proved and what is true at the square case.

### S23. Higher row moments (rows_Z*.bin) — NUMERICAL
With x = |R(m)|²/rowdiag: E[x²]/E[x]² = 2.07, 1.89, 2.02, 2.04 and E[x³]/E[x]³ = 6.35, 5.18, 6.10, 6.47 at Z = 1e5, 3e5, 1e6, 3e6
(complex Gaussian: 2 and 6). So the 4th and 6th moments over the thin family also match the Gaussian model
(E[x^k] = k!), with no sign of extra structure at any moment order tested.
- Correction to S18: the reciprocity sign s(c,c′) = (−1)^{((Nc−1)/4)((Nc′−1)/4)} is not a product; it equals ½(1+λ(c)+λ(c′)−λ(c)λ(c′))
  with λ(c) = (−1)^{(Nc−1)/4} (checked on all 4 classes). So the off-diagonal is a finite sum of CROSS terms, and (R) must be assumed
  for pairs (ν, ν′) of class twists: Σ_k V/(L(½+u,ξ_{k,ν})L(½+v,conj ξ_{k,ν′})). Report §10 updated accordingly.

### S24. Squaring the rows: Rédei-type weights (derivation; structural remark)
For coprime c, c′: T(c)T(c′) = T(cc′)·(c/c′)₂·s(c,c′). Hence
   R(m)² = Σ_n T(n) χ̄_n(m) τ̃(n),   τ̃(n) = Σ_{cc′=n} (c/c′)₂ s(c,c′) W(N(c)/Z) W(N(c′)/Z).
For squarefree n = π₁⋯π_k, τ̃(n) is a sum over subsets S of Π_{i∈S, j∉S}(π_i/π_j)₂: the Rédei-matrix sums that govern 4-ranks of
class groups (Fouvry–Klüners). So the 4th moment Σ_m|R(m)|⁴ is a (G)-type mean square with Rédei-weighted coefficients of length Z².
This is a structural link to Fouvry–Klüners-type oscillation of Jacobi symbols in many variables, but it addresses the 4th moment,
not the 2nd, and gives no handle on (L³).

### S25. Where the "degree 3" comes from (HEURISTIC, from memory of Kazhdan–Patterson)
Kazhdan–Patterson: the theta representation of the n-fold cover of GL_r is distinguished (unique Whittaker model) iff n = r or n = r+1.
For n = 4 this allows GL4 and GL3. The GL3 case (n = r+1) is the natural home for quartic Gauss sums as Whittaker coefficients. This
would explain subagent A's degree-3 count (twists by (m/·)_4 have conductor N(m)³): the Kubota reflection is then a GL3-type Voronoi
formula. Under this picture (L³) is literally a second moment of L(½, Θ₃ ⊗ (m/·)_4) over Q twists of a degree-3 metaplectic theta form,
the metaplectic twin of the open GL3 problem in S6. A delta-method attack (Munshi-style, which gives single-twist subconvexity for genuine
GL3) would need Voronoi for Θ₃ with additive twists. That is plausible from the GL3 automorphy, but a family-average saving beyond
Q^{3/2} is not known even for genuine GL3. Citation from memory; not verified this session.
- (S23 addendum) corr(|R_tw(m)|², |R_untw(m)|²) over rows = +0.0007 (Z = 1e6) and −0.0073 (Z = 3e6), with s.e. 0.02. The α-twisted and
  untwisted rows of the same m are uncorrelated: no common "arithmetic size factor" depending on m.

### S26. Parity-blind sums and quartic-theta terms at small indices (l3/gen_T2.c dump to N(c) = 1e6, l3/theta_fit.py) — NUMERICAL
|K_1(m,Z)|/Z^{3/4} (K_1 = Σ_c g̃(χ̄_c,m)W, index m) on 13 log-spaced Z in [1e4, 1e6]; noise level ≈ 0.155·Z^{−1/4} (0.016 → 0.005):
- index 1: ≈ 0.010 (0.013, 0.015, …, 0.008, 0.009), roughly constant, ~2× noise at the top: a theta term is plausibly present but small.
- index π5² (a square): ≈ 0.018 (0.019, 0.026, 0.024, 0.025, 0.022, …), constant: clear Z^{3/4} theta term.
- squarefree non-square indices π5, π13, π5·π13, and π5³: decay toward the noise floor (π5: 0.016 → 0.001), so the theta
  coefficients there are small or zero at this resolution. [The "concentrated on squares" reading is WRONG; see S29.]
Möbius rows at small index (S13 follow-up, 121 Z-points, decade bins), mean |K|²/Z^{3/4}:
  index 1: 0.28, 0.42, 0.16, 0.48 (stable, fits Z^{3/8}); π5: 0.16, 0.36, 0.62, 0.62; π13: 0.22, 0.44, 0.44, 0.57 (closer to √Z).
  INCONCLUSIVE between the Z^{3/8} zero-term picture and √Z for π5, π13. S13's heuristic is supported at index 1 only.
- (S23 addendum) Row phases arg R(m) at Z = 3e5 (1759 rows): 8-bin histogram [239 212 205 225 231 205 212 230] (expected 220);
  |mean e^{2i·arg}| = 0.020, |mean e^{4i·arg}| = 0.021 (random 0.024). Odd moments vanish automatically, since R(m̄) = conj R(m).
  Phases are uniform, so rows look like complex Gaussians in modulus and phase.

### S27. Theorem Z re-derived (outline)
PO(η,λ): Σ_{N(k)≤H} |Σ_{c~D} μ(c)χ_k(c)|² ≤ (H + D)D^{1+λ+ε}, H = D^η. Take k = k₀j⁴, N(j) ≤ (H/N(k₀))^{1/4}: χ_k = χ_{k₀} on c coprime to j.
Positivity gives |S_{k₀}(D)|² ≲ (H + D)D^{1+λ}/H^{1/4}, so |S_{k₀}(D)| ≲ D^{λ/2 + max(1 − η/8, ½ + 3η/8)}. Hence L(s, χ_{k₀}) ≠ 0 beyond that
exponent. Caveat: for the finite sums, the coprimality to j changes S by shorter sums, which need an inductive treatment (outline only). The optimum is η = 1, giving 7/8 + λ/2. Matches NOTES §8 exactly.

### S28. The Z = 3e6 maximum (NUMERICAL)
Row max |R|²/diag = 13.14 at N(m) = 5777 = 53·109. Under Exp(1) with 4086 independent pairs, P(max ≥ 13.14) = 0.008, while maxima at
Z = 1e5, 3e5, 1e6 had p = 0.53, 0.44, 0.64. Two of the top three norms contain 53, but rows with 53 | N(m) are not inflated on average.
Per-prime means over rows with p | N(m): the largest is z = 3.15 among ~100 primes with ≥ 20 rows, which is the expected maximum of
~100 normal draws. No arithmetic structure; the single large row is a look-elsewhere event. The Z = 6e6 run checks the tail again.

### S29. Theta terms at prime indices follow N(π) mod 8 (l3/theta_fit2.py, theta_fit3.py; Z ∈ [1e5, 1e6], 9 points) — NUMERICAL
Mean of K_1(π, Z)/Z^{3/4} (complex; noise ≈ ±0.002 on the mean):
  N(π) = 17 (≡ 1 mod 8): 0.0244 ± 0.0104i (|·| ≈ 0.027)   N(π) = 41 (≡ 1 mod 8): 0.0145 ∓ 0.0205i (|·| ≈ 0.025)
  N(π) = 5, 13, 29, 37, 53 (≡ 5 mod 8): all |mean| ≤ 0.002, i.e. zero within noise.
  Conjugate primes give complex-conjugate coefficients (as G(π̄) = conj G(π) predicts).
  Other indices: 1 → ≈ 0.010; π5² → 0.015; π13² → 0.017; 5 = π5π̄5 → 0.013; π17² → 0.007; 3, 3², 7² (inert) → ≤ 0.011, unclear.
Observed rule: the quartic-theta (s = 3/4) residue at a prime index π is nonzero iff N(π) ≡ 1 (mod 8), i.e. iff χ_π(−1) = +1. That is a
local vanishing like the classical theta vanishing off squares. It corrects S26 (π17 is squarefree, non-square, and carries the largest
term), and it explains the untwisted bias row N(m) = 17³ in S1/S11 (index m³ = π17⁹ ≡ π17 modulo fourth powers). The α-twisted
family used for (L³) has no s = 3/4 pole at all, so none of this enters (G).

### S30. Sharp cutoff W = 1 on (Z/2, Z] (l3/gen_rows_sharp.c, l3/sharp.txt) — NUMERICAL
Z = 3e5: S_tw/diag = 1.010, S_untw/diag = 1.129. Z = 1e6: S_tw/diag = 0.984, S_untw/diag = 1.128.
(G) is insensitive to the weight. The untwisted control is relatively more biased with a sharp cutoff: the theta residue is weighted by
W̃(3/4), which is larger relative to ‖W‖₂ for the indicator.
- (S29 addendum, l3/theta_cls.py) Split by the class of N(c) mod 8: for π5 and π13 (N ≡ 5 mod 8) BOTH c-classes give ≈ 0 (≤ 0.0025),
  so the vanishing is not a cancellation between classes. For π17 the classes give 0.0072+0.0023i and 0.0172+0.0081i; for π41,
  0.004−0.0069i and 0.0105−0.0136i; for index 1, +0.0192 and −0.0081. The residue vanishes classwise at N(π) ≡ 5 mod 8.
- (S29 literature) A web search found no statement of this vanishing rule. The natural reference is Eckhardt–Patterson, "On the Fourier
  coefficients of biquadratic theta series" (Proc. LMS 1992), which was unreachable here (arXiv mirrors blocked), as was DDHL
  (arXiv 2306.11875), which uses Suzuki's evaluation at squares. The rule may be known there; it is recorded only as an observation.

### S29 confirmation (dump with test primes up to norm 300, l3/T2_q300.bin, l3/theta_rule.py) — NUMERICAL
All 61 split primes with N(π) ≤ 300: for N ≡ 5 mod 8 (34 primes) |mean K_1/Z^{3/4}| ≤ 0.0042 (median 0.0015); for N ≡ 1 mod 8
(27 primes) ≥ 0.0066 (median 0.0102). There is no overlap, so the rule "the theta residue at π vanishes iff N(π) ≡ 5 mod 8" holds on all of them.
Inert prime squares 3², 7², 11² (N ≡ 1 mod 8) give 0.0085, 0.0066, 0.0114, also nonzero.
Finer structure within N ≡ 1 mod 8: the residue is LARGE (0.019–0.027) for p = 17, 41, 97, 137, 193, 241 and SMALL (0.008–0.010) for
p = 73, 89, 113, 233, 257, 281. The small group is exactly the p with p = x² + 64y², i.e. where 2 is a biquadratic residue (Gauss):
the size of the quartic-theta coefficient at π depends on the quartic character of 2 at π. NUMERICAL observations; no proof attempted.

### S31. Windows up to Z = 2.5e6 (l3/multi25e5.txt) — NUMERICAL
Z = 2.5e6: 0.990; 1.25e6: 1.007; 6.25e5: 1.044; 3.1e5: 1.031; 1.6e5: 1.060; 7.8e4: 0.994. Combined with S15 (36 windows):
weighted mean S_tw/diag = 1.0081 ± 0.0068 (s.e. underestimated by up to √3 because windows overlap). Consistent with exactly 1.
- (S29 phases, l3/theta_phase*.py) For the 24 primes with N ≡ 1 mod 8, the phase of the theta term shows no simple relation to arg π, arg π/2,
  or the Gauss sum g̃(χ̄_π) (tested v/g, v·g, v/conj g, v·α^{1/2}/g). Conjugate primes give conjugate values, as they must. Unlike the cubic
  case (Patterson: τ(π) ∝ conj g(π)), the quartic theta coefficient at π is not a simple function of g(π). That is consistent with the
  4-fold GL2 theta not being distinguished (Kazhdan–Patterson).
- (S29 explained as a 2-adic condition, l3/theta_units.py) Residues at index i^u·h, u = 0..3:
    h = 1: 0.0111, 0.0121−0.0118i, 0.0273, 0.0121+0.0118i;   π17: all four nonzero;
    π5 = −1−2i: u=0 → 0.0004−0.0005i ≈ 0, u=2 (−π5 = 1+2i) → ≈ 0, but u=1 (2−i) → 0.015+0.0035i and u=3 → 0.0043−0.0136i (nonzero).
  For primary π, N(π) ≡ 5 mod 8 ⟺ π ≡ 3+2i mod 4. So the observed rule is a 2-adic congruence on the INDEX: the residue vanishes for odd
  h ≡ 1+2i or 3+2i mod 4 and not for h ≡ ±1 or ±i-type classes. That is a local condition at 2 of the kind theta coefficients
  usually satisfy, not a deep global phenomenon. The x² + 64y² size split within N ≡ 1 mod 8 is also 2-adic in nature (quartic character of 2).
- (S6 check) #E_D(F_p) = p + 1 − 2·Re(conj((D/π)_4)·π) for E_D : y² = x³ − Dx and primary π | p: 125/125 cases (14 primes p ≤ 113,
  9 values of D). The alternative normalizations fail (69/125, 56/125, 0/125). So α(c)χ̄_c(m) = ψ_{E_m}(c)/|c| with the stated convention.
- (Identity 1 at scale) T(c)²/α(c) is real ±1 to 5e-12 for ALL 347 518 squarefree primary c with N(c) ≤ 1e6. Its sign equals μ(c) on
  20 000 random c (ω(c) computed from N(c): split p count v_p, inert q count v_q/2), with 0 mismatches. Note this checks the
  twisted-multiplicativity pipeline against μ, not an independent Gauss-sum evaluation (that is S17, N ≤ 600).

### S32. Prime-only columns break (G) for the α-twisted family (l3/gen_rows_prime.c) — NUMERICAL
R_P(m) = Σ_{π prime, N(π)~Z} T(π)χ̄_π(m)W at Q = Z^{2/3}:
   Z = 1e5: S_tw/diag = 1.228, S_untw/diag = 1.004;   Z = 3e5: S_tw/diag = 1.286, S_untw/diag = 0.979.
Mechanism (PROVED identity S2): at primes T(π) = −conj g̃(χ̄_π), so R_P(m) = −conj Σ_π g̃(χ̄_π)χ_π(m)W. That is a parity-blind prime sum at
INDEX m, which carries the quartic-theta term (Patterson-type, ~Z^{3/4}/log Z, with a 2-adic support rule, S29). The excess grows like
Z^{1/2}/log Z (predicted 0.23 → 0.36, observed 0.23 → 0.29).
Consequence for the quartic route (S8): the α-twist removes the theta pole only for Möbius-structured squarefree sums. Any quartic
analogue of the prime-slot (marked) moments carries quartic-theta main terms, like the exceptional Gram terms of the sextic proof, and
needs the theta coefficients τ(m), which are not explicitly known at non-square indices (S29 phases).
- (S32 split, Z = 3e5) Prime-column twisted rows: N(m) ≡ 1 mod 8 → mean |R_P|²/diag = 1.570 (875 rows); N(m) ≡ 5 mod 8 → 1.005
  (884 rows). The top 10 rows all have N(m) ≡ 1 mod 8. The bias is exactly the theta support rule of S29 (index m ≡ 3+2i mod 4 has zero
  residue), which confirms the mechanism.
- (contrast) In the full squarefree (G) rows there is no such split: at Z = 3e6, N(m) ≡ 1 mod 8 rows average 1.037 (classes 1, 9 mod 16)
  and N(m) ≡ 5 mod 8 rows average 1.047 (classes 5, 13 mod 16) (S11). The Möbius structure removes the theta term, as S9/S10 predict.

### S33. Theta residues at cube indices vanish (l3/theta_cube.py) — NUMERICAL
|mean K_1/Z^{3/4}|: π17: 0.0265, π17²: 0.0056, π17³: 0.0015; π41: 0.0251, π41³: 0.0010; π5³: 0.0007; π17·π41: 0.0094; (π17π41)³: 0.0015
(noise ≈ 0.002). The residue at exponent 3 is zero within noise even for N ≡ 1 mod 8.
This explains S32's asymmetry. Untwisted prime rows correspond to index m³ (zero residue): unbiased. α-twisted prime rows correspond to
index m (nonzero for m ≡ 1 mod 4): biased. It also fits the untwisted (G) bias rows N(m) = 17³: index π17⁹ = π17·(π17²)⁴ has exponent ≡ 1 mod 4.
Observed local pattern for the parity-blind residue at π^e (primary π): e = 1 nonzero iff π ≡ 1 mod 4; e = 2 small but nonzero; e = 3 zero.
- (S33 addendum, l3/theta_sq.py) Index π² over all primes ≤ 300: N ≡ 5 mod 8 → nonzero (median 0.0087, range 0.007–0.016, decaying in p);
  N ≡ 1 mod 8 → small (median 0.0039, max 0.0056, near the noise). So the pattern REVERSES between e = 1 and e = 2:
  residue at π^e nonzero for (e=1, N≡1 mod 8), (e=2, N≡5 mod 8), never for e = 3 (S33). This is a local (Hecke-type) structure at π
  depending on χ_π(−1) and e. It is recorded as a curiosity about the parity-blind quartic theta and is not used for (L³).
- (S32 at Z = 1e6) Prime columns: S_tw/diag = 1.442, S_untw/diag = 0.935. Excess 0.228 → 0.286 → 0.442 for Z = 1e5, 3e5, 1e6;
  Z^{1/2}/log Z scaling from 3e5 predicts 0.48 at 1e6. Confirmed: the prime-column α-twisted mean square grows like Z^{1/2}/log Z
  relative to the diagonal, so it is NOT ≪ Q·Z^{1+ε}/log.
- (S32 split across Z) prime columns, mean |R_P|²/diag: N(m) ≡ 1 mod 8 → 1.408, 1.570, 1.911; N(m) ≡ 5 mod 8 → 1.046, 1.005, 0.975
  (Z = 1e5, 3e5, 1e6). The growth is entirely in the theta-supported class.
- (S29 final, all split primes with N(π) ≤ 1000, l3/T2_q1000.bin, l3/rule_q1000.txt) |residue/Z^{3/4}| in three disjoint groups:
    N ≡ 5 mod 8 (86 primes): ≤ 0.0042 (zero within noise);
    N ≡ 1 mod 8, 2 a biquadratic residue (p = x² + 64y², 28 primes): 0.0061 – 0.0102 (median 0.0082);
    N ≡ 1 mod 8, 2 not a biquadratic residue (46 primes): 0.0146 – 0.0265 (median 0.0191).
  The ranges do not overlap. Robust numerical structure of the parity-blind quartic theta residue at prime index.
- (S30 addendum) A wider smooth weight (bump on (1/4, 1]) at Z = 1e6, Q = Z^{2/3}: S_tw/diag = 0.984, S_untw/diag = 1.035.

### Next steps if this line is resumed
1. Literature: get Eckhardt–Patterson (Proc. LMS 1992) and Suzuki on biquadratic theta coefficients, to compare with the S29/S33
   pattern (vanishing at π ≡ 3+2i mod 4 for e = 1, the reversal at e = 2, zero at e = 3, the x² + 64y² size split). If it matches a known
   local formula, the parity-blind main terms are explicit, which is what a quartic prime-slot moment (S32) would need.
2. Try to prove the square-case optimal large sieve for T-coefficients (sub-problem 1 in L3_REPORT §12). Numerically 1.006 at Q = Z = 3e5.
   With Theorem Z it would give a 7/8 zero-free region for all quartic Hecke L-functions over Q(i), a self-contained new result.
3. A ratios-type theorem for the quartic family in a restricted range (Čech-style MDS methods for quadratic characters are the model),
   aiming at the S20 trade-off curve.
4. For θ in the proven sextic architecture: a supercritical marked inverse moment (MM′). Every unit of capacity τ is worth ≈ 0.043τ in θ.
- (S29 symbols) Over the 74 split primes N ≡ 1 mod 8, N ≤ 1000, grouped by ((1+i)/π)_4 and (i/π)_4: the residue is large iff
  ((1+i)/π)_4²·(i/π)_4 = (−2/π)_4 = (2/π)_4 = −1. This is the same statement as the x² + 64y² criterion, now as a quartic symbol of 2 at π.

### S34. Size of the implied Möbius off-diagonal versus the (R) error H^{1/2} — NUMERICAL
Z=1e+05: Off_S/diag = +0.0055;  |Off_S|/(Q*H^{1/2}*diag/(QZ)) = 0.254  (|.| <~ 1 is consistent with (R))
Z=3e+05: Off_S/diag = +0.0226;  |Off_S|/(Q*H^{1/2}*diag/(QZ)) = 1.516  (|.| <~ 1 is consistent with (R))
Z=1e+06: Off_S/diag = -0.0042;  |Off_S|/(Q*H^{1/2}*diag/(QZ)) = 0.417  (|.| <~ 1 is consistent with (R))
Z=3e+06: Off_S/diag = +0.0421;  |Off_S|/(Q*H^{1/2}*diag/(QZ)) = 6.068  (|.| <~ 1 is consistent with (R))
(Interpretation: Q·H^{1/2}·(diag/QZ) = diag/√Q, the same order as the random-model fluctuation of S. In units of the actual standard
error √(var/(NR/2))·diag the four values are 0.1σ, 0.6σ, 0.2σ, 2.6σ. So the implied off-diagonal is at the (R)/random-model scale
H^{1/2}, with the 3e6 value the largest. Z = 6e6 decides whether that is a fluctuation.)

### S35. No theta term in the Möbius-weighted (L³) family (l3/mob_rule.py, all 160 split primes N(π) ≤ 1000) — NUMERICAL
|mean K_{−1}(π,Z)/Z^{3/4}| over Z ∈ [1e5, 1e6]: N ≡ 1 mod 8 → median 0.0018, max 0.0040; N ≡ 5 mod 8 → median 0.0020, max 0.0046.
All at noise level, including the primes where the parity-blind sum has residues 0.006–0.027 (S29). Direct confirmation that the
Möbius weight (equivalently the α-twist on the index-m³ side) kills the s = 3/4 pole at every prime index tested.
- (S35 composite indices, l3/mob_rule2.py) 60 random products of two primes (each N ≤ 200): Möbius-weighted median 0.0021, max 0.0045
  (noise); parity-blind median 0.0020, max 0.0265 (theta terms at some composites). Same conclusion.
- (S33 table, split primes N ≤ 400, l3/theta_e3.py) |residue|: e=2, N≡1 mod 8: median 0.0039 (max 0.0056); e=2, N≡5 mod 8: median 0.0083
  (min 0.0062); e=3: ≤ 0.0038 for both classes (noise). The e=2 groups separate (max 0.0056 < min 0.0062), and e=3 is zero.

### Retrospective (written 11:03Z)
- The decisive structural facts are the exact identities T² = μα and T = μ·conj g̃. They make (G) both self-dual (Poisson and reflection
  fix M = 2/3) and parity-sensitive (Möbius-weighted Kubota sums). Every unconditional tool tried either preserves the fixed point or
  is parity-blind.
- What can be proved now: the conditional theorem from a ratios-type hypothesis; the reformulations; the numerical picture. What
  would be new and maybe within reach: the square-case optimal large sieve for T-coefficients, which via Theorem Z gives a 7/8 zero-free
  region for the quartic Hecke family. That is not a θ improvement for Dirichlet L-functions, though.
- For θ itself (the tracker), every remaining lever in the proven sextic architecture is either negligible (B, slot cap: ≤ 4e-6) or
  sits at a theta/self-dual threshold (marked inverse capacity). So the next real gain needs a new idea of the same difficulty class.
- (S17 extended) T(c) from the pipeline vs brute force, matched element-by-element (c itself, from the gen_T2 dump), for all 416
  squarefree primary c with N(c) ≤ 1200: max error 4.9e-13.
- (S17 extended) plus all 624 squarefree primary c with 1200 < N(c) ≤ 3000: max error 9.7e-13. In total 1040 c checked element-wise.
- (S17 extended) quartic-symbol exponents: 31 720 checks (all squarefree c with N ≤ 1500 × 61 test primes up to norm 300), 0 mismatches.
- (identity 2 at scale) T(ab) = T(a)T(b)(a/b)_2 s(a,b) on 1344 random coprime prime pairs (products up to N ≤ 1e6): 0 mismatches. Partly
  circular, since the pipeline builds composites with the equivalent rule G(c′π) = s(c′/π)_2 G(c′)G(π). The independent check is the
  element-wise brute force of S17 for all composites with N ≤ 3000.
- (S17 extended) plus 1050 c with 3000 < N(c) ≤ 6000: max error 9.1e-13. Total 2090 squarefree primary c with N(c) ≤ 6000 checked element-wise.
- (S17 extended) plus 1393 c with 6000 < N(c) ≤ 10000: max error 7.1e-13. Total 3483 squarefree primary c with N(c) ≤ 10 000 checked.
- (S17 extended) quartic-symbol exponents: 288 010 checks (all squarefree c with N ≤ 5000 × 166 test primes up to norm 1000), 0 mismatches.
- (S17 extended) plus 854 c with 10000 < N(c) ≤ 12500: max error 1.0e-12. Total 4337 squarefree primary c with N(c) ≤ 12 500 checked element-wise.
- (S29 inert primes) index = the primary associate −q of an inert prime q (≡ 1 mod 4): |residue| 0.0085, 0.0066, 0.0114, 0.0078, 0.0071,
  0.0056 for q = 3, 7, 11, 19, 23, 31. All nonzero, consistent with the 2-adic rule (index ≡ 1 mod 4 ⇒ nonzero).
- (S17 extended) plus 693 c with 12500 < N(c) ≤ 14500: max error 1.3e-12. Total 5030 squarefree primary c with N(c) ≤ 14 500 checked.

### S36. Z = 6e6 (l3/scale6e6.txt, rows_Z6000000.bin; 6757 s) — NUMERICAL
Q = 33019, 12 961 rows: S_tw/diag = 1.0101 ± 0.0128, S_untw/diag = 1.0445. Var 1.058, E[x²]/E[x]² = 2.04, E[x³]/E[x]³ = 6.21 (Gaussian 2, 6).
P(>3) = 0.0556 (Exp 0.0498), P(>5) = 0.0069 (Exp 0.0067). Max 9.09 at N = 28289, P(max ≥ obs) = 0.52 (typical).
N(m) ≡ 1 mod 8: 1.005; ≡ 5 mod 8: 1.013 (no theta split). The 3e6 excess and max (S11, S28) were fluctuations.
- (S17 extended) plus 353 c with 14500 < N(c) ≤ 15500: max error 9.7e-13. Total 5383 squarefree primary c with N(c) ≤ 15 500 checked.
- (S32 at Z = 3e6) Prime columns: S_tw/diag = 1.672, S_untw/diag = 1.002. Excess 0.228, 0.286, 0.442, 0.672 for Z = 1e5, 3e5, 1e6, 3e6.
  Z^{1/2}/log Z scaling from 1e6 predicts 0.71 at 3e6, observed 0.67. The prime-column bias grows exactly like the theta term predicts.

### S37. Reflection pair at Z = 1e6 (l3/mscan1e6.txt) — NUMERICAL
M = 0.6 (Q = 3981, 1566 rows): S_tw/diag = 1.012;  M = 0.75 (Q = 31623, 12 419 rows): S_tw/diag = 1.004. Both ≈ 1, as at Z = 3e5 (S14).
M = 0.85 at Z = 1e6 (Q = 125 893, 49 439 rows): S_tw/diag = 0.998 (s.e. ≈ 0.0065), S_untw/diag = 1.005.

## Round 2 (resumed 2026-10-15, after "then keep finding L3 proves")

### S38. GRH ⟹ optimal large-sieve shape for the T-family — PROVED (conditional on GRH)
**Theorem G1.** Assume GRH for ζ_{Q(i)} and for the Hecke L-functions L(s, ξ) of the finite-order characters
ξ(c) = (r/c)₄·ν(c) of Q(i) (r ∈ Z[i]∖0, ν a class character mod 16). Let w, W be smooth, supported in [½, 2] and [½, 1].
Then for all Q, Z ≥ 2,
  S_w(Q,Z) = Σ_{m primary} w(N(m)/Q) |Σ_{c sqfree primary} T(c) χ̄_c(m) W(N(c)/Z)|²  ≪_ε (QZ)^ε (QZ + Z²).
Consequences:
* (G) holds for every M ≥ 1. In particular sub-problem 1 of L3_REPORT §12 (the square case Q = Z) holds on GRH.
* At M = 2/3 this gives S ≪ Z^{2+ε} = Q^{1/2}·QZ^{1+ε}. So on GRH, (L³) holds up to a factor Q^{1/2}, i.e.
  Σ_{N(m)≤Q} |ψ_α(m³, ½+it)|² ≪ Q^{3/2+ε}(1+|t|)^A, via the standard approximate functional equation of ψ_α (lengths Q^{3/2}).
  Compare: unconditional Q^{15/8} (Gao–Zhao) and Q^{7/4} (recursion); target Q^{1+ε}.
**Proof.**
1. Open the square: S = Σ_{c,c'} T(c) conj T(c') W W' Σ_{m ≡ 1 (λ³)} w(N(m)/Q) ψ(m), with ψ = χ̄_c χ_{c'} periodic mod L = lcm(c,c').
   Write g = gcd(c,c'), c = g c₁, c' = g c₁'. Then ψ = χ̄_{c₁} χ_{c₁'} · 1_{(·,g)=1}, and L = g c₁ c₁' is squarefree.
2. Poisson over the coset 1 + λ³Z[i] (Z[i] is self-dual for (x,y) ↦ e(Im xy)), and CRT across λ³ and the three coprime moduli g, c₁, c₁':
   Σ_m w(N(m)/Q)ψ(m) = (Q/(8N(L))) Σ_{k ∈ Z[i]} φ₀(k)·(unit/class factors)·r_g(k)·g(χ̄_{c₁}, k)·g(χ_{c₁'}, k)·Ŵ_w(Q N(k)/(8N(L))),
   where r_g is the Ramanujan sum mod g, |φ₀| = 1 depends on k mod λ³, the class factors depend on c₁, c₁', g mod 16, and Ŵ_w(N(ξ)) is the
   2-D Fourier transform of w(N(·)), a Schwartz function of N(ξ) ∈ [0,∞).
3. Collapse the Gauss sums (identity checked in S19 for all k, coprime or not): T(c₁) g(χ̄_{c₁}, k) = μ(c₁) χ_{c₁}(k) √N(c₁), and
   conj T(c₁') g(χ_{c₁'}, k) = μ(c₁') χ_{c₁'}(−1) χ̄_{c₁'}(k) √N(c₁'). Also T(g c₁) = T(g)T(c₁)(g/c₁)₂ s(g,c₁) (twisted multiplicativity), and |T(g)| = 1.
   After quartic reciprocity, c₁ ↦ χ_{c₁}(k)(g/c₁)₂ = (k g²/c₁)₄^{±1} times a class factor, i.e. a Hecke character ξ of c₁ of the type in the hypothesis.
   Splitting the reciprocity signs into products (as in L3_REPORT §10, step 3), the (g, k) term of S becomes a bounded combination of
     (Q/(8N(g))) r_g(k) Σ_{c₁, c₁'} μ(c₁)μ(c₁') ξ(c₁) ξ'(c₁') N(c₁c₁')^{−1/2} W(N(gc₁)/Z) W(N(gc₁')/Z) Ŵ_w(Q N(k)/(8N(g)N(c₁c₁'))),
   with c₁, c₁' coprime to g and to each other. (Coprimality between c₁ and c₁' is removed by Möbius inversion over their gcd at the cost of
   another divisor sum of the same shape.)
4. Separate variables: Ŵ_w(x) = (2πi)^{−1} ∫_{(1+ε)} M(z) x^{−z} dz, with M(z) decaying faster than any power of |Im z|. The (g, k) term becomes
   ∫ M(z) (Q/(8N(g)))^{−z} N(k)^{−z} P_z(k) P'_z(k) dz, where P_z(k) = Σ_{c₁, (c₁,g)=1} μ(c₁) ξ(c₁) N(c₁)^{−1/2+z} W(N(gc₁)/Z).
5. GRH: Σ_{N(c)≤y, (c,g)=1} μ(c) ξ(c) ≪ y^{1/2} (y N(k g))^ε uniformly (Perron plus the GRH bound for 1/L on Re s = ½ + 1/log log, with the finitely
   many Euler factors at primes of g absorbed into (N g)^ε). Partial summation gives |P_z(k)| ≪ (Z/N(g))^{1+ε} (Z N(k))^ε (1+|z|)².
6. Sum over k with Re z = 1 + ε: Σ_k N(k)^{−1−ε} |r_g(k)| ≤ Σ_{d|g} N(d) Σ_{d|k} N(k)^{−1−ε} ≪ τ(g)/ε. So the g-part of S is
     ≪ (Q/N(g)) · (N(g)/Q) · (Z/N(g))² · τ(g) · Z^{O(ε)} = τ(g) Z^{2+O(ε)} / N(g).
   Summing over squarefree g with N(g) ≤ Z gives ≪ Z^{2+O(ε)}. The k = 0 term occurs only for c = c' (ψ principal) and is the diagonal
   Σ_c |T(c)|² W² Σ_m w·1_{(m,c)=1} ≪ QZ. Hence S ≪ (QZ + Z²) Z^{O(ε)}. ∎
**Remarks.**
* This is the quartic T-analogue of the GRH dispersion estimates of Dunn–Radziwiłł (cubic case, Annals 2024, Prop. 9.1–9.2), whose error terms also
  have the large-sieve shape B^{2+2η}. Here the result is cleaner: T = μ·conj g̃ carries a Möbius factor, so the Poisson dual has no theta main term
  (in DR the dual of g̃(b)β_b produces the cubic-theta correction term).
* It shows that the obstruction to the optimal quartic large sieve found for general coefficients (dFDH, arXiv 2607.07911, and DR Theorem 1.4 in
  the cubic case: the extremal sequences are Gauss sums g̃(b) themselves) disappears for T-coefficients, conditionally on GRH.
* It does NOT give (L³). GRH bounds each |P(k)| by Z^ε, so the off-diagonal is ≪ H^{1+ε}, H = Z²/Q. (L³) needs Σ_k (|P(k)|² − C) ≪ Z = H^{3/4}:
  a power saving in the family variance, which GRH does not control (S12). The remaining gap is exactly Q^{1/2} in the (L³) normalization.

### S39. The multiple-Dirichlet-series route for the second moment has an infinite (hyperbolic) group — PROVED (linear algebra; l3b/mds_group.py)
Three-variable series Z_T(s₁,s₂,w) = Σ_m N(m)^{−w} ψ_α(m³,s₁) conj ψ_α(m³, s̄₂). In centred coordinates x_j = s_j − ½, y = w − ½ the functional equations act by
* K₁, K₂: Kubota FE in s₁ or s₂ (conductor N(m)³): x_j ↦ −x_j, y ↦ y + 3x_j;
* P: Poisson in m, T-form ↔ Möbius form Σ_k N(k)^{−½−b}/(L(½+a₁,ξ_k)L(½+a₂,ξ̄_k)): (x₁,x₂,y) ↦ (x₁+y, x₂+y, −y);
* J: the joint FE of the two 1/L factors (root numbers cancel to a class sign): (a₁,a₂,b) ↦ (−a₂, −a₁, b − a₁ − a₂).
Conjugating J into T-coordinates gives J_T = [[−1,−2,−4],[−2,−1,−4],[1,1,3]]. K₁, K₂, J_T are reflections; their unique invariant symmetric form has
eigenvalues −5/2, (3 ± √19)/2, so it is indefinite (signature (1,2)). K₁J_T has eigenvalue −9.899 (hyperbolic, infinite order), and the number of group
elements of word length ≤ n grows like 1.63ⁿ (4, 9, 17, 30, 51, 85, 140, 229). Single FEs of one 1/L factor add Gauss-sum-weighted forms and more
generators, so the full group is also infinite. Consequence: Bochner-type continuation from a finite group (the Diaconu–Goldfeld–Hoffstein mechanism)
is unavailable. This is the same phenomenon as for the sixth moment of quadratic Dirichlet L-functions (hyperbolic Weyl group), which is consistent
with the GL3-sixth-moment analogy of S12.
Note: this re-derives, with a different generating set, the earlier finding of subagent F (NOTES.md §8: "second-moment MDS group is hyperbolic
Coxeter (2,6,6) ... Under GRH these methods give exactly Q^{3/2}"). Here the subgroup ⟨K₁, K₂, J_T⟩ is of type (2, ∞, ∞). Both are hyperbolic,
with the same conclusion. S38 is the written-out proof of F's GRH remark.

### S40. Literature located this round (arXiv PDFs read in ../lit/)
* David–Dunn–Hamieh–Lin, "Quartic Gauss sums over primes and metaplectic theta functions", arXiv 2306.11875 (v5, Jan 2026).
  - Theorem 1.1: Σ g̃₄(c)(c/|c|)^ℓ Λ(c) R(N(c)/X) ≪ X^{5/6+ε} + X^{3/4+ε}|ℓ|^{3/2+ε}. Conjecture 1.2: bias b_β X^{3/4} for ℓ = 0, X^{1/2+ε} for ℓ ≠ 0.
    (Our S32 prime-column growth ∝ Z^{1/2}/log Z is consistent with this X^{3/4} bias at prime moduli; our S26/S29 prime-sum bias, +1838 in class
    1 mod 4 and −760 in class 3+2i, matches the shape of Conjecture 1.2 with two classes β of the modulus.)
  - Type II via g̃₄(ab) = g̃₄(a)g̃₄(b)(a/b)₂ and the quadratic large sieve over Q(i) (Onodera; Goldmakher–Louvel). Same cocycle as our identity 2.
  - Suzuki: theta coefficient ψ^{(4)}(a²) = g̃₄(a)/N(a)^{1/4} for squarefree primary a; other theta coefficients undetermined.
  - (1.11): Σ_{N(a)≤A} μ²(a)|ψ̃^{(4)}(½+ε+it, a², 0)|² ≪ A^{3/2+ε}(|t|+1)^{3+ε}; Remark 9.7 says Lindelöf on average A^{1+ε} "seems likely with
    additional work". That is the index-a² analogue of (L³): quadratic characters (optimal large sieve available) and conductor N(a)² (ratio 2).
    (L³) is the index-m³, ℓ = 1 case: quartic characters and conductor N(m)³ (ratio 3), so strictly harder.
  - In their notation ψ_α(m³, s) = ψ̃^{(4)}(s, m³, 1) up to normalization; ℓ = 1 has no pole, which matches our S35 (no theta term in the T-family).
* Dunn–Radziwiłł, "Bias in cubic Gauss sums: Patterson's conjecture", Annals 200 (2024), arXiv 2109.07463. Under GRH, Heath-Brown's cubic large
  sieve is sharp (extremal sequence β_b = g̃(b)); a corrected large sieve (their (1.13), Prop. 9.1–9.2) for sequences β_b g̃(b) with β cancelling
  against cubic characters. S38 is the quartic T-family analogue; the Möbius factor removes the correction term.
* de Faveri–Dunn–Hoffstein, arXiv 2607.07911 (July 2026): the cubic and quartic large sieves are not perfectly orthogonal (unconditional), and the
  obstruction is Gauss-sum bias.
* Gao–Zhao, arXiv 1101.3646: quartic large sieve with (MN)^{3/4}, both variables squarefree.

### S41. Euler-product control: (L³) rows carry no m-dependent size (l3/gen_corr.c, l3/corrstats.py, l3/corr1e6.txt) — NUMERICAL
Fix Q = 9999 and compute the rows for four disjoint dyadic c-windows Z ∈ {1e6, 5e5, 2.5e5, 1.25e5}, on the 3482 squarefree m ≠ 1 with N(m) ≤ Q.
Two families are compared:
* the T-family R(m);
* the control Σ_{c primary} d₃(c) χ̄_c(m) W(N(c)/Z), the partial sums of L(s, ξ_m)³. This is an Euler product with the same conductor N(m)³.
The statistic is x = |R|²/rowdiag.

| statistic | T-family | d₃ control (Euler product) |
|---|---|---|
| S/diag per window | 1.007, 1.034, 1.011, 1.047 | 0.80, 0.88, 0.92, 0.92 |
| E[x²]/E[x]² (Gaussian: 2) | 2.00, 2.04, 2.00, 1.93 | 4.09, 4.04, 3.62, 3.71 |
| E[x³]/E[x]³ (Gaussian: 6) | 5.97, 6.47, 6.02, 5.28 | 34.8, 33.8, 25.1, 28.9 |
| corr(x_j, x_k) over the 6 window pairs | −0.040 … +0.020 (mean −0.008; s.e. 0.017) | 0.40 … 0.52 (≈ 28 s.e.) |

At Z = 1e5 (748 rows) the picture is the same: T-family corr ≤ 0.06 (s.e. 0.037) with Gaussian moments; d₃ corr 0.31–0.48 and E[x²]/E[x]² ≈ 3.4.

Interpretation:
* For an Euler product the local mean of |L(½+it)|² over |t| ≲ 1 is m-dependent. It is lognormal, set by the small-prime values of ξ_m (Selberg's CLT). So disjoint windows are strongly correlated and the moments are heavy-tailed.
* For ψ_α(m³,·) there is no such m-dependent size: different windows behave as independent complex Gaussians.
* This is strong numerical evidence that ψ_α(m³, s) is neither an Euler product nor a Waldspurger-type square root of one (|ψ|² ∝ L(½, Π⊗χ_m) would inherit the lognormal size). That is consistent with the non-uniqueness of Whittaker models on the 4-fold cover.
* Consequence for proof strategies: the Soundararajan/Harper GRH method, which turns GRH into sharp moment bounds via the Euler product, has nothing to act on. A proof of (L³) must exploit cancellation in sums over moduli (the "random" structure), not multiplicative structure. Conversely, the Gaussian behaviour is the heuristic reason S/diag → 1 exactly.

## Round 3 (from ~14:05Z, 9-hour push, (L³) prioritised)

### S42. Zero-density route to the optimal large-sieve shape fails with Gao–Zhao — PROVED (exponent algebra; l3b/zd.py)
Unconditionally, |P(k)| ≪ Z^{β_k−1/2+ε}, with β_k the rightmost zero of L(s,χ_k) in a box of height Z^ε (standard: 1/L ≪ (qT)^ε to the right
of all zeros). So Σ_{N(k)≤H}|P(k)|² ≪ Z^ε·max_σ N(σ,H)·Z^{2σ−1}. With N(σ,H) ≪ H^{A(σ)(1−σ)} and Z = H^{3/4}, this beats the recursion bound
Q^{7/4} iff A(σ)(1−σ) < 15/8 − 3σ/2 for all σ (e.g. A < 2.44 at σ = 0.6, A < 2.75 at σ = 0.7).
Montgomery's zero-detector with the available large sieves, optimised over mollifier and smoothing lengths, gives:

| σ | 0.55 | 0.60 | 0.65 | 0.70 | 0.75 | 0.80 |
|---|---|---|---|---|---|---|
| A, if the optimal quartic LS (H+N) held | 2.10 | 2.20 | 2.31 | 2.41 | 2.53 | 2.65 |
| A, with Gao–Zhao + classical LS | trivial | trivial | trivial | 3.28 | 3.67 | 4.23 |

With the optimal LS the requirement is just met; with the proven Gao–Zhao loss (HN)^{3/4} the density is trivial up to σ ≈ 0.7. So this
route also needs the optimal quartic large sieve, which is false for general coefficients (dFDH 2026).

### S43. A side remark on the small-conductor rows (PROVED)
* On the Möbius side the rows k = k₀j⁴ (j ≤ (H/N(k₀))^{1/4}) all carry the same Möbius sum P(k₀), up to coprimality with j.
* For k₀ = 1 this is the untwisted Möbius sum over Z[i], so any Möbius-side proof of (G), even for M ≥ 1, needs a power-saving prime number
  theorem in Z[i]. The OpenAI/tracker quasi-RH (ζ_{Q(i)} = ζ·L(·,χ₋₄), zero-free for Re s > 0.8749) supplies exactly enough for M ≥ 1.
* For small k₀ ≠ 1 it needs a quasi-RH for the quartic Hecke L-functions L(s,(k₀/·)₄) over Q(i). That is not known; by Theorem Z it is
  implied by (L³).
* On the T-side the same rows are harmless, because ψ_α(m₀³, s) is entire (the K-type twist removes the theta pole) and convexity suffices.

### S44. Multiple-Dirichlet-series type of the (L³) sum versus nearby sums (structural; group theory PROVED, identification HEURISTIC)
* With BBF conventions H(c₁,c₂,c₃) = g(c₁)g(c₂)g(c₃)(c₁/c₂)^{-1}(c₂/c₃)^{-1} for type A3 on the 4-fold cover, summing the outer nodes gives
  Σ_{c₁} g(c₁)(m/c₁)^{-1}N(c₁)^{-s} = ψ(m³, s). So Σ_m g̃(χ̄_m)·ψ(m³,s₁)·ψ(m³,s₂)·N(m)^{-w} is (up to coprimality and the α-twist) an A3
  Weyl group multiple Dirichlet series: a finite group (order 24), meromorphic continuation, asymptotics.
* (L³) is instead Σ_m ψ(m³,s₁)·conj ψ(m³,s̄₂)·N(m)^{-w}. Here conj ψ uses the conjugate-orientation Gauss sums (equivalently R(m̄) = conj R(m)),
  and the middle node carries weight 1, not a Gauss sum. Its functional-equation group is hyperbolic: (2,6,6) per subagent F, and S39's
  subgroup (2,∞,∞). So (L³) is the "wrong" second moment. The Gauss-sum-weighted one is finite type; the absolute square is hyperbolic.
* (L³) would follow from the A3 object if conj ψ_α(m³,½) were a root-number multiple of ψ_α(m³,½) with root number ∝ g̃(χ̄_m)^{-1}. It is not:
  the Kubota functional equation preserves the index and does not conjugate, so |ψ|² ≠ ε·ψ².
* Analogy check for the quadratic family (n = 2): there the k-th moment corresponds to a star diagram with k leaves (D4 finite for k = 3,
  affine D̃4 for k = 4, hyperbolic for k ≥ 5). (L³)'s type (2,6,6) puts it with the quadratic 5th-and-higher moments: no MDS result of this
  type is known.

### S45. No phase locking between ψ_α(m³,·) and the Gauss sum of m (l3/gen_phase.c, l3/phase_stats.py) — NUMERICAL
Hypothesis tested: if conj ψ_α(m³,½) were a root-number multiple of ψ_α(m³,½) with root number ∝ g̃(χ̄_m)^{-1}, then |ψ|² = g̃·ψ², and
(L³) would reduce to the finite-type A3 series of S44. gen_phase.c outputs complex rows R(m) and the row Gauss sums g̃(χ̄_m), computed by
twisted multiplicativity and validated against brute force on all 138 squarefree rows with N(m) ≤ 400 (max error 5e-12).
Statistic |Σ_m w(m)R(m)²|/Σ_m|R(m)|² over squarefree m ≠ 1, for w ∈ {1, g̃, conj g̃, T, conj T, α, ᾱ, g̃², conj g̃²}:
* Z = 1e5 (748 rows): 0.001–0.118, null scale 0.05;
* Z = 1e6 (3482 rows, two windows): 0.003–0.054, null scale 0.024.
All are at noise level (largest 2.2σ among 18 statistics). There is no coupling, consistent with theory: the Kubota functional equation
preserves the index and sends α to ᾱ; it never conjugates. Route closed.

### S46. The precise obstruction: a quadratic toy version of (L³) follows from GRH; the quartic one does not (PROVED, exponent algebra)
Quadratic toy. Over the quadratic family, Σ_{d~D}|Σ_{n~N}μ(n)χ_d(n)n^{-1/2}W|² = D·C + O(N^{1+ε}) with N = D^{3/4}, the same exponents as (L³)
on its Möbius side. Poisson in d maps the off-diagonal to D·Σ_{0<k≤N²/D}ŵ_k(|Σ_{n~N}μ(n)(k/n)n^{-1}|² − Σ_n n^{-2}), because quadratic Gauss
sums are characters times √n. Here the dual has fewer rows (N²/D = D^{1/2}) than its length N. So GRH's pointwise bound
|Σ_{n~N}μ(n)(k/n)| ≪ N^{1/2+ε} gives off-diagonal ≪ D·(N²/D)·N^{-1+ε} = N^{1+ε}. GRH proves the toy.
Quartic case. The Poisson dual of the Möbius side is the T-side (T(c)g(χ̄_c,k) = μ(c)χ_c(k)√N(c) runs backwards): Gauss sums, not Möbius.
* The side with fewer rows than length (rows m ≤ Q = Z^{2/3} < Z) is the T-side. There a pointwise bound would suffice, but its objects ψ_α(m³,·)
  have no Euler product, so GRH says nothing (and S41 shows no hidden Euler product).
* The Möbius side, where GRH applies, has rows H = Z^{4/3} > Z. There pointwise bounds only give the large-sieve shape (S38).
So (L³) is exactly "GRH-strength information for the T-side", i.e. Lindelöf on average for Kubota series. The swap Möbius ↔ Gauss sums under
Poisson, trivial for quadratic characters and nontrivial for quartic ones, is the whole difficulty.

### S47. A naive large-sieve recursion is invalid; the multiplicity floor (l3b/recur.py) — PROVED (counterexample) + lesson
* Rules tried, in log units with g(r) = log_D Δ(D^r, D): Gao–Zhao, classical, row/column symmetry, Poisson in rows (Gauss sums with modulus cc'),
  and Hölder/squaring.
* Their fixed point is g = max(r,1) except at r = 1 (excess ½). That would be the optimal quartic large sieve almost everywhere, which
  contradicts dFDH's non-orthogonality. So at least one rule is misapplied.
* The error: Gao–Zhao (and the classical sieve's "distinct characters") hold for squarefree rows only, while the Poisson identity produces
  all rows k ∈ Z[i]∖0. Rows k = k₀j⁴ repeat the character of k₀, so for constant coefficients Δ(H,D) ≥ H^{1/4}·D, i.e. g(r) ≥ 1 + r/4.
  The fixed point violates this at r < 4/3.
* A correct recursion must carry square-class bookkeeping: quartic × quadratic families k₁k₂², plus the j⁴ multiplicity. That is what
  subagent C's "by square class" recursion did (λ = 1/6, floor 1/12). It remains unverified here.
* For (L³) itself the multiplicity is harmless on the T-side for small m₀ (S43: ψ_α(m₀³,·) is entire and convexity covers N(m₀) ≤ Q^{1/3}).
  For N(m₀) > Q^{1/3} it is part of the main problem.

### S48. Function-field analogue: (L³) is a Weil–Deligne / Selberg-integral statement (ff/ffgauss.py, ff/ffL3.py) — PROVED identities + NUMERICAL
Setting: F_q[T], q ≡ 1 mod 4, c monic squarefree of degree n, quartic symbol (·/c)₄, additive character ψ(x/c) = e_p(res_∞(x/c)).
1. **Prime Gauss sums are explicit (PROVED).** By the residue theorem, res_∞(x/π) = −Tr(x(θ)/π′(θ)). By Hasse–Davenport, the quartic Gauss sum of
   F_{q^k} is −(−G_q)^k. Hence g(χ_π) = χ(−π′(θ))·G_{q^k} = (−π′/π)₄·(−1)^{k+1}G_q^k: the argument of a quartic Gauss sum is the quartic residue
   symbol of the derivative.
2. **Composite form (PROVED from 1 + CRT; NUMERICAL check).** For squarefree c, g̃(χ_c) = μ(c)·u_n·(c′/c)₄ with u_n a unit depending only on n
   (the pairwise CRT cocycles Π_{i≠j}(π_j/π_i)₄ recombine with the (π_i′/π_i)₄ into (c′/c)₄ by the product rule).
   Checked for all squarefree c of degree ≤ 3 over F_5 (spread 1e-15).
3. **Pellet (classical; NUMERICAL check).** μ(c) = (−1)^n·(c′/c)₂, i.e. Möbius is the quadratic character of the discriminant. Checked for all tested
   squarefree c of degree ≤ 4 over F_5. Combining with 2: g̃(χ_c)·(c′/c)₄ = u_n, a Möbius-free formula.
4. **Consequence.** T(c) = μ(c)·conj g̃(χ̄_c) = v_n·χ_c(−1)·(c′/c)₄ exactly (checked: one value of the ratio for all 100 (q = 5) and 2028 (q = 13)
   squarefree cubics). So R(m) = v_n Σ_{c∈A^n(F_q)} χ₄(±Disc c)·χ̄₄(Res(c,m)), a multiplicative character sum of polynomials over affine n-space
   (non-squarefree c drop out automatically since Disc = 0). On the split locus it is a finite-field Selberg/Aomoto sum
   Σ_{θ∈F_q^n} χ₂(Π_{i<j}(θ_i−θ_j))Π_i χ̄₄(m(θ_i)), which Anderson–Evans evaluate exactly as products of Jacobi sums (|·| ≍ q^{n/2}) when m has two roots.
5. **Numerics (large q, fixed degrees n = 3, d = 2).** S/diag = 1.20 (q = 5), 1.077 (q = 13). The excess ≈ 1/q comes from the q non-squarefree rows
   m = (T−a)², where χ_c(m) is quadratic and |R|² ≈ q·diag.
**Meaning.** In function fields (L³) is a statement about character sums of discriminants and resultants, accessible to Weil–Deligne (pointwise
square-root cancellation in the large-q limit, for nondegenerate m) and in the split case even exactly evaluable. Over Z[i] there is no analogue of
item 1: the argument of a quartic Gauss sum is not a residue symbol of anything (Kummer–Patterson problem), and μ is not a character of a
discriminant (no Pellet). These two identities are exactly what makes the Poisson dual Möbius ↔ Gauss sums trivial in F_q[T] and hard in Z[i]
(compare S46). So (L³) over Z[i] is a problem about the joint distribution of quartic Gauss-sum arguments and quartic residue symbols.
- (S48 numerics, ff/ffbig.c, ff/ffbig_runs.txt) Large q at n = 3, d = 2 (q = 29, 37, 41, 53): total S/diag = 1.035, 1.027, 1.024, 1.019. The squarefree
  rows have S/diag ≈ 1/q (0.037, 0.029, 0.026, 0.020), i.e. |R(m)|² ≈ q^{n−1}: beyond square-root cancellation by √q, as expected from the
  Selberg-integral evaluation. The q non-squarefree rows m = (T−a)² have |R|² ≈ q·diag (max row ratio 28.96, 36.97, 40.98, 52.98 ≈ q).
  n = 4, d = 3 (q = 13, 17): 0.988, 0.993; squarefree 0.087, 0.065; max ratio ≈ 0.85q². Fixed q = 5, growing degree (n,d) = (6,4), (7,5):
  0.70, 0.68, with squarefree rows 0.18, 0.22.
  Caveat on the model: in F_q[T] the analogue of α(c) = c/|c| is trivial for monic c, so the function-field T-family is the untwisted family
  (T ∝ g̃). That is why it shows theta-type large non-squarefree rows. The Z[i] α-twisted family has none (S35, S36). So F_q[T] models the
  arithmetic identities, not the absence of theta terms.

### S49. One Poisson step + BGL on the dual: quartic large sieve in unbalanced ranges (l3b/recur3.py) — superseded by the S49 resolution below (PROVED; full proof in l3c/S49_PROOF.md)
Claim (to be verified). For B ≤ A ≤ B², Ξ₄(A,B) ≪_ε (AB)^ε (A + A^{1/3}B + A^{3/4}B^{1/2}). For A ≥ B^{6/5} this is A + A^{3/4}B^{1/2}, i.e. it would
match the dFDH lower bound (Thm 1.2) and prove dFDH Conjecture 1.3 in that range. It does not improve BGL at A = B (both give B^{4/3}).
Sketch.
1. Expand Σ_a W(N(a)/A)|Σ_b β_b χ_b(a)|² and apply Poisson in a modulo bb' (conductor norm ≤ B²). The zero frequency is A‖β‖² (b = b').
   For coprime squarefree b, b': g(h; χ_b χ̄_{b'}) = χ_b(b')χ̄_{b'}(b) g(h;χ_b) ḡ(h;χ_{b'}), and by quartic reciprocity the cross factor is a sign
   depending on residue classes mod 16. Splitting β into those classes gives (A/B²) Σ_{0<N(h)≲B²/A} |Σ_b β_b g(h;χ_b)|². Non-coprime pairs (b,b')
   are handled with Ramanujan sums.
2. Since g(h;χ_b) = χ̄_b(h)g(1;χ_b) for (h,b) = 1 and |g(1;χ_b)|² = N(b), the dual is ≤ (A/B) Ξ^{all}(B²/A, B). Here "all" means all rows h,
   not only squarefree ones.
3. Write h = h₁h₂²h₃³j⁴ and split into dyadic boxes. Bound each box by the cheapest of four estimates: BGL on h₁, the optimal quadratic
   large sieve on h₂, BGL on h₃, and the j-multiplicity (rows j⁴ are trivial characters).
   - The h₁-box gives (A/B)·(B³/A)^{2/3} = A^{1/3}B.
   - The j-box gives (A/B)·(B²/A)^{1/4}·B = A^{3/4}B^{1/2}.
   - recur3.py checks all mixed boxes. Applied to the conjectured truth L(r) these rules reproduce L exactly (consistency). Applied to BGL they give L(r) for r = log A/log B ≥ 1.25.
4. (L³) consequence (to be verified). Rows m = m₀j⁴ with small m₀ need the T-specific convexity input
   max_K min(Q^{1/4}ZK^{9/4}, (Q/K)^{1/4}Z²) = Z^{2+1/15} < Z^{25/12}. Combined with step 3 through duality, this would give the unconditional
   Σ_{N(m)≤Q}|ψ_α(m³,½+it)|² ≪ Q^{13/8+ε}(1+|t|)^{O(1)}. That improves Q^{7/4} (round 2), Q^{5/3} (BGL) and Q^{15/8} (Gao–Zhao), but is far from (L³) = Q^{1+ε}.
   Q^{13/8} is the dFDH lower bound A^{1/2}B^{3/4} at (Q, Q^{3/2}), so no large-sieve argument can go below it. (L³) needs T-specific cancellation (S43, S46).
Open verification items (why this is CANDIDATE, not PROVED):
   (a) Step 1 for b, b' sharing prime factors. Also check that the sign split is uniform in h, including h divisible by λ and units.
   (b) Whether BGL Thm 1.3 applies to the non-squarefree row boxes after factoring χ_b(h₂²) = (χ_b)²(h₂). The quadratic sieve then needs the
       column coefficients to be twisted by χ_b(h₁h₃³), which depend on the row. So each box needs Cauchy–Schwarz over h₁, h₃ before the h₂ sieve.
       recur3 assumes this costs nothing beyond the multiplicities x1 + x3 + y. That has to be checked carefully.
   (c) Literature: dFDH write that "any power-saving improvement on the (AB)^{2/3+ε} term … for n ≥ 4 would be a substantial advance". That
       suggests either the remark is meant only at A ≍ B, where this candidate gives nothing, or that (b) fails. Heath-Brown's cubic sieve
       [HB00] uses the same Poisson recursion, so check whether [HB00]/[BGL14] already record unbalanced-range bounds of this shape.
       BGL14 is not yet in lit/; fetch it next.
S49 verification progress (22:13Z, lit/bgl/bgl.txt = BGL14, arXiv 1112.1650, §3):
- Step 1 is exactly BGL Lemma 3.3: B₃(M,N) ≪ M + (M/N) max_{K≤N²/M} B₂(K,N). Lemma 3.2 handles non-coprime b,b' (d = (b₁,b₂), Möbius over r | sd).
  Lemma 3.1/(3.4) is the reciprocity and duality step. So items (a) and the Poisson normalisation are covered by BGL as published.
- Step 3 refines BGL Lemma 3.4. BGL decompose a = a₁a₂²…a_n^n but apply the sieve only in a_ℓ, ℓ ∈ {1,2}, and they bound X/X_ℓ by
  min(M^{1/3}, M^{2/3}/X_ℓ). That bound loses M^{1/3}N where the j-box truth is M^{1/4}N, and it never uses the sieve in a₃. For n = 4 and j ∈ {1,3},
  jℓ ≢ 0 (mod 4) also holds for ℓ = 3. So their own argument (λ*_b = λ_b χ_b(a/a_ℓ^ℓ), |λ*| ≤ |λ|) allows ℓ = 3, and for j = 2 the GL quadratic sieve is optimal.
  The row-dependent twist in (b) is handled exactly as in BGL: absorb it into λ* and sum the multiplicity over the other a_i.
  So recur3's c1, c2, c3 are BGL's (3.6) with ℓ = 1, 2, 3. c4 (fix a₄, classical large sieve over the characters (·/a₁)(·/a₂)²(·/a₃)³,
  distinct triples ↦ distinct primitive characters) is a standard addition.
- Why BGL did not see it. At M ≍ N their fixed point is (MN)^{2/3}, and the M^{1/3}N loss in Lemma 3.4 is harmless there. The gain only appears
  for M > N, where the new bound is M + M^{1/3}N + M^{3/4}N^{1/2} < (MN)^{2/3} for every r = log M/log N > 1.
  Item (c) therefore reads: dFDH's "substantial advance" remark concerns the balanced range, which this does not touch.
- Status: still CANDIDATE. Remaining before PROVED:
  (i) write the n = 4 version of Lemma 3.4 with ℓ ∈ {1,2,3} plus the c4 box, and the max-min over boxes as an explicit inequality chain
      (recur3 is a numerical grid check, not a proof);
  (ii) the (L³) transfer with T-specific m₀j⁴ rows (S49 item 4);
  (iii) re-derive the dual range K ≤ (MN)^ε N²/M including the s | c^∞ factor in BGL Lemma 3.3.

S49 resolution (22:18–22:33Z) — PROVED. Full write-up in l3c/S49_PROOF.md.
- Lemma A (sharpened BGL Lemma 3.4, n = 4): B₂(K,N) ≪ (KN)^ε(K + K^{1/4}N + (KN)^{2/3}) for all K, N.
  The new inputs are the classical large sieve on mixed boxes (rows a₁a₂²a₃³ give distinct primitive characters of conductor a₁a₂a₃) and the
  sieve in a₃. The proof is a 4-case analysis; the key step is min(X₄Nσ, X₄P²) ≤ X₄N^{2/3}(mσ²)^{2/3} ≤ K^{2/3}N^{2/3}, by rearrangement.
- Theorem 1: Ξ₄(A,B) ≪ (AB)^ε(A + B + A^{3/4}B^{1/2} + A^{1/2}B^{3/4} + min(A^{1/3}B, AB^{1/3})).
  This is BGL Lemmas 3.1–3.3 as published plus Lemma A. It proves dFDH Conjecture 1.3 when max ≥ min^{6/5}, where it is sharp.
- Lemma B (Poisson in the column variable for box-structured columns m₁m₂²m₃³): the off-diagonal is (Y/P)·X₂X₃·B₂(P²/Y, X₁).
  The cross factors of τ(ψ̄_mψ_{m'}) split into functions of (m₁,n'), (m₁',n) and (n,n'), so Cauchy over n = (m₂,m₃) is taken after Poisson.
  The case of coprime radicals is written out. Shared primes follow the BGL Lemma 3.2 pattern (STANDARD; each term is shown dominated by N(d)^{-1/3}).
- Theorem 2 (T-family, unconditional): Σ_{N(m)≤Q}|R_Y(m)|² ≪ (QY)^ε(Y + QY^{1/3} + Q^{1/2}Y^{3/4})·Y for Y ≤ Q^{5/3}.
  Hence Σ_{N(m)≤Q}|ψ_α(m³,½+it)|² ≪ Q^{13/8+ε}. Inputs: Lemma A, Lemma B and Kubota convexity.
  The convexity input is applied to index ν₀d², via T(dc') = ±χ_{c'}(d²)T(d)T(c'). Its N(d) cancels, so |R_Y(m)| ≪ Y^{1/2+ε}N(m₁³m₂²m₃)^{1/4}.
  Case I (P² ≤ Y) uses the classical sieve or convexity. Case II (P² > Y) uses Lemma B with Lemma A, which gives four exponent inequalities.
  Grid confirmation: l3c/boxes_L3b.py and l3c/boxes_Y.py; the worst box is the main box, at 1/2 + 3η/4.
- Bugs found and fixed while checking (logged honestly):
  - p4 first omitted the j⁴ multiplicity of dual rows h and "beat" the dFDH lower bound (Q^{1.6}). With H^{1/4}·Q′ it is correct.
  - The c1 row-sieve alone gives only Q^{5/3}, because of the m₁m₂² box (Q^{2/3}, Q^{1/6}). Column Poisson (Lemma B) is what removes it.
- Consequence for θ: none. Q^{13/8} is the large-sieve limit (dFDH), and §7/S7 need (G) beyond the optimal shape near M = 1.
  At M = 1/2 (η = 2) Theorem 2's Case I leaves a Q^{y} loss for j⁴-rows with 0 < y < 1/18, and the Möbius-side j⁴ rows need quartic zero-free regions (S43).

### S50. Beyond Q^{13/8}: the obstruction is exactly the T–theta correlation (l3d/S50_REDUCTION.md, l3d/reduction.py, l3/lin_twist.py) — reduction PROVED (STANDARD details), input (P4) OPEN/NUMERICAL
* **Route.** Cauchy route S² ≤ Z‖A*R‖², then Poisson in c with γ = R̄. The dual rows h = h₀j⁴ give X(h₀,j) = Σ_m R̄(m)g̃(χ̄_m)ψ_m(h₀).
* **Double reflection** (Kubota functional equations in m with ℓ = 0, then in c with ℓ = 1):
  X = Q^{3/4}V̂(3/4)·𝒞(h₀) + X_dual, where 𝒞(h₀) = Σ_c T(c)W·τ₄(h₀c) (quartic theta coefficient) and |X_dual| ≪ Z^{1+ε}N(h₀)^{1/2} (small h₀).
* **Bookkeeping.** With the Cauchy/Rankin–Selberg bound 𝒞 ≪ Z^{7/8} this returns exactly Q^{13/8}, the dFDH large-sieve limit, consistent with S49.
  (P4) 𝒞 ≪ Z^{3/8+ε} gives Q^{1.543}. Any power saving over Z^{7/8} improves 13/8.
* **Numerics.** |X(h₀)|/√S = O(1) for h₀ = 1 and nine twists at Z = 10⁵, 10⁶; Cauchy-size 𝒞 would give ≈ 10³.
* **Theory.** Suzuki's τ₄ results vanish or evaluate at indices m³δ and m²δ, but the duality forces index exponent 1.
  Averaging over h₀ helps only for N(h₀) ≥ Z^{0.03}, while bounded h₀ alone carry 13/8.
* **Correction to S44.** In BBF orientation the needed double series pairs Gauss sums of opposite orientation, so it is not a standard
  (W-invariant) A2 WMDS. S44's A2/A3 identification needs its orientation re-derived.
* **Interpretation.** (P4) is a quartic Heath-Brown–Patterson-type equidistribution: the relative sign of T(π) and τ₄(π), both square
  roots of Gauss-sum quantities under Patterson–Eckhardt.

### S51. The S50 obstruction series is a G₂ WMDS (l3d/S51_G2.md, lit/cg) — structure PROVED (FE group), identification with Chinta–Gunnells G₂ CHECKED, analytic extension to the α-twist STANDARD/unwritten
* The two Kubota FEs (index d³ in c, index h₀c in d) generate W(G₂) (order 12). This matches CG's G₂ series with symbol ε³:
  c is the short root (χ̄), d the long root (χ), coupling ε^{−1}.
* With the α-twist only the theta divisor survives. Its residue r(s) = Σ_c T(c)τ₄(h₀c)N(c)^{−s} has an FE about s = 3/8.
* Convexity gives X ≪ Q^{3/4}Z^{3/8}h^{B} + (QZ)^{1/2}h^{B′}, hence the (L³) sum ≪ Q^{3/2+ε}. This is conditional on CG extended to the
  archimedean twist ℓ = 1 and to Q(i) ⊉ μ₈, with standard growth bounds.
* Numerics: the smooth first moment is real and of square-root size. The "negligible remainder" prediction is withdrawn: theta_const.py
  shows no constancy at Q ≈ 10³–10⁴.
* Classification: index exponent r = 1, 2, 3 gives second-moment Coxeter groups (2,3,3), (2,4,4), (2,6,6). (L³) is hyperbolic;
  S44's A3 claim was for the index-m family.

### S51 addendum (23:57–00:08Z): §17–§20
* §17: DDHL (4.28) confirms the twist rule and the W(G₂)-invariant Gamma normaliser (root twists 0,1,1,2,3,3).
* §18–§19: μ₈ ⊄ Q(i). An exact check of CG's local G₂ action for n = 4 shows raw Gauss sums fail at N ϖ ≡ 5 mod 8. The Brylinski–Deligne
  sign χ_ϖ(−1)^{C(i,2)} restores every relation; globally it becomes χ_c(−1) on the c-node, absorbed into Ψ.
* §20 (correction): CG's normalising L-factors L(·,α⁴) and L(·,α¹²) sit on their critical line at the residue centre.
  * The α⁴ factor is removed by convolution with T*.
  * The α¹² factor is not, so the CG-conditional bound is Q^{1.51}, and Q^{3/2} additionally needs that factor (GRH for one Hecke
    L-function suffices).

### S52 (00:50–01:35Z, 16 Oct): mean-value Bochner — an unconditional route below 13/8 (l3f/S52_MEANVALUE_BOCHNER.md, S52_PROOF.md)
* **Idea.** Bochner's tube theorem with mean-value regions instead of absolute convergence. The regions are
  * LoA in d (quartic large sieve);
  * S49-type mean squares in c;
  * blockwise "best of" the large sieve over c and the mean square over d, with c-FE reflection (fix (a)).
  They are taken over all 12 G₂ chambers, with DDHL normalisers restricting evaluation to x, y ≥ 0.
* **Growth.** An envelope-of-holomorphy sup argument gives polynomial growth and convex h₀-exponents. The vector-valued version gives
  RMS-in-h₀ bounds, which is what S50 needs.
* **Corrections made along the way.**
  * S50's j⁴-coprimality turns twists into h₀e² (§10).
  * The true quartic large sieve is needed in h₀ (§11).
  * Joint twist exponents are not additive (§11).
  * Non-main row boxes are handled by c-FE reflection (§8–§9).
* **Result.** A recursion on the T*-family exponent table S(q) converges to (L³) ≪ Q^{1.5156+ε} (§12), modulo Lemma 4 bookkeeping,
  which has an argument but untracked constants. Finer-grid and e-tail-cut runs are in progress.
* **Limit of the method.** The j⁴-multiplicity at large effective twists. In particular S(Y,Y) ≈ Y^{2.25}, rows = length, which caps the
  residue-line point. Below Q^{3/2} the Cauchy route is impossible, and the second-moment MDS hull gives nothing (m-FE → 1/L).
