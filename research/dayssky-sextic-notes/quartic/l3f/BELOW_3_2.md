# Below Q^{3/2}: what a proof of (L³) has to supply (DaysSky with Claude assistance, 16 Oct 2026)

Notation as in L3_REPORT §10b: Z = c-length, Q = Z^{2/3}, P(k) = Σ_c μ(c)ξ_k(c)W(N(c)/Z) (Möbius sums with quartic Hecke characters
ξ_k), and S(Q′,Z) = Z^{s(q′)} for Q′ = Z^{q′}.

## 1. Dual identity (L3_REPORT §10b, steps 1–3; PROVED up to the gcd/Ramanujan bookkeeping written there)
  S(Q,Z) − c₀QZ = (Q/Z)·Σ_{N(k) ≲ Z²/Q} (|P(k)|² − C_k·Z).
(L³), i.e. S ≪ QZ^{1+ε}, is therefore equivalent to the k-variance being ≪ Z^{2+ε}.

## 2. Which k-ranges are free (PROVED, except where marked)
* For a dyadic range k ∼ K′ with K′ ≤ Z^{5/6}, S49 Theorem 1 (the dFDH shape when max ≥ min^{6/5}) gives
  Σ_{k∼K′}|P(k)|² ≪ (K′ + Z)Z^{1+ε} ≪ Z^{2+ε}. The main term K′Z ≤ Z² is of the same size, so these ranges cost (Q/Z)Z^{2+ε}log Z =
  QZ^{1+ε}, which is admissible.
* For K′ ∈ [Z^{5/6}, Z], Gao–Zhao's (K′Z)^{2/3} term gives only Z^{7/3}. These ranges are free under the dFDH conjecture.
* The real difficulty is Z < N(k) ≤ Z^{4/3}. There the main term K′Z exceeds Z², and an asymptotic with relative error Z/K′ is needed:
  the variance must equal its diagonal up to Z².

## 3. That range is self-dual (PROVED as exponent bookkeeping)
* Poisson in k on a range k ∼ K′, followed by T(c)g̃(χ̄_c) = μ(c) in reverse, turns the k-variance into (K′/Z)·(S_{Q′} − c₀Q′Z) with
  Q′ = Z²/K′ ∈ [Q, Z]. Here S_{Q′} is the T-family second moment with Q′ rows, with a different smooth weight.
* So (L³) at q = ⅔ needs Lindelöf on average, S(q′) ≤ 1 + q′, at every q′ ∈ [⅔, 1]. The q′ = ⅔ piece is the original problem.
* What the S52 method achieves on this whole interval is s(q′) = q′ + 4/3 (S52_LIMIT_PROOF for q′ ≤ 5/6, S49 above). That is exactly
  Z^{1/3} above Lindelöf at every q′.

## 4. Consequence for methods (PROVED for the class described)
* Every argument built only from these inputs stops at s(q′) = q′ + 4/3, i.e. at Q^{3/2}:
  * Cauchy–Schwarz plus Poisson (S50);
  * first-moment multiple Dirichlet series with mean-value regions (S52);
  * the second-moment hull (second.py: its binding constraint at the (L³) direction q = ⅔ is invariant under the c-functional equations);
  * large sieves, including GRH for each P(k) (L3_REPORT §10b).
* The reason: each of these inputs is invariant under the duality of §3. Q^{3/2} is also the GRH value, and the value the conjectured
  optimal quartic large sieve (de Faveri–Dunn–Hoffstein) would give.

## 5. What would break the duality (candidates, not yet attempted in full)
* **(i) Spectral input.** Write the off-diagonal of the k-variance through metaplectic Kloosterman sums on the 4-fold cover of
  GL₂(Q(i)), and use a Kuznetsov-type formula. The residual spectrum (the quartic theta function) should give the main term, and
  cancellation would come from spectral large sieve bounds for the cuspidal part. This is the analogue of how Conrey–Iwaniec and
  Petrow–Young break self-dual cubic moments.
* **(ii) A ratios-type asymptotic for Σ_k |1/L(½+it, ξ_k)|²** with power saving Z^{1/3} in the thin quartic family (the conditional
  route of L3_REPORT §10). Unconditionally even the complete-family Möbius variance (Barban–Davenport–Halberstam / Montgomery–Hooley
  for μ) has only logarithmic savings, so this route needs new zero-density input. By Theorem Z, (L³) itself implies a 7/8 zero-free
  region for the family.
* **(iii) Type I/II decomposition of μ in the small-conductor form R(m) = Σ_c μ(c)χ_c(−1)g̃(m;χ_c)** (conductor N(m), length
  Z ≫ N(m)^{1/2}), in the style of Heath-Brown–Patterson.
  * Type I pieces are long Kubota sums: theta main terms τ(ma²) plus negligible duals. Their sum over a must cancel; this is the
    T–theta correlation of S50.
  * Type II pieces are trilinear in (m, a, b) with the quadratic coupling (b/a)₂. Bounding them by a large sieve for sparse modulus
    sets {b²m³} is again self-dual, unless the quadratic coupling supplies extra cancellation.

## 6. Literature check (16 Oct 2026; arXiv texts read in scratchpad/dl_arxiv1)
* **de Faveri–Dunn–Hoffstein, "Non-orthogonality of the cubic and quartic large sieves via Rankin–Selberg"** (arXiv 2607.07911,
  8 Jul 2026).
  * They use Σ_{N(a)≤A}|ψ⁽⁴⁾(a, 1+ε+it)|² ≪ A^{1+ε}(1+|t|)^R, over all numerators a, proved by Zagier's Rankin–Selberg regularisation.
    They remark that BGL's large sieve gives it too.
  * That family is not thin (size A > √conductor A^{1/2}).
  * Restricted to a = m³ by positivity it gives only Q^{3+ε}, so it does not touch (L³).
  * Their theorem says the quartic large sieve is not perfectly orthogonal, the excess coming from Gauss-sum (theta) bias. The dual form
    of (L³) in §1 is exactly a statement that the Möbius sequence is orthogonal with power saving Z^{−1/3}: the variance equals its
    diagonal up to Z². §10b's observation that μ removes the theta term is consistent with that.
* **Diaconu–Ion–Paşol–Popa, "On the second moment … of Hecke L-functions of r-th order characters"** (arXiv 2607.27131, 29 Jul 2026).
  * They prove asymptotics for Σ_a μ²(a)|L(½,χ_a)|², which are r-th order Hecke L-functions, not Kubota series.
  * They state that Lindelöf for the Kubota-series family is conjectural.
  * Their error term is limited by a large-sieve barrier that GRH would remove: the same kind of barrier as §4 here.
* **Classical analogues.** Searches found no improvement on either:
  * the single-modulus sixth moment of Dirichlet L-functions beyond q^{3/2} (CIS and Chandee–Li–Matomäki–Radziwiłł need an extra average
    over moduli);
  * Heath-Brown's twelfth-moment bound T^{2+ε}.
  (L³) has the same shape: family size Q, conductor Q³, so the large sieve gives Q^{3/2}.
* **Conclusion.** No published unconditional input found so far breaks the §3 duality for a thin family of this shape.

## 7. The hyperbolic second-moment series reaches only convexity (06:05Z; l3h/tits.py) — PROVED as geometry, NUMERICAL as LP
* **The series.** Z(s₁,s₂,w) = Σ_m ψ(m³,s₁)·conj ψ(m³,s₂)·N(m)^{−w}, in centred coordinates u = (s₁−½, s₂−½, w−½).
* **Its functional equations.**
  * Two Kubota FEs with index m³ shift w by 3u_i, so c(3,i) = −3.
  * The Hecke FE of L(w, χ_{c₁}χ̄_{c₂}) shifts s₁ and s₂ by u₃, so c(i,3) = −1.
  * The group is the hyperbolic Coxeter group (2,6,6) (orders checked).
* **Where (L³) sits.** The target point (0,0,½) lies in the interior of the Tits cone, since its stabiliser ⟨σ₁,σ₂⟩ is finite. So a
  continuation to the Tits-cone interior would in principle reach it.
* **What Bochner actually reaches.** Starting from absolute convergence C₀ = ρ + {u > 0}, with ρ = (½,½,½), Bochner reaches
  conv(W·C₀) = conv(W·ρ) + X.
* **Computation.** An LP over all 5634 group elements of length ≤ 16 gives min{u₃ : (0,0,u₃) ∈ conv(W·ρ) + X} = 2. It is already attained
  by ρ and σ₁σ₂ρ, and does not move between length 10 and 16.
* **Meaning.** The second-moment series is reached on the central line only for Re w ≥ 5/2. That is Σ_{N(m)≤Q}|ψ|² ≪ Q^{5/2},
  plain convexity.
  * Unlike a finite group, the hyperbolic orbit W·ρ does not surround the origin. So the finite-group mechanism (whose hull of W·ρ
    contains 0) has no hyperbolic substitute, even granting exact functional equations.
  * This quantifies L3_REPORT §10c.
