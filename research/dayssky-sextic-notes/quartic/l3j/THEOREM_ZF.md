# Theorem ZF: any bound below Q^{3/2} for (L³) gives a zero-free strip for the quartic Hecke family over Q(i)

DaysSky with Claude assistance, 16 Oct 2026.

Labels:
* PROVED: argued in full here.
* STANDARD: routine and sketched only.
* CITED: taken from the named file, where it was checked numerically or proved.

## Statement

For θ′ ≥ 0 let (L³)_{θ′} be the bound

  Σ_{N(m)≤Q} |ψ_α(m³, ½+it)|² ≪_ε Q^{1+θ′+ε}(1+|t|)^A  for all Q ≥ 1, t ∈ R.

* (L³) is (L³)₀.
* THEOREM_Q32 proves (L³)_{1/2}.

For k₀ ∈ Z[i]∖0 and e ∈ {0,1}, let ξ = ξ_{k₀,e} be the Hecke character c ↦ (k₀/c)₄·λ(c)^e on primary c coprime to 2. Here
λ(c) = (−1)^{(N(c)−1)/4} = χ_c(−1).

(G)^{tw}_{θ′} is the squarefree large-sieve form of the same bound (PROBLEM.md's (G)), uniform in index twists h and for all scales:

  Σ_m V(m/√Q)|Σ_{c sqfree} T(c)W(N(c)/Z)χ̄_c(mh)|² ≪ Z·Q^{1+θ′+ε}·N(h)^ε  for all Z, Q, h,

for the narrow-annulus weights W of the proof.
* It follows from the Dirichlet-series form for the squarefree Kubota series by Mellin.
* It follows from (L³)^{tw} for ψ_α (all c) given a local 2×2 inversion at the primes dividing m. That inversion is not done here; see
  THEOREM_ZF_PROOF Lemma 1.
* Lindelöf for each ψ^{sf}_α(ν,·) implies it.
* Every unconditional method in these notes (S49–S53, "twisted main-box problems" in S52_FULLSUM) treats such twisted families, so this
  is the form in which a proof of (L³) would naturally come.

**Theorem ZF.** Assume (G)^{tw}_{θ′} with 0 ≤ θ′ < ½. Then for every k₀ and e, L(s, ξ_{k₀,e}) has no zeros in

  Re s > σ(θ′) := (7 + 10θ′)/(8(1 + θ′)) < 1.

**Values of σ.**

| θ′ | 0 | 1/8 | 1/4 | 3/8 | ½ − δ |
|---|---|---|---|---|---|
| σ(θ′) | 7/8 | 11/12 | 19/20 | 43/44 | 1 − δ/6 + O(δ²) |

**Corollaries.**
1. (L³) implies Re ρ ≤ 7/8 for every zero of every quartic Hecke L-function L(s, (k₀/·)₄λ^e) over Q(i).
   * For k₀ ∈ {±1, ±i} these are ζ_{Q(i)} and Dirichlet L-functions mod 16. There 7/8 is OpenAI's quasi-RH theorem (Sept 2026).
   * For general k₀ (e.g. k₀ = 1+i, 2, 3, 1+2i) they are Kummer and dihedral Hecke L-functions, which that theorem does not cover
     (PROBLEM.md). No zero-free strip is known for them.
2. Any bound Σ_{N(m)≤Q}|ψ_α(m³,½)|² ≪ Q^{3/2−δ} with δ > 0 (and polynomial t-dependence) implies a zero-free strip of width about δ/6
   for the whole family.
   * So the Q^{3/2} reached by every method so far (S49–S53, the GRH argument of L3_REPORT §10b, BELOW_3_2) is the limit of any method
     that does not also prove a quasi-Riemann hypothesis for quartic Hecke L-functions over Q(i).
   * This turns BELOW_3_2 §4 (the "self-duality" barrier) into a theorem about zeros.

3. **The square-case large sieve is just as hard.** Assume only the optimal large-sieve SHAPE for T-coefficients:
   S(Q,Z) ≪ (Q+Z)Z^{1+ε} for all Q and Z, uniform in index twists. This is L3_REPORT §12 sub-problem 1, listed there as "most approachable",
   and §10b proves it on GRH.
   * Running the proof at x = 1 gives 𝒜 ≪ (Z + Z²/Q)Z^ε at Q = Z.
   * The gcd pieces, at row scale Q/N(h)² and length Z/N(h), are bounded by the same shape. Their contribution is
     Σ_h N(h)^{−1}(Z/N(h) + Z²/Q), which is fine.
   * So that sub-problem alone implies zero-freeness in Re s > 7/8 for every quartic Hecke L-function over Q(i). It is not a "bonus" but a
     quasi-RH-strength statement, sitting between GRH and QRH(7/8) in strength.

## Proof (summary; the complete argument with all lemmas is in l3j/THEOREM_ZF_PROOF.md)

Notation:
* W is a smooth bump supported in N(c)/Z ∈ [1, e^η], with η > 0 small and fixed.
* R(m) = Σ_c T(c)χ̄_c(m)W(N(c)/Z) over squarefree primary c coprime to 2, with T(c) = α(c)g̃(χ̄_c).
* Rows are all m ∈ Z[i], weighted by V(m/√Q), where V ≥ 0 is the radial Bessel potential with V̂(ξ) = (1+|ξ|²)^{−B}, B = 2.
  V is continuous, positive and exponentially decaying.

**Step 1 (T-side bound). Direct from (G)^{tw}_{θ′}. The derivation from the Dirichlet-series form (Mellin, plus local inversion) is sketched below; see the caveat in THEOREM_ZF_PROOF Lemma 1.**
* For N(c) in the support, T(c)χ̄_c(m) = α(c)g̃(χ̄_c, m³) (CITED: PROBLEM.md identities).
* So R(m) is a Mellin integral of ψ_α(m³, ½+it) against W̃, apart from the non-squarefree c in ψ_α. Those terms are handled as follows.
  * Write c = c_sf·d, where d is the powerful part. By CRT, g̃(χ̄_c, m³) = (real cross sign)·g̃(χ̄_{c_sf}, m³d²)·g̃(χ̄_d, m³).
  * g̃(χ̄_{π^e}, m³) = 0 for e ≥ 2 unless π^{e−1} | m³. So for each m only the d with rad(d) | m and e ≤ 3v_π(m) + 1 occur: at most
    τ(m)^{O(1)} of them.
  * Each local factor satisfies |g̃(χ̄_{π^e}, m³)|·N(π^e)^{−½} ≤ 1. Example cases: e = 2 gives ≤ N(π)^{−½}; e = 4 (a Ramanujan sum)
    gives ≤ 1.
  * Inverting this finite triangular system expresses the squarefree part through ψ_α(m³h_d, ·) for at most τ(m)^{O(1)} indices h_d,
    with bounded coefficients.
  * Cauchy–Schwarz and the Dirichlet-series form on the rows m = rad(d)·m′ then give the bound below. The sum over d costs Q^ε, provided
    the local 2×2 systems are boundedly invertible on the line; that is not checked.
  * This argument is sketched here, not written out in full.
* Cauchy–Schwarz in t gives |R(m)|² ≤ Z·(∫|W̃|)·∫|ψ_α(m³,½+it)|²|W̃(½+it)|dt. Summing against V and using (L³)_{θ′} on dyadic shells:

  S_V := Σ_m V(m/√Q)|R(m)|² ≪_η Z·Q^{1+θ′+ε}.

* The same holds for each weight W·(u/η)^n, with u = log(N(c)/Z) and constants polynomial in n.
* It also holds for the class-restricted rows ½(R(m) + R(−m)). These come from replacing T(c) by ½(1 + λ(c))T(c), because
  λ(c)χ̄_c(m) = χ̄_c(−m).
* From now on c is restricted to one class N(c) ≡ 1 (mod 8), or to the class N(c) ≡ 5 (mod 8).

**Step 2 (exact dual form). PROVED algebra; Gauss-sum identities CITED from L3_REPORT §10 step 3.**
* Open the square and Poisson-sum over m ∈ Z[i], which is self-dual for e(Im ·).
* For a coprime pair c₁ ≠ c₂ with q = c₁c₂, CRT gives g(k; χ̄_{c₁}χ_{c₂}) = χ̄_{c₁}(c₂)χ_{c₂}(c₁)g(k;χ̄_{c₁})g(k;χ_{c₂}).
* Two identities simplify this:
  * T(c)g(k;χ̄_c) = μ(c)χ_c(k)√N(c);
  * g(k;χ) = χ(−1)·conj g(k;χ̄).
* The overall sign is (−1)^{e₁e₂ + e₂}, with e = (N(c)−1)/4 mod 2. It equals +1 when c₁, c₂ lie in the same class mod 8.
  * Quartic reciprocity gives (−1)^{e₁e₂}.
  * χ_{c₂}(−1) gives (−1)^{e₂}.
* Hence the (c₁,c₂,k) term is
  Q·a_{c₁}(k)·conj a_{c₂}(k)·(1 + Q N(k)/(N(c₁)N(c₂)))^{−B},  where a_c(k) = μ(c)χ_c(k)W(N(c)/Z)/√N(c).
* Pairs with gcd g > 1 are not handled through their exact dual form. Step 4 instead treats them on the T-side, plus a formal form at
  length Z/N(h).
* For the record, the exact dual of a gcd-g pair, c_i = g·a_i, has the structure r_g(k)·μ(a₁)χ_{a₁}(gk)·conj(μ(a₂)χ_{a₂}(gk)).
  * The sign there is (−1)^{e(g)(e(a₁)+e(a₂)) + e(a₁)e(a₂) + e(a₂)}.
  * Inside one class mod 8 we have e(a₁) = e(a₂) = f, so the sign equals (−1)^{f(f+1)} = +1.
* Numerical check (l3j/signcheck.py): for coprime pairs, T(c₁)·conj T(c₂)·g(k; χ̄_{c₁}χ_{c₂}) = (−1)^{e₁e₂+e₂}·μμ′χ_{c₁}(k)·conj χ_{c₂}(k)·√N(c₁c₂).
  156 random (c₁,c₂,k) with N(c_i) ≤ 90 were checked by brute force over residues mod c₁c₂, with 0 mismatches. All four sign
  patterns (e₁,e₂) occurred.

**Step 3 (approximate positivity in a narrow annulus). PROVED.**
* Write N(c_i) = Ze^{u_i} with u_i ∈ [0, η], and r_k² = QN(k)/Z². The weight is G_k(u₁+u₂) with
  G_k(s) = (1 + r_k²e^{−s})^{−B}.
* All derivatives of log G_k are bounded by C_j·B, uniformly in r_k. The first two are g₁ ∈ [0,B] and g₂ ∈ [−B/4, 0].
* Taylor-expand log G_k at s = 0 and split exp(g₁(u₁+u₂) + g₂(u₁²+u₂²)/2) into product factors. What remains is
  exp(g₂u₁u₂ + R(u₁+u₂)) = 1 + H(u₁,u₂), where H is a real-analytic symmetric function with no constant or linear term, so H = O(η²).
* Expanding H in monomials u₁^a u₂^b and using |xȳ| ≤ ½(|x|²+|y|²), the k-th dual form satisfies

  D_k ≥ G_k(0)·[(1 − O(η³))|A₀(k)|² − Σ_{a≥1} c_a |Ã_a(k)|²],  with Σ c_a = O(η²).

  Here A₀ and Ã_a are μ-side sums with weights W·e^{g₁u+g₂u²/2}·(u/η)^a.
* Expand e^{g₁u+g₂u²/2} = Σ_n e_n(k)(u/η)^n, with e₀ = 1 and Σ_{n≥1}|e_n| = O(Bη). This rewrites everything through the k-independent
  sums B_n(k) = Σ_c μ(c)χ_c(k)W(N(c)/Z)(u_c/η)^n/√N(c).
* Let 𝒜 := max_n (n+1)^{−A−2} Σ_k G_k(0)|B_n(k)|². Running the inequality for each n, i.e. for T-side weight W(u/η)^n, gives

  𝒜 ≤ C·max_n (n+1)^{−A−2} Σ_k D_k^{(n)} + O(Bη + η²)·𝒜.

* Fix η ≪ 1/B so that the last coefficient is ≤ ½. Then 𝒜 ≤ 2C·max_n(…).

**Step 4 (sum the dual forms). PROVED, by induction on Z.**
* Let D_k be the "formal" Hermitian form: the coprime-pair formula of Step 2 summed over all pairs (c₁,c₂), which is what Step 3 bounds
  from below. Then

  Σ_k D_k = (1/Q)·S_V^{cop} + Σ_k D_k^{noncop}.

  * S_V^{cop} is the T-side sum restricted to coprime pairs (c₁ ≠ c₂); the c₁ = c₂ terms are included in D^{noncop}.
  * D_k^{noncop} is the formal form restricted to pairs with gcd > 1.
* **The T-side.**
  * S_V^{cop} = S_V + Σ_{h>1} μ(h)·S^{(h|)}, where S^{(h|)} = Σ_{(m,h)=1} V(m/√Q)|Σ_{h|c} T(c)χ̄_c(m)W|² ≥ 0.
  * By twisted multiplicativity, T(ha)χ̄_{ha}(m) = T(h)χ̄_h(m)·T(a)χ̄_a(h²m)·(class sign). So S^{(h|)} is the T-family at length Z/N(h),
    evaluated at the rows h²m.
  * By (G)^{tw}, S^{(h|)} ≪ (Z/N(h))·Q^{1+θ′+ε}N(h)^ε. Summing over N(h) ≤ Z gives S_V^{cop} ≪ Z^{1+ε}Q^{1+θ′}.
* **The formal non-coprime part.**
  * By inclusion–exclusion over h | gcd, Σ_k D_k^{noncop} = Σ_{h>1} −μ(h)·F_h.
  * F_h is the formal form over pairs both divisible by h. Since |χ_h(k)|² = 1_{(k,h)=1}, it is (1/N(h)) times a Step-3 form for the
    Möbius sums Σ_{(a,h)=1} μ(a)χ_a(k)W of length Z/N(h), over k coprime to h. Dropping the k-restriction only adds nonnegative terms,
    up to the O(η²) of Step 3.
  * The a-restriction is a twist: 1_{(a,h)=1}χ_a(k) = χ_a(kh⁴). On the T-side it corresponds to the family
    T(a)1_{(a,h)=1}χ̄_a(m) = T(a)χ̄_a(mh⁴), i.e. the twisted index m³h¹².
  * So the induction runs over the class of twisted families with indices m³h′, all controlled by (G)^{tw} with N(h′)^ε uniformity.
    The diagonal for frequencies coprime to h is ≤ (Z²/Q)·η/N(h).
  * By the induction hypothesis at length Z/N(h), with the same frequency weights, F_h ≪ (1/N(h))·[(Z/N(h))Q^{θ′} + Z²/Q]·Z^ε. Summing
    over h costs a factor log² Z.
* **The c₁ = c₂ terms** contribute Σ_k G_k(0)Σ_c|a_c|² ≪ (Z²/Q)·η.
* **Result.** Altogether, in the normalisation of Step 2 (dual forms divided by Q):

  𝒜 ≪ (Z·Q^{θ′} + Z²/Q)·Z^ε.

  This is also the induction hypothesis at length Z. The base case (Z bounded) is trivial.

**Step 5 (repeated frequencies). PROVED.**
* Take k = k₀ϖ⁴ with ϖ a prime of Z[i], N(ϖ) ∈ [J/2, J], and J = (Z²/(QN(k₀)))^{1/4}. There are ≫ J/log J such ϖ.
* For these, G_k(0) ≥ 2^{−B}, and χ_c(k₀ϖ⁴) = χ_c(k₀)·1_{ϖ∤c}. So B₀(k₀ϖ⁴) = M − E_ϖ, where:
  * M := Σ_c μ(c)(k₀/c)₄W(N(c)/Z)/√N(c), over c in the chosen class mod 8. Since 1_{N(c)≡1 (8)} = ½(1+λ(c)), M is ½(M_{k₀,0} + M_{k₀,1})
    for the class 1 mod 8 and ½(M_{k₀,0} − M_{k₀,1}) for the class 5 mod 8, where M_{k₀,e} is the full sum twisted by ξ_{k₀,e}. Running the
    argument for both classes bounds each M_{k₀,e}.
  * E_ϖ collects the c divisible by ϖ. It has at most ηZ/N(ϖ) terms of size ≤ Z^{−1/2}, so |E_ϖ| ≪ Z^{1/2}/J trivially.
* If |M| ≥ 2|E_ϖ|, positivity of every term in 𝒜 gives (J/log J)|M|²/4 ≤ 𝒜. Otherwise |M| ≪ Z^{1/2}/J = Z^{x/4}, which is smaller than the
  bound below because x/4 ≤ ¼ < σ − ½.
* With Q = Z^x:

  |M|² ≪ (Z^{1+xθ′} + Z^{2−x})·Z^{−(2−x)/4+ε}.

* Choosing x = 1/(1+θ′) balances the two terms. Then 2 − x = (1+2θ′)/(1+θ′), so |M|² ≪ Z^{(3/4)(1+2θ′)/(1+θ′)+ε}.
* Unnormalised, this is Σ_c μ(c)ξ(c)W(N(c)/Z) = √Z·M ≪ Z^{½ + 3(1+2θ′)/(8(1+θ′)) + ε} = Z^{σ(θ′)+ε}.

**Step 6 (zeros).**
* The bound holds for every narrow-annulus weight W, uniformly. Smooth dyadic weights are finite sums of such W.
* So Σ_c μ(c)ξ(c)φ(N(c)/X) ≪ X^{σ+ε} for every smooth φ. Then φ̃(s)/L(s, ξ) is holomorphic in Re s > σ.
* For a putative zero ρ, choose φ with φ̃(ρ) ≠ 0. ∎

## Numerical check of the mechanism (l3j/j4pred.py, j4obs.py, j4var.py) — NUMERICAL

* **The control family.** Take β(c) = μ(c)T(c). On its Möbius side μ is replaced by μ², so P(k₀ϖ⁴) has size about Z. This is the
  extreme case of a "zero at s = 1", and the mechanism of Steps 2 and 5 should then inflate S/diag by about K^{1/4} = Z^{1/4} at Q = Z.
* **Setup.** Rows are all m ∈ Z[i] with V = e^{−πN(m)/Q}. The c are squarefree primary with split or inert prime factors and
  N(c) ∈ (Z/2, Z]. Q = Z.

| Z | observed S_{μT}/diag − 1 | predicted from the k = u·j⁴ frequencies alone | observed S_T/diag − 1 |
|---|---|---|---|
| 500 | 0.058 | 0.099 | 0.029 |
| 2000 | 0.159 | 0.173 | 0.004 |
| 4000 | 0.194 | 0.222 | 0.016 |

* **Prediction.** The predicted column evaluates the exact dual terms at k = u·j⁴ (u a unit), summed over coprime pairs.
* **Reading.**
  * The repeated frequencies account for 60–90% of the observed excess. The excess grows at the predicted K^{1/4} rate, and the other
    frequencies do not compensate it at M = 1.
  * For T itself, whose repeated frequencies carry genuine Möbius sums, there is no excess.
  * This matches Steps 2 and 5.
* **A clarification of older numerics.** The "untwisted control" of l3/gen_rows.c is g̃(χ̄_c) = ᾱT, not μT. Its Möbius side is μᾱ,
  which has no pole, which is why it showed no inflation (j4var.py reproduces both gen_rows ratios exactly: 1.05291 and 1.02713 at
  Z = Q = 4000).

## Why the sextic (OpenAI) analogue has no such barrier (HEURISTIC, structural)

* The quartic dual coefficient is T(c)g̃(k;χ̄_c) = μ(c)χ_c(k). This rests on g̃₄² being a Hecke character times the trivial quadratic
  Gauss sum, which is where the Möbius function enters.
* For the sextic family over Q(√−3), the Poisson dual coefficient is ᾱ(c)·g̃₃(c)·χ_c(k), a cubic Gauss sum (PROBLEM.md).
* There the repeated frequencies k₀ϖ⁶ carry Σ_c ᾱ(c)g̃₃(c)ξ_{k₀}(c). These are cubic Kubota sums, controlled by Patterson's cubic theta
  function, not by zeros of L-functions.
* So the sextic (G)-type lemma has no quasi-RH content and can be proved. The quartic (L³) cannot be proved without proving a quasi-RH
  for the quartic Hecke family.
* This is a structural reason why the quartic route is harder than the sextic one, beyond the self-duality of BELOW_3_2.

## Relation to OpenAI's quasi-RH argument (from scratchpad/qrh78_intro.tex, the saved introduction of their preprint)

* **Their mechanism is the same as Step 5.** Their base sum averages cubic-theta coefficients against sextic characters. Poisson in the
  averaging variable produces rows u·a⁶, and "for u = 1, this quotient contains the reciprocal of the target function". So their
  argument is the sextic analogue of Step 5, with sixth-power repetition giving 1/L.
* **Where their direct estimate comes from.** The direct estimate for their base sum comes from cubic-theta reflection (Patterson) and
  the quadratic large sieve.
* **Theorem ZF is the quartic analogue of their stage 1.** Here the direct estimate is (L³) itself, and the targets are the quartic
  Hecke characters over Q(i).
* **Matching numbers.**
  * The output 7/8 at θ′ = 0 matches PROBLEM.md's "k = 4: σ₀ = 7/8 − 3ℓ/8" at ℓ = 0.
  * OpenAI's stage 2 (prime compensation) is what would lower it further, e.g. to 13/16 or 5/6.
* **Consequence for the quartic route.**
  * In the sextic case the direct estimate is free, because cubic theta coefficients are explicit.
  * In the quartic case Poisson replaces the Gauss sums by μ, so no theta is left. The direct estimate, (L³), is then itself as strong
    as the quasi-RH it is meant to produce, at least for the quartic Hecke family over Q(i).
  * So the quartic route cannot reach any θ < 1 for Dirichlet L-functions without first proving a quasi-RH-strength statement for the
    quartic family.

## What is and is not new

* **NOTES §8G** ("Theorem Z") states PO(η,λ) ⟹ zero-free regions, and records "(L³) with loss Q^{θ′} ⟹ 7/8 + θ′/2" in outline.
  * At θ′ = ½ that outline gives 9/8, i.e. nothing.
* **Changes here.**
  1. The positivity of the dual form is supplied by the narrow-annulus expansion (Step 3). Without it the T-side bound does not give a
     positive k-sum: a log-convex weight cannot be flat and then decay (l3j/SESSION2.md §2).
  2. The row scale is optimised, Q = Z^{1/(1+θ′)} instead of Q = Z. This is what makes every θ′ < ½ nontrivial.
  3. Theorem Z's coprimality caveat is removed by taking ϖ prime: the error E_ϖ is trivially smaller than the bound.
* **The gcd > 1 pairs** are handled in Step 4. The twisted form of the hypothesis is used there.
* **The non-squarefree c of ψ_α** are handled in Step 1 by a finite local inversion, sketched there.
* **Still STANDARD, not written out:** the 2-adic normalisations (primary representatives, the units in χ_c(k₀)).
* An independent check should go through Steps 1–5 before any external use.
