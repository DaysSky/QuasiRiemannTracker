# S49. An unbalanced quartic large sieve, and (G)/(L³) at Q^{13/8}

Status labels: PROVED means a complete argument from published lemmas, cited precisely. Steps marked STANDARD are routine adaptations of a
published proof. They are spelled out below but have not been checked line by line against a referee-level write-up.

Notation (BGL14 §3, n = 4, K = Q(i), S = {(1+i)}).
* I(S): ideals coprime to 2, identified with primary generators.
* χ_b(a) = (a/b)₄, the quartic residue symbol.
* a ∼ M means M < N(a) ≤ 2M. Σ* means a sum over squarefree ideals.
* B₁^j(M,N) = sup_λ Σ*_{a∼M}|Σ*_{b∼N} λ_b χ_a^j(b)|² / ‖λ‖²: squarefree rows.
* B₂^j(M,N) = the same with all rows a ∈ I(S) and χ_b^j(a).
* B₃^j(M,N, C) = the smoothed form with coprime column pairs in a class C (BGL (3.1)).
* Ξ₄(A,B) as in dFDH (1.4). By BGL Lemma 3.1, Ξ₄(A,B) ≍ B₁¹(A,B) up to (AB)^ε and dyadic splitting.

Published inputs used:
* (BGL-T) B₁^j(M,N) ≪ (MN)^ε(M + N + (MN)^{2/3}) for j ∈ {1,3} [BGL14 Thm 1.3; j = 3 by complex conjugation].
* (GL) B₁²(M,N) ≪ (MN)^ε(M + N) [GL, Cor. 1.5; quadratic symbols].
* (L3.1) B₁^j(M,N) ≍ B₁^j(N,M) [BGL Lemma 3.1, duality and reciprocity].
* (3.4) B₁^j ≪ B₂^j [BGL (3.4)].
* (L3.2) B₂^j(M,N) ≪ N^ε B₃^j(M/Δ₁, N/Δ₂) for some 1 ≤ Δ₁, Δ₂ [BGL Lemma 3.2].
* (L3.3) B₃^j(M,N) ≪ M + (MN)^ε (M/N) max_{K ≤ (MN)^ε N²/M} B₂^j(K,N) [BGL Lemma 3.3]. This is the Poisson step.
* (CLS) Classical multiplicative large sieve over Z[i]: Σ_{N(q)≤R} Σ*_{ψ mod q} |Σ_{N(b)≤N} a_b ψ(b)|² ≪ (R² + N)‖a‖²
  [Huxley 1968, "The large sieve inequality for algebraic number fields"; here ψ runs over primitive characters of (Z[i]/q)^*].

## Lemma A (sharpened BGL Lemma 3.4 for n = 4) — PROVED
For j ∈ {1,3} and all K, N ≥ 1/2:   **B₂^j(K,N) ≪ (KN)^ε (K + K^{1/4}N + (KN)^{2/3}).**

Proof.
1. **Decomposition.** Write a = a₁a₂²a₃³a₄⁴ with a₁, a₂, a₃ squarefree and pairwise coprime, and a₄ ∈ I(S); this representation is unique.
   Split into (log 2K)⁴ dyadic boxes X = (X₁,X₂,X₃,X₄) with a_i ∼ X_i, so that X₁X₂²X₃³X₄⁴ ≍ K. Put P = X₁X₂X₃.
2. **Four bounds per box.** For λ with ‖λ‖ = 1, the box contribution Σ_{a∈X}|Σ_b λ_b χ_b^j(a)|² is bounded by each of the following.
   * (c_ℓ), ℓ ∈ {1,2,3}. Fix the a_i with i ≠ ℓ and put λ*_b = λ_b χ_b^j(a/a_ℓ^ℓ), so |λ*| ≤ |λ|. The inner sum is
     Σ*_{a_ℓ∼X_ℓ}|Σ_b λ*_b χ_b^{jℓ}(a_ℓ)|². By reciprocity, after splitting b into its two classes N(b) mod 8 as in BGL (3.4), this is
     ≪ B₁^{jℓ}(X_ℓ,N). Summing over the fixed variables gives (X/X_ℓ)·B₁^{jℓ}(X_ℓ,N), where X = X₁X₂X₃X₄ (BGL (3.6)).
     Here jℓ ≢ 0 (mod 4) for ℓ = 1,2,3 and j odd. BGL use only ℓ ∈ {1,2}; ℓ = 3 is equally valid for n = 4.
     For ℓ = 2 the character is quadratic, so (GL) applies.
   * (c4). Fix a₄. Then χ_b(a₄⁴) = 1_{(b,a₄)=1}, which can be absorbed into λ. By quartic reciprocity for primary elements,
     χ_b(a₁a₂²a₃³) = ±ψ_{a₁,a₂,a₃}(b), where ψ_{a₁,a₂,a₃}(b) = (b/a₁)₄(b/a₂)₄²(b/a₃)₄³. The sign is constant once b is fixed in a class N(b) mod 8.
     ψ_{a₁,a₂,a₃} is primitive modulo a₁a₂a₃: each local factor (·/p)^e with e ∈ {1,2,3} is nontrivial on (Z[i]/p)^*, whose order N(p)−1 is
     divisible by 4. Distinct triples give distinct characters, because the conductor determines a₁a₂a₃ and the local order of the character at
     p determines which a_i contains p. For j = 3 conjugate. By (CLS) with R ≍ P, the inner sum is ≪ (P² + N).
     Summing over a₄ gives X₄(P² + N).
3. **Case analysis.** Fix a box.
   * (i) If X₁ ≥ N^{1/2}, then (BGL-T) gives (c1) ≪ X₂X₃X₄(X₁ + (X₁N)^{2/3}).
     First, X₁X₂X₃X₄ ≪ K. Second, X₂X₃X₄ ≤ (X₂²X₃³X₄⁴)^{2/3} up to constants, because the exponents satisfy 1 ≤ 4/3, 1 ≤ 2 and 1 ≤ 8/3.
     Hence (c1) ≪ K + (KN)^{2/3}.
   * (ii) If X₁ < N^{1/2} ≤ X₃, then (c3) ≪ X₁X₂X₄(X₃ + N + (X₃N)^{2/3}) ≪ K + X₁X₂X₄(X₃N)^{2/3}.
     Next, X₁X₂X₄X₃^{2/3} ≤ K^{2/3} ⟺ X₁^{1/3} ≤ X₂^{1/3}X₃^{4/3}X₄^{5/3}, which holds since X₁ < X₃.
   * (iii) If X₂ ≥ N, then (c2) ≪ X₁X₃X₄(X₂ + N) ≪ K.
   * (iv) Otherwise X₁, X₃ < N^{1/2} and X₂ < N, so B₁^{jℓ}(X_ℓ,N) ≪ N for ℓ = 1, 2, 3.
     Let m = max(X₁,X₂,X₃) and σ = P/m. Choosing ℓ as the index attaining m gives min(c1,c2,c3) ≪ X₄Nσ.
     If P² ≤ N, then (c4) ≪ X₄N ≪ K^{1/4}N.
     Otherwise min(X₄Nσ, X₄P²) ≤ (X₄Nσ)^{2/3}(X₄m²σ²)^{1/3} = X₄N^{2/3}(mσ²)^{2/3}.
     By the rearrangement inequality, X₁X₂²X₃³ ≥ m·(mid)²·(min)³ ≥ mσ². Hence X₄(mσ²)^{2/3} ≤ X₄(K/X₄⁴)^{2/3} ≤ K^{2/3}.
     The box is therefore ≪ (KN)^{2/3}.
4. Summing over the (log)⁴ boxes proves the lemma. ∎

Remarks.
* The bound is sharp in the K^{1/4}N term. The rows a = a₄⁴ with λ ≡ 1 give B₂ ≥ K^{1/4}N.
* BGL's Lemma 3.4 gives K^{1/3}N in place of K^{1/4}N + (KN)^{2/3}. The improvement comes from (c4), the classical sieve on mixed boxes, and
  from (c3). Without (c4) the case (iv) analysis returns exactly K^{1/3}N.

## Theorem 1 (quartic large sieve, unbalanced ranges) — PROVED
For all A, B ≥ 1:

  **Ξ₄(A,B) ≪_ε (AB)^ε (A + B + A^{3/4}B^{1/2} + A^{1/2}B^{3/4} + min(A^{1/3}B, AB^{1/3})).**

In particular Ξ₄(A,B) ≪ (AB)^ε(A + B + A^{3/4}B^{1/2} + A^{1/2}B^{3/4}) whenever max(A,B) ≥ min(A,B)^{6/5}. Thus dFDH Conjecture 1.3 holds in
that range, and Theorem 1 is sharp there by dFDH Theorem 1.2. For A ≠ B it is a power saving over the (AB)^{2/3} of BGL14 Thm 1.3. At A ≍ B it
coincides with BGL.

Proof.
1. By (L3.3) and Lemma A, for any M, N:
   B₃(M,N) ≪ M + (MN)^ε (M/N)·[K + K^{1/4}N + (KN)^{2/3}] evaluated at K = (MN)^ε N²/M. The bracket is increasing in K.
   This gives (MN)^{2ε}(M + N + M^{3/4}N^{1/2} + M^{1/3}N). If N²/M < 1/2, the max in (L3.3) is empty and B₃ ≪ M.
2. The right-hand side is increasing in M and N. So (3.4) and (L3.2) give
   B₁(M,N) ≪ B₂(M,N) ≪ N^ε max_Δ B₃(M/Δ₁,N/Δ₂) ≪ (MN)^{3ε}(M + N + M^{3/4}N^{1/2} + M^{1/3}N).
3. Apply (L3.1) to swap M and N, take the minimum of the two bounds, and use Ξ₄ ≍ B₁. ∎

Exponent check (r = log A / log B ≥ 1):
* BGL: 2(1+r)/3.
* Theorem 1: max(r, 1 + r/3, 1/2 + 3r/4).
* dFDH lower bound: max(r, 1/2 + 3r/4).

For example, r = 1.5 gives 5/3, 1.625 and 1.625; r = 1.2 gives 1.467, 1.4 and 1.4.

## Lemma B (Poisson in the column variable, box-structured columns) — PROVED for coprime-radical pairs; the remaining pairs are STANDARD (BGL Lemma 3.2 pattern)
Let the box 𝔅 consist of m = m₁m₂²m₃³ with m_i ∈ I(S) squarefree, pairwise coprime and m_i ∼ X_i. Put P = X₁X₂X₃ and n = (m₂, m₃).
For Y ≥ 1 and any coefficients γ:

  Σ*_{c∼Y} |Σ_{m∈𝔅} γ_m χ_c(m)^{-1}|² ≪ (PY)^ε ( Y + (Y/P)·X₂X₃·max_{K≤(PY)^εP²/Y} B₂(K,X₁) ) ‖γ‖².

Proof.
1. **Smoothing and reciprocity.** By positivity, replace the sum over squarefree c ∼ Y by Σ_{c∈I(S)} W(N(c)/Y) with a smooth majorant W ≥ 0.
   Here (m/c)₄ is the Jacobi-type symbol. By quartic reciprocity for coprime primary elements, (m/c)₄ = ε·ψ_m(c) with
   ψ_m(c) = (c/m₁)(c/m₂)²(c/m₃)³. The sign ε is constant once N(m) and N(c) are fixed mod 8, so split into ≤ 4 classes by Cauchy–Schwarz.
2. **Opening the square.** Σ_{m,m'} γ_m γ̄_{m'} Σ_c W(N(c)/Y) 1_{c∈κ} ψ̄_m(c) ψ_{m'}(c).
   * Pairs with rad m = rad m' and equal labels at every prime are exactly m = m'. They give the diagonal ≪ Y‖γ‖².
   * Pairs whose radicals are coprime: χ := ψ̄_m ψ_{m'} is primitive mod qq' (q = m₁m₂m₃) and nontrivial. Poisson summation over the class κ
     mod 16 and mod qq' gives (Y/N(qq'))·Σ_{h≠0} Ŵ(h·(Y)^{1/2}/(16qq'))·e_κ(h)·χ̄(h)τ(χ).
     Ŵ is Schwartz, so the h-sum may be truncated at N(h) ≤ (PY)^ε N(qq')/Y ≤ (PY)^ε P²/Y, up to O((PY)^{-A}).
3. **Factoring τ(χ).** By CRT, τ(χ) = ψ̄_m(q')ψ_{m'}(q)·τ(ψ̄_m)τ(ψ_{m'}), and |τ(ψ_m)| = N(q)^{1/2}.
   Expanding ψ̄_m(q')ψ_{m'}(q) = Π_{i,j}(m_j'/m_i)^{-i}(m_i/m_j')^{j}:
   * the (1,1) factor is a sign by reciprocity;
   * factors with i = 1 < j depend only on (m₁, n');
   * factors with j = 1 < i depend only on (m₁', n);
   * factors with i, j ≥ 2 depend only on (n, n').
   Likewise τ(ψ_m) is N(q)^{1/2} times a unimodular function of m, and χ̄(h) = ψ̄_m(h)ψ_{m'}(h).
4. **Separating variables.** Separate N(q)N(q') in Ŵ by Mellin inversion, with |Im s| ≤ (PY)^ε (BGL p. 9). Remove the coprimality conditions
   (m₁,m₁') = 1, (m₁,n') = 1 and (n,m₁') = 1 by Möbius inversion; this keeps the expression separable, as in [HB2, p.121] and BGL Lemma 3.3.
   For each fixed (n, n') the off-diagonal is then
   ≪ (Y/P)(PY)^ε Σ_{h≤H} |Σ_{m₁} a_{m₁} ψ_{m₁}(h)|·|Σ_{m₁'} b_{m₁'} ψ_{m₁'}(h)|, with |a_{m₁}| ≤ |γ_{m₁,n}| and |b_{m₁'}| ≤ |γ_{m₁',n'}|.
   The factors ψ_n(h) and ψ_{n'}(h) are unimodular in h and drop out under absolute values.
   Cauchy–Schwarz in h bounds this by B₂(H,X₁)‖γ_{·,n}‖‖γ_{·,n'}‖. ψ_{m₁}(h) = (h/m₁)₄ is χ_{m₁}(h) in BGL notation, and the dual h runs over a
   lattice commensurable with O. A unit or λ² rescaling of h multiplies by a character of m₁ that is constant on classes, so it is absorbed.
   Summing over (n, n') gives (Σ_n ‖γ_{·,n}‖)² ≤ X₂X₃‖γ‖².
5. **Pairs sharing primes (STANDARD).** Sort by d = gcd(rad m, rad m') and by the labels at p | d, as in BGL Lemma 3.2.
   * Primes with equal labels give 1_{(c,p)=1}. Möbius over r | d shortens c to Y/N(r).
   * Primes with different labels i ≠ i' contribute a fixed primitive character (·/p)^{i'−i} of conductor d_≠. It multiplies the Poisson
     modulus by N(d_≠), contributes |τ| = N(d_≠)^{1/2}, and is common to both sides of the Cauchy–Schwarz.
   * With X_i' = X_i/N(d_i) and P' = P/N(d), the prefactor (Y/P')X₂'X₃' = Y/X₁' and the dual length H' ≤ (PY)^ε P²N(d_≠)/(Y N(d)²) ≤ H/N(d).
     Each term of Y/X₁'·[H' + H'^{1/4}X₁' + (H'X₁')^{2/3}] is then at most the d = 1 value times N(d)^{-1/3}. The d-sum costs (PY)^ε. ∎

## Theorem 2 ((G) at M = 2/3 to Q^{13/8}; (L³) to Q^{13/8}) — PROVED, given Lemma B step 5 (STANDARD) and Kubota convexity (input K below)
Let T(c) = α(c) g̃(χ̄_c) on squarefree c (0 otherwise) and R_Y(m) = Σ_c T(c)χ̄_c(m)W(N(c)/Y). For 1 ≤ Y ≤ Q^{5/3}:

  **Σ_{N(m)≤Q} |R_Y(m)|² ≪ (QY)^ε (Y + QY^{1/3} + Q^{1/2}Y^{3/4})·Y.**

At Y = Z = Q^{3/2} this is Q^{13/8}·Z^{1+ε}, versus the target QZ^{1+ε}. Via the approximate functional equation of ψ_α (both pieces have
length ≤ N(m³)^{1/2+ε} = Q^{3/2+ε}), this gives

  Σ_{N(m)≤Q} |ψ_α(m³, ½+it)|² ≪ Q^{13/8+ε}(1+|t|)^{O(1)}.

**Input K (Kubota convexity).** ψ_α(ν,s) = Σ_c T(c)χ_c(ν)N(c)^{-s} (normalised, centre ½), including its twists by characters mod 16, is entire.
It satisfies |ψ_α(ν,½+it)| ≪ N(ν)^{1/4+ε}(1+|t|)^{O(1)}. This is entire because the α-twist means ℓ ≠ 0, and it follows from the functional
equation DDHL (4.28) together with Phragmén–Lindelöf; this is the same input as S43.

Fine print (STANDARD). Our series runs over squarefree c. Squarefree c with (c,ν) ≠ 1 contribute 0, since g(χ̄_p) summed against e(νx/p)
vanishes for p | ν. Non-squarefree c with (c,ν) = 1 also contribute 0, since g(ν,p^k) = 0 for k ≥ 2 and p ∤ ν. So the series equals the full
Kubota series restricted to (c, rad ν) = 1. That differs from the full series by a finite local factor at the primes p | ν. This is the
analogue of dFDH (9.6) and [CdFD26, Lemma 6.1(iii)], stated there for squarefree ν. On Re s = ½ + ε the factor is ≍ N(ν)^{±ε}: each local
term is 1 + O(N(p)^{−1/2}) with N(p) ≥ 5. The non-squarefree-ν version of that lemma is assumed here in the same form; it is a local computation.

**Corollary for rows (PROVED from K).** For m = m₁m₂²m₃³m₄⁴, put ν₀ = m₁³m₂²m₃ and x = m₂m₃²m₄³. Then m³ = ν₀x⁴ and χ_c(x⁴) = 1_{(c,x)=1}.
Möbius over d | rad x, together with the twisted multiplicativity T(dc') = ±(c'/d)₂ T(d)T(c') = ±χ_{c'}(d²)T(d)T(c') (CRT for Gauss sums plus
quartic reciprocity; the sign is constant on classes mod 16), gives

  |R_Y(m)| ≤ Σ_{d|rad x} |Σ_{c'} ±T(c')χ_{c'}(ν₀d²)W(N(d)N(c')/Y)| ≪ Σ_d (Y/N(d))^{1/2}(N(ν₀)N(d)²)^{1/4}Y^ε ≪ Y^{1/2+2ε}N(ν₀)^{1/4}.

N(d) cancels between the shorter length and the larger index. The condition (c',d) = 1 is automatic because χ_{c'}(d²) = 0 otherwise.

**Proof of Theorem 2.** Work in Q-units: Y = Q^η with η ≤ 5/3, and a box is (x₁,x₂,x₃,y) with x₁ + 2x₂ + 3x₃ + 4y ≤ 1 and p = x₁ + x₂ + x₃.
The target is the exponent τ(η) = 1 + max(η − 1, η/3, 3η/4 − 1/2) in units of Q, relative to Y. Box contributions are normalised by
‖T‖² ≍ Y. Let Δ be the normalised box bound.
* **Case I (P² ≤ Y).**
  * Classical large sieve on rows with m₄ fixed: Δ ≤ X₄(Y + P²) ≤ 2X₄Y.
  * Input K summed over the box: Δ ≤ X₁X₂X₃X₄·N(ν₀)^{1/2}, i.e. exponent y + p + (3x₁+2x₂+x₃)/2 ≤ y + (5/2)(x₁+2x₂+3x₃) ≤ 5/2 − 9y.
  * The first bound wins if y ≤ τ − η; otherwise the second gives < 5/2 − 9(τ − η) ≤ τ for η ≤ 5/3 (checked at the binding point).
* **Case II (P² > Y).** Apply the transpose of the box matrix (rows c, columns m). First Cauchy over m₄ (factor X₄), then drop
  1_{(c,m₄)=1} by positivity, then apply Lemma B with Lemma A for B₂(H, X₁), where H = P²/Y, i.e. h = 2p − η.
  The exponent terms are:
  * (a) y + η ≤ (1−p)/4 + η < 1/4 + 7η/8 ≤ τ.
  * (b) y + η − p + x₂ + x₃ + h = y + p + x₂ + x₃ ≤ 1.
  * (c) y + η + h/4 = y + 3η/4 + p/2 ≤ 1/4 + 3η/4 + p/4 ≤ 1/2 + 3η/4.
  * (d) y + η/3 + x₁ + (4/3)(x₂+x₃) ≤ 1 + η/3.
  All four are ≤ τ(η). ∎
* Numerical grid confirmation of the box bookkeeping: l3c/boxes_L3b.py (Y = Z) and l3c/boxes_Y.py (η ∈ [0.25, 1.5]). Both have worst box
  exactly the main box, with value 1/2 + 3η/4.

## What this does and does not give
* It improves the unconditional state of the art for Σ_{N(m)≤Q}|ψ_α(m³,½)|² from Q^{5/3} (BGL; Q^{7/4} by our earlier recursion) to Q^{13/8}.
  That is the limit of any large-sieve argument: dFDH Thm 1.2 shows Ξ₄(Q, Q^{3/2}) ≫ Q^{13/8−ε}.
* (L³) needs Q^{1+ε}. The remaining gap Q^{5/8} is the T-specific "Lindelöf on average for Kubota series" problem (S46), and it is unchanged.
* θ: no change. §7 of L3_REPORT: any θ gain needs (G) near M = 1/2 beyond the optimal large-sieve shape, which no large sieve can give.
* Theorem 1 is a self-contained result about the quartic large sieve. As far as the literature in lit/ shows, it is new: dFDH (2026) state
  Conjecture 1.3 for all A, B and record no range where it is known. It uses only published lemmas plus Lemma A.

## Remark: general n (same proof; status STANDARD, since BGL's general-n lemmas are used with ray-class bookkeeping)
For n ≥ 4 over K ⊇ Q(ζ_n) (BGL §3 setting), decompose a = Π_{i=1}^{n} a_i^i with a₁, …, a_{n−1} squarefree and pairwise coprime.
* Sieving in a_ℓ uses χ^ℓ, of order n/gcd(ℓ,n) ≥ 2. BGL (E_{4/3}) applies when that order is ≥ 3, and GL when it is 2.
* The classical box fixes a_n and gives X_n(N + P²) with P = Π_{i<n} X_i: distinct tuples give distinct primitive characters of conductor
  Π a_i, because the local order at p determines the label.
* Case (iv) uses rearrangement: Σ i·x_i ≥ x_max + 2·Σ_{others} x_i, so mσ² ≤ Π X_i^i.
Together these give Lemma A_n: B₂(K,N) ≪ (KN)^ε(K + K^{1/n}N + (KN)^{2/3}). With BGL Lemmas 3.1–3.3 this yields

  Ξ_n(A,B) ≪ (AB)^ε (A + B + A^{1−1/n}B^{2/n} + A^{2/n}B^{1−1/n} + min(A^{1/3}B, AB^{1/3})).

This is the dFDH prediction (1.8) whenever max(A,B) ≥ min(A,B)^{(3n−6)/(2n−3)} (6/5 for n = 4, 4/3 for n = 6, → 3/2 as n → ∞).
For n = 3 the same argument returns Heath-Brown's A + B + (AB)^{2/3}, which is sharp by dFDH Thm 1.1. That is a consistency check, since the
cubic j-box term A^{2/3}B^{2/3} coincides with the BGL term.
