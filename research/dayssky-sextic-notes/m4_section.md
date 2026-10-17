
## 14. Attempt on the Möbius fourth moment (M4)

M4(r) means Σ_{q_u≤U} |M_u(U^r)|⁴ ≪ U^{1+ε}, with M_u(D) = D^{−1/2} Σ_n μ(n)ν(n)χ_n(u)^{±1} W(q_n/D).

### 14.1 What is actually proved

- **r ≤ 1/3, rigorous.** Write the second factor M_u(D) as a finite sum of products of prime slots of total length r (μ factorises over prime dyads). The paper's marked inverse moment (Lemma marked: r + 2z ≤ m, 2r + 8z ≤ 3m, with z = r, m = 1) then gives r ≤ 1/3.
- **Not proved for r ∈ (1/3, 1/2].**
  - Writing M_u² as one polynomial of length U^{2r} meets the BGL large-sieve term (U·U^{2r})^{2/3}, which exceeds U once r > 1/4.
  - Avoiding that term needs the canonical marked estimate (Lemma canonical-moment) for columns ᾱγ₂(c)·t(c), with t(c) = Σ_{c₁c₂=c} W W (c₁/c₂)₃. That is the A2 cubic double series. The paper's lemma explicitly excludes any residual column coefficient.
  - Even granted, the canonical regime requires dual rows ≤ columns: D⁴/U ≤ D², i.e. r ≤ 1/2.

### 14.2 Why r = 1/2 is a wall for the Poisson/theta toolkit (PROVED, exponent algebra)

- **Self-duality.** Poisson in u sends rows U and length D² to dual rows 𝓗 = D⁴/U and columns D². The two sides exchange roles exactly at D² = U.
- **Reflection lengthens past r = 1/2.**
  - Cubic theta reflection of a column of length C at row conductor q_k has length q_k²/C. In the dual this is 𝓗²/D² = D⁶/U², which exceeds D² iff r > 1/2.
  - The A2 functional equations reflecting c₁ at fixed c₂ have conductor q_k·q_{c₂}, so they lengthen too.
- **Bilinear (Type II) pieces.** For general coefficients, Σ_u |A_u(M)B_u(N)|² ≪ U||α||²||β||² is false once MN > U. Take α, β conjugate-aligned with ψ_{u₀}. So any argument must use the arithmetic of μ inside the bilinear pieces.
- **Classical analogue.** For all characters mod a single q, M4 at r = 3/4 is the Lindelöf-on-average sixth moment for an individual modulus, which is open. M4 at r slightly above 1/2 is the individual-modulus variance of μ*μ in residue classes at x = q^{1+δ}, which is a power-saving averaged-Chowla statement over shifts kq.

### 14.3 Range needed (Cycle61 framework, §13.4)
