I did not prove (L³), and I did not disprove it either.

What I produced instead:
1. A conditional proof: a ratios-type hypothesis (R) for the quartic Hecke family implies (G) with a power saving. Weakest sufficient form: error H^{3/4 − 3σ₀/2} at Re s = ½ + σ₀. Its algebra was checked by brute force.
2. A sharper picture of the problem:
   - Rows are Möbius-weighted Kubota sums.
   - It is a parity-sensitive statement.
   - It is the GL3-type thin-family second moment, which is open even in the classical setting.
   - Every Poisson, reflection or large-sieve route stops at the self-dual point M = 2/3.
3. Numerics, validated against brute force: S/diag = 1.008 ± 0.007 over 36 windows up to Z = 2.5e6; Z = 6e6: [ ].
4. A new obstruction for the quartic route: with prime-only columns the mean square is biased and growing (1.23 → 1.29 → 1.44 for Z = 1e5 → 1e6), so quartic prime-slot moments need the quartic theta coefficients explicitly.
5. Side observation: the parity-blind theta terms follow a 2-adic rule (zero at π ≡ 3+2i mod 4, at e = 3, and so on).

θ / PR: upstream merged #53 (not lower). I rebased onto 24fc1ccd and re-ran all the checks, which pass; the commit is still DaysSky's and not pushed. Levers analysed.
