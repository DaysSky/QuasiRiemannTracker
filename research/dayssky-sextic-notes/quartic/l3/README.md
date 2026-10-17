# l3/: numerics for the quartic Lindelöf-on-average problem (L³)

All programs use the conventions of `../PROBLEM.md`:

* c is a squarefree primary element of Z[i]; χ_c = (·/c)₄; g̃ is the normalized Gauss sum with ψ(z) = e(Im z).
* T(c) = α(c) g̃(χ̄_c), with α(c) = c/|c|.
* W(x) = exp(−1/((x−½)(1−x)) + 16) on (½, 1).

Build each program with `gcc -O3 -march=native -fopenmp -o X X.c -lm`. `helpers.inc` holds the Z[i] arithmetic, extracted from `../barrier/gen_big.c`.

| Program | Purpose | Output |
|---|---|---|
| `gen_rows.c` | rows R(m) = Σ_c T(c)χ̄_c(m)W(N(c)/Z) for all primary N(m) ≤ Q, plus the untwisted control (T → g̃) | JSON summary line; binary per-row (N, \|R_tw\|², \|R_untw\|², rowdiag) |
| `gen_multi.c` | same, over J disjoint dyadic Z-windows in one pass, with Q_j = Z_j^M (args: Z Q out M J) | one JSON line per window |
| `gen_k.c` | index-m rows: Möbius-weighted (Σ μ g̃(χ̄_c)χ_c(m)W) and parity-blind (Σ g̃(χ̄_c)χ_c(m)W) | as gen_rows |
| `gen_T.c` | dump (N, T) and the exponent bytes e_k = log_i (m_k/c)₄ for every squarefree c ≤ Zmax | binary |
| `gen_T2.c` | as gen_T, also dumping c = a + bi | binary |
| `primes_g.c` | g(χ̄_π) at all split primes π with N(π) ≤ X | binary (a, b, p, Re G, Im G) |
| `gen_rows_prime.c` | gen_rows restricted to prime columns c = π (S32) | as gen_rows |
| `gen_rows_sharp.c` | gen_rows with the sharp weight W = 1 on (½, 1] (S30) | as gen_rows |
| `gen_rows_c.c` | gen_rows that also writes complex R_tw(m) (phase test) | binary with 6 doubles per row |
| `gen_corr.c` | fixed Q, J disjoint dyadic c-windows; mode 0 = T-family, mode 1 = d₃ Euler-product control (all primary c) (S41); args `Z Q out mode J` | binary: NR, J, ±N(m) (negative = not squarefree), then (\|R\|², rowdiag) per window and row |

| Script | Purpose |
|---|---|
| `rowstats.py` | per-row distribution against Exp(1), bands, top rows |
| `corrstats.py` | row moments and cross-window correlations from gen_corr output (S41) |
| `sigtest.py`, `sigtest2.py` | prime-sign and Hecke-twist correlation tests |
| `rz.py` | small-index rows across many Z, from a gen_T dump |
| `tasp.py` | t-aspect test (rows twisted by N(c)^{−it}) |
| `theta_fit.py` | Z^{3/4} fit of parity-blind index-m sums, from a gen_T2 dump |
| `theta_fit2.py`, `theta_fit3.py`, `theta_rule.py`, `theta_cls.py`, `theta_units.py`, `theta_cube.py`, `theta_phase*.py` | quartic-theta residue pattern at small indices (S29, S33) |

Validation: T(c) and the symbol exponents match brute-force `../gauss4.py` for all c with N(c) ≤ 600 (S17 in `../L3_ATTEMPT.md`).

Results: `scale1.txt`, `scale3e6.txt`, `scale6e6.txt`, `mscan3e5.txt`, `mscan_M1.txt`, `multi1e6.txt`, `multi25e5.txt`, `sharp.txt`, `prime1e6.txt`. Analysis is in `../L3_ATTEMPT.md` (S1–S41) and `../L3_REPORT.md`.

`../l3b/mds_group.py`: the functional-equation group of the second-moment multiple Dirichlet series (S39).
