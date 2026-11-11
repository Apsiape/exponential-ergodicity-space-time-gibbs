# Checks added in version 1.2.0

These programs check the algebraic and numerical ingredients of the material added in version 1.2.0:
Section 11.6 (the critical heat-bath PCA, Lemma 11.7 and Theorem 11.8) and Section 13 (the prefix-cap
theorem, Theorem 13.3). No proof in the paper uses a computation; these checks support the proofs and
are not part of them. They were written independently of the derivations they check, for an internal
review. The comments were edited for this release, and parts of `check_algebra.py` that concern material
outside this paper were removed; the remaining outputs are unchanged.

Requirements: Python 3 with NumPy (`check_cap_identity_and_messages.py`, `check_critical_mc.py`);
`check_algebra.py` needs only the standard library. The recorded outputs were produced with
Python 3.13 and NumPy 2.4 on Windows, running each program alone with `python -I -B`.

## Files

| File | What it checks | Time |
|---|---|---|
| `check_cap_identity_and_messages.py` | Two exact ingredients of Theorem 13.3, for binary PCA with random strictly positive tables and a random exterior configuration, on the minimal expanding volumes above one root. (1) The cap identity: the root marginal of the finite-volume specification, computed by brute force from the full product of transition factors over every site whose factor involves the volume, equals `p(a|parents) Z_a / sum_b p(b|parents) Z_b` computed by a forward transfer; and every exterior child of the volume lies in the cap row. (2) The message bound of Lemma 13.4: at every intermediate row, for every pair of rows at Hamming distance one, the forward expectation of the cap likelihood changes by a factor at most `K^{|N|}`; the program prints the largest observed log-ratio in units of `|N| log K` (it must be at most 1). Cases: `d=1`, `N={-1,0,1}`, 3 rows; `d=1`, 5 rows (transfer and message bound only, no brute force); `d=2`, `N={0,e1,e2}`, 3 rows; `d=2`, `N={0,±e1,±e2}`, 2 rows; three seeds each. | a few seconds |
| `check_algebra.py` | The inputs of Lemma 11.7 and Theorem 11.8: the nine table values and `sinh(2J_h) sinh(2J_v) = 1` in exact rationals; Cauchy's identity `det[1/(j-i+1/2)] = 2^N prod_k [k^2/(k^2-1/4)]^(N-k)` for `N <= 14` in exact rationals; the Fourier coefficients of the limiting symbol `i exp(-i theta/2)` by quadrature; the Toeplitz determinants of (I2) at `k = 0.5, 0.9, 0.999` for `N = 1, 2, 4, 8`, which approach `(1-k^2)^(1/4)` at `k = 0.5` and the critical product near `k = 1`; and the constants of Lemma 11.7(c), including `min_{n <= 2e5} (log C_n + (1/4) log n + 1) = 0.548 > 0` and `min_t (1/2)(2t+1)^(-1/8) / (t+1)^(-1/8) = 0.4585 >= 1/4`. | about a second |
| `check_critical_mc.py` | A Monte Carlo refutation attempt on the lower bound (47) of Theorem 11.8, with the expected outcome registered in the program's header before it was first run. Synchronous heat-bath dynamics on a 1024 x 1024 torus from all +1, two seeds, up to `t = 400`; and a control with both couplings multiplied by 0.9. | about a minute, under 1 GB |
| `*.out.txt` | Recorded outputs. `check_cap_identity_and_messages.out.txt` concatenates the four single-case runs `0`, `1`, `2`, `3`. | |
| `*.run.json` | Resource receipts written by the guard: process id, script name, free memory before the run, the memory cap, peak memory and wall time. | |
| `resource_guard.py` | Memory guard imported by each program. On Windows it caps the process at 1 GB, refuses to start below 9 GiB of free memory and writes `<script>.run.json` at exit; on other systems it does nothing. To run without it, delete the two lines of the program that import and call it. | |

## Run

From this directory, one program at a time:

    python -I -B check_cap_identity_and_messages.py 0     # cases 0..3; without an argument, all four
    python -I -B check_algebra.py
    python -I -B check_critical_mc.py

Each program prints its results to standard output. The programs are deterministic: the random tables,
exterior configurations and Monte Carlo runs use fixed NumPy seeds, so with the same NumPy version the
outputs should agree with the recorded ones line for line.

## Expected results

- `check_cap_identity_and_messages.py`: in all twelve instances, `exterior_children_not_in_cap` is 0, the
  brute-force and transfer root marginals agree to about `1e-16` where both are computed, and
  `lemma1_max_logratio_over_|N|logK` lies between 0.50 and 0.89; each case ends with `ALL CHECKS PASS`.
- `check_algebra.py`: the table values `1/37, 1/10, 1/5, 4/13, 1/2, 9/13, 4/5, 9/10, 36/37`; the product of
  the two sinh values equal to 1; `True` for Cauchy's identity; the Fourier coefficients equal to
  `1/(pi(1/2-n))` to six digits; the determinants at `k = 0.5` approaching `0.930605`; and the constants
  listed above.
- `check_critical_mc.py`: at the critical couplings the estimated `m_t` stays above `(1/2)(2t+1)^(-1/8)` at
  every printed time (`m_400 = 0.683` against `0.217`), with a log-log slope of about `-0.054` over
  `t = 50..400`; the control falls below the bound by `t = 200` (`m_200 = 0.166`, `m_400 = 0.037`).

## Conventions

In `check_algebra.py` the Fourier coefficient of the limiting symbol is taken with `exp(+i n theta)`, so
the printed `a_n` is the paper's `a_{-n}(1) = 1/(pi(1/2-n))`, and the matrix `[1/(j-i+1/2)]` is the
transpose of the paper's `[1/(i-j+1/2)]`; neither affects a determinant. The Toeplitz coefficients of the
symbol of (I2) for `k < 1` are taken with `exp(-i n theta)`, as in the paper.

In `check_cap_identity_and_messages.py` the expanding volume is the smallest one above a single root that
contains every child of each of its rows. The cap identity needs only that property; the volumes of
Section 13 are larger, with a buffer and a lateral speed at least the range of the neighbourhood, and
have it for the same reason.
