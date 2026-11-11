# Verification: Exponential ergodicity and the space–time Gibbs states of one-dimensional PCA (version 1.2.0)

The computational content of the paper is Table 2 (Section 11.3, noisy majority certificates) with the
three polynomials q_1, q_2, q_3 quoted there, and, since version 1.2.0, the checks described in
`v1.2/VERIFY.md` for Sections 11.6 and 13. No proof uses a computation. This directory reproduces
every number in Table 2 in exact rational arithmetic, with two independent programs, and checks the
recorded outputs by hash.

## Added in version 1.2.0
- v1.2/: checks of the critical example (Section 11.6) and of the prefix-cap theorem (Section 13):
  three programs, their recorded outputs and resource receipts. Start with `v1.2/VERIFY.md`.

## Contents
- code/majority_certificates.py           the program that produced Table 2 (Python 3; needs sympy and numpy)
- code/majority_certificates_output.txt   recorded output of `python3 majority_certificates.py 7`
- code/majority_certificates_output.json  the same numbers, machine-readable
- code/SHA256SUMS                         digests of the three files above
- independent/table2_check.py             a second, independent program that re-derives Table 2 (seconds)
- VERIFY.md                               this file

The repository's root `SHA256SUMS` also covers every file here.

## Check the recorded files (seconds)
    cd code && sha256sum -c SHA256SUMS
Expected: three lines ending in "OK".

## Reproduce (depth 5: ~5 s; depth 6: ~30 s; depth 7: ~12–20 min, all in exact rationals)
The script writes `majority_certificates_output.txt` and `majority_certificates_output.json` into the
directory it is run from, replacing files of those names, so run it from a fresh copy:

    pip install sympy numpy
    mkdir fresh && cp code/majority_certificates.py fresh/ && cd fresh
    python3 majority_certificates.py 7
    diff majority_certificates_output.txt ../code/majority_certificates_output.txt

Apart from the timing annotations "[... s]" the two files should agree. A fresh run to depth 5
reproduced every number through L = 5 exactly.

## What the script does
Part 1 (L ≤ 3): enumerates the 3^{L^2} innovation configurations of the depth-L cone (each update
site is Deterministic / Noisy-0 / Noisy-1), records integer counts by number of noisy updates,
forms the exact polynomials q_L(ε), evaluates them at rational points, and finds the maximizing
input pair for θ_L by exhaustive comparison over all 2^{2L+1} × 2^{2L+1} pairs of input rows.
Part 2 (L ≤ 7): computes q_L = P(F_L(1)(0)=1) − P(F_L(0)(0)=1) by an exact row transfer,
cross-checks against Part 1 at ε = 2/3 for L ≤ 3 (they agree), locates the first passing grid point
k/1000 in floating point, and certifies it in exact rational arithmetic (fractions.Fraction).

## What it certifies (Table 2)
(2L+1) q_L(ε) < 1, hence all-DLR uniqueness and exponential ergodicity (Theorem 3.1 with
Propositions 4.1–4.2 and 4.4), holds exactly at
    L=1: ε=667/1000   L=2: 627/1000   L=3: 587/1000   L=4: 277/500
    L=5: 527/1000     L=6: 63/125     L=7: 97/200   (so at noise 1/2 by depth 7)
and (2L+1) q_L at ε = 3/5 and ε = 1/2 are the exact rationals listed in the output file.
For the monotone majority rule θ_L = q_L (Proposition 4.4); Part 1 confirms this for L ≤ 3.
Every noise above 1/2 is also covered by the analytic high-noise bound of Lebowitz, Maes and Speer
recalled in Section 11.3; the depth-seven certificates at 1/2 and 97/200 go beyond that range.

## Independent second program
`independent/table2_check.py` shares no code with `majority_certificates.py` and uses a different algorithm: an
exact integer transfer over the backward light cone, one update site at a time, for each rational noise. It
re-derives every entry of Table 2 at all seven depths, checks that each certified noise is the first passing grid
point k/1000, compares the results with the paper's table and with the recorded outputs, and checks the Part 1
polynomials and the analytic high-noise bound at sample noises. It needs SymPy and runs in a few seconds:

    cd independent && python -B table2_check.py

Every check passes on the shipped files.

## Not computational
Everything else in the paper is proved in the text; the checks in `v1.2/` support proofs and are not
part of them. The defect bound of Section 10 (Theorem 10.3)
was also tested numerically on every subset of small space–time windows, as "How this work was
done" in the paper describes; that test is a sanity check, not a certificate, and is not part of
this package.
