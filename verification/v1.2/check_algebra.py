"""Checks of the algebraic inputs of the critical example (Section 11.6 of the
paper: Lemma 11.7 and Theorem 11.8), written independently of the derivation.

1. The heat-bath table at J_h=(log 3)/2, J_v=(log 2)/2: its nine values, and
   sinh(2J_h) sinh(2J_v) = 1, in exact rational arithmetic.
2. Cauchy's determinant: det[1/(j-i+1/2)] = 2^N prod_k [k^2/(k^2-1/4)]^(N-k)
   for N = 1..14 in exact rationals; the Fourier coefficients of the limiting
   symbol i exp(-i theta/2) by quadrature; the low-temperature Toeplitz
   determinants of (I2) at k = 0.5, 0.9, 0.999, compared with (1-k^2)^(1/4)
   (the square of the magnetization in (I1)) and with the critical product.
3. The constants of Lemma 11.7(c): log C_n >= -(1/4) log n - 1 for n <= 2e5,
   the elementary inequality log(1-x) >= -x - (2/3)x^2 on [0,1/4], and
   (1/2)(2t+1)^(-1/8) >= (1/4)(t+1)^(-1/8).

Conventions.  In part 2 the Fourier coefficient is taken with exp(+i n theta),
so the printed a_n equals the paper's a_{-n}(1) = 1/(pi(1/2-n)), and the Cauchy
matrix [1/(j-i+1/2)] is the transpose of the paper's [1/(i-j+1/2)]; neither
change affects a determinant.  The Toeplitz coefficients of the k < 1 symbol
are taken with exp(-i n theta), as in (I2).

Run: python -I -B check_algebra.py   (single process, about a second)
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from resource_guard import guard
guard()

from fractions import Fraction as Fr
from math import log, exp, pi, cos, sin

out = []
def say(*a):
    s = ' '.join(str(x) for x in a)
    print(s, flush=True); out.append(s)

# ---------- 1. Heat-bath table ----------
say('== 1. critical heat-bath table')
vals = []
for H in (-2, 0, 2):
    for V in (-2, 0, 2):
        x = Fr(3)**H * Fr(2)**V
        vals.append(x / (1 + x))
ps = sorted(set(vals))
say('p values', [str(p) for p in ps])
sh = (Fr(3) - Fr(1, 3)) / 2 * (Fr(2) - Fr(1, 2)) / 2
say('sinh(2Jh)sinh(2Jv) =', sh)

# ---------- 2. Cauchy determinant ----------
say('== 2. Cauchy determinant')
def det(M):
    M = [row[:] for row in M]; n = len(M); d = Fr(1)
    for c in range(n):
        piv = next((r for r in range(c, n) if M[r][c] != 0), None)
        if piv is None: return Fr(0)
        if piv != c: M[c], M[piv] = M[piv], M[c]; d = -d
        d *= M[c][c]
        for r in range(c + 1, n):
            f = M[r][c] / M[c][c]
            for k in range(c, n): M[r][k] -= f * M[c][k]
    return d
ok = True
for N in range(1, 15):
    M = [[1 / (Fr(j - i) + Fr(1, 2)) for j in range(N)] for i in range(N)]
    lhs = det(M)
    rhs = Fr(2) ** N
    for k in range(1, N):
        rhs *= (Fr(k * k) / (Fr(k * k) - Fr(1, 4))) ** (N - k)
    ok &= (lhs == rhs)
say('det[1/(j-i+1/2)] == 2^N prod[k^2/(k^2-1/4)]^(N-k) for N=1..14:', ok)
# Fourier coefficient of i*exp(-i theta/2) on (0,2pi): midpoint rule
def fourier(n, L=200000):
    re = im = 0.0
    for s in range(L):
        th = (s + 0.5) * 2 * pi / L
        # exp(i n th) * i * exp(-i th/2) = i*exp(i (n-1/2) th)
        ang = (n - 0.5) * th
        re += -sin(ang); im += cos(ang)
    return re / L, im / L
for n in (-2, -1, 0, 1, 2, 3):
    re, im = fourier(n, 20000)
    say('a_%d numeric (%.6f, %.6f) vs 1/(pi(1/2-n)) = %.6f' % (n, re, im, 1 / (pi * (0.5 - n))))

# low-temperature symbol near k=1: Toeplitz det vs critical product, and Szego check
def symbol_coeffs(k, nmax, L=40000):
    # phi(theta) = [(1-k e^{-i th})/(1-k e^{i th})]^{1/2}, the continuous branch
    # equal to 1 at th=0, which is conj(z)/|z| with z = 1-k e^{i th}
    co = {}
    pts = []
    for s in range(L):
        th = (s + 0.5) * 2 * pi / L
        zr, zi = 1 - k * cos(th), -k * sin(th)
        r = (zr * zr + zi * zi) ** 0.5
        pts.append((th, zr / r, -zi / r))  # phi = conj(z)/|z|
    for n in range(-nmax, nmax + 1):
        re = im = 0.0
        for th, pr, pim in pts:
            c, s_ = cos(-n * th), sin(-n * th)
            re += pr * c - pim * s_; im += pr * s_ + pim * c
        co[n] = (re / L, im / L)
    return co
def fdet(M):
    M = [row[:] for row in M]; n = len(M); d = 1.0
    for c in range(n):
        piv = max(range(c, n), key=lambda r: abs(M[r][c]))
        if piv != c: M[c], M[piv] = M[piv], M[c]; d = -d
        d *= M[c][c]
        for r in range(c + 1, n):
            f = M[r][c] / M[c][c]
            for kk in range(c, n): M[r][kk] -= f * M[c][kk]
    return d
for k in (0.5, 0.9, 0.999):
    co = symbol_coeffs(k, 8)
    maxim = max(abs(co[n][1]) for n in co)
    dets = []
    for N in (1, 2, 4, 8):
        M = [[co[i - j][0] for j in range(N)] for i in range(N)]
        dets.append(round(fdet(M), 6))
    say('k=%.3f: max|Im a_n|=%.1e, D_N (N=1,2,4,8) =' % (k, maxim), dets, ' (1-k^2)^(1/4)=%.6f' % ((1 - k * k) ** 0.25))
crit = []
for N in (1, 2, 4, 8):
    v = (2 / pi) ** N
    for kk in range(1, N): v *= (kk * kk / (kk * kk - 0.25)) ** (N - kk)
    crit.append(round(v, 6))
say('critical product C_N (N=1,2,4,8) =', crit)

# ---------- 3. Lower bound C_n >= e^{-1} n^{-1/4} ----------
say('== 3. C_n lower bound')
worst = 1e9; logC = 0.0
# log C_n = n log(2/pi) + sum_{j<n} (j-n) L_j ; compute incrementally: logC_{n+1}-logC_n = log(2/pi) - sum_{j<=n} L_j
Lsum = 0.0
logC = log(2 / pi)  # n=1
for n in range(1, 200001):
    if n >= 1:
        ratio = logC - (-0.25 * log(n) - 1.0)
        worst = min(worst, ratio)
        if n in (1, 2, 10, 1000, 200000):
            say('n=%d C_n=%.6g  C_n*n^(1/4)=%.6f' % (n, exp(logC), exp(logC + 0.25 * log(n))))
    Lsum += log(1 - 1 / (4 * n * n))
    logC = logC + log(2 / pi) - Lsum
say('min over n<=2e5 of log C_n - (-(1/4)log n - 1) =', worst, '(must be >0)')
say('constant 1/2+1/36+1/4+1/16 =', 0.5 + 1 / 36 + 0.25 + 1 / 16)
say('min_t (1/2)(2t+1)^(-1/8)/((t+1)^(-1/8)) =', min(0.5 * (2 * t + 1) ** (-1 / 8) * (t + 1) ** (1 / 8) for t in range(0, 100000)))
# check log(1-x) >= -x - (2/3)x^2 on [0,1/4]
say('min over grid of log(1-x)+x+(2/3)x^2 on [0,1/4]:', min(log(1 - x) + x + (2 / 3) * x * x for x in [i / 40000 for i in range(0, 10001)]))
say('DONE')
