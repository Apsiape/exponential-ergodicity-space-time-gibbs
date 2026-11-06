"""Independent exact re-derivation of Table 2 (noisy majority certificates).

This program shares no code with code/majority_certificates.py and uses a different algorithm. For the noisy
majority rule, p(1 | u) = (1 - eps) maj(u) + eps/2, it computes

    q_L(eps) = P_1(origin = 1 after L steps) - P_0(origin = 1 after L steps),

where P_c is the law of the automaton started from the constant row c (the monotone identity of Proposition 4.4).
Each probability is computed exactly for rational eps = a/b by an integer transfer over the backward light cone,
one update site at a time, with every weight scaled by 2b per update site; the 0/1 symmetry
P_0(1) = 1 - P_1(1) is asserted on every call.

It checks, for every depth L = 1..7:
  - (2L+1) q_L < 1 exactly at the certified noise eps_L of Table 2, and (2L+1) q_L >= 1 at eps_L - 1/1000,
    so eps_L is the first passing grid point k/1000;
  - the three rounded values in the row of Table 2 in paper/main.tex;
  - the corresponding six-decimal lines of code/majority_certificates_output.txt and the certified points in
    code/majority_certificates_output.json;
  - the analytic high-noise bound q_L <= (2(1-eps))^L recalled in Section 11.3, at sample noises above 1/2;
and, for L = 1..3, the exact polynomials of the .json file and the cross-check values q_L(2/3), and the noise list
of Theorem D in the paper.

Run from this directory (Python 3 with SymPy; a few seconds):
    python -B table2_check.py
or with explicit paths:
    python -B table2_check.py PAPER_TEX OUTPUT_TXT OUTPUT_JSON
It prints one line per check and exits with status 0 when every check passes.
"""
import json
import re
import sys
from fractions import Fraction as Fr
from pathlib import Path


def fmt(x, d):
    """Decimal string of the nonnegative rational x rounded half up to d places."""
    n = (x * 10**d + Fr(1, 2)).__floor__()
    s = str(n).rjust(d + 1, '0')
    return s[:-d] + '.' + s[-d:]


def out_prob(L, eps, c):
    """Exact P_c(origin = 1 after L steps) as a Fraction."""
    a, b = eps.numerator, eps.denominator
    D = 2 * b
    n1 = [a, 2 * (b - a) + a]          # weight of output 1 given majority value m = 0, 1
    n0 = [D - n1[0], D - n1[1]]
    w = 2 * L - 1
    P = [0] * (1 << w)                  # time-1 row: independent sites, the majority of a constant row is the constant
    for idx in range(1 << w):
        v = 1
        for i in range(w):
            v *= n1[c] if (idx >> i) & 1 else n0[c]
        P[idx] = v
    while w > 1:                        # one more row: new site j reads old sites j, j+1, j+2
        A = P
        for j in range(w - 2):
            B = [0] * (1 << w)
            for idx in range(1 << w):
                v = A[idx]
                if not v:
                    continue
                m = 1 if ((idx >> j) & 1) + ((idx >> (j + 1)) & 1) + ((idx >> (j + 2)) & 1) >= 2 else 0
                base = idx & ~(1 << j)
                B[base | (1 << j)] += v * n1[m]
                B[base] += v * n0[m]
            A = B
        mask = (1 << (w - 2)) - 1
        P = [0] * (1 << (w - 2))
        for idx in range(1 << w):
            P[idx & mask] += A[idx]
        w -= 2
    return Fr(P[1], D ** (L * L))


def q(L, eps):
    p1, p0 = out_prob(L, eps, 1), out_prob(L, eps, 0)
    assert p0 == 1 - p1, ('symmetry', L, eps)
    return p1 - p0


def main():
    here = Path(__file__).resolve().parent
    if len(sys.argv) == 4:
        tex, txt, js = (Path(x) for x in sys.argv[1:4])
    else:
        tex = here.parent.parent / 'paper' / 'main.tex'
        txt = here.parent / 'code' / 'majority_certificates_output.txt'
        js = here.parent / 'code' / 'majority_certificates_output.json'
    src = tex.read_text(encoding='utf-8')
    rows = {}
    for m in re.finditer(r'^(\d) & \$(\d+/\d+)\$ & \$([\d.]+)\$ & \$([\d.]+)\$ & \$([\d.]+)\$\\\\$', src, re.M):
        rows[int(m.group(1))] = (Fr(m.group(2)), m.group(3), m.group(4), m.group(5))
    m7 = re.search(r'\\newcommand\{\\MAJSEVEN\}\{7 & \$(\d+/\d+)\$ & \$([\d.]+)\$', src)
    a7 = re.search(r'\\newcommand\{\\MAJSEVENA\}\{([\d.]+)\}', src).group(1)
    b7 = re.search(r'\\newcommand\{\\MAJSEVENB\}\{([\d.]+)\}', src).group(1)
    rows[7] = (Fr(m7.group(1)), m7.group(2), a7, b7)
    assert sorted(rows) == list(range(1, 8)), rows
    output = txt.read_text(encoding='utf-8')
    data = json.loads(js.read_text(encoding='utf-8'))
    bad = 0

    def check(cond, msg):
        nonlocal bad
        print(('ok   ' if cond else 'FAIL ') + msg)
        bad += 0 if cond else 1

    for L in range(1, 8):
        eL, s_cert, s35, s12 = rows[L]
        k = 2 * L + 1
        vc, v35, v12 = k * q(L, eL), k * q(L, Fr(3, 5)), k * q(L, Fr(1, 2))
        below = k * q(L, eL - Fr(1, 1000))
        check(vc < 1, f'L={L}: (2L+1)q_L({eL}) = {float(vc):.9f} < 1 exactly')
        check(below >= 1, f'L={L}: (2L+1)q_L({eL - Fr(1, 1000)}) = {float(below):.9f} >= 1, so {eL} is the first passing grid point')
        check(fmt(vc, 4) == s_cert and fmt(v35, 4) == s35 and fmt(v12, 4) == s12,
              f'L={L}: paper row {s_cert} {s35} {s12} vs recomputed {fmt(vc, 4)} {fmt(v35, 4)} {fmt(v12, 4)}')
        line = f'L={L}: (2L+1) q_L < 1 holds exactly at eps = {eL}:  (2L+1) q_L = {fmt(vc, 6)}'
        check(line in output, f'L={L}: recorded output has "{line}"')
        line = f'L={L}: eps=3/5: {fmt(v35, 6)}   eps=1/2: {fmt(v12, 6)}'
        check(line in output, f'L={L}: recorded output has "{line}"')
        check(data['deeper_certs'][str(L)] == str(eL), f'L={L}: json certified point {data["deeper_certs"][str(L)]}')
        for e in sorted({Fr(1, 2), Fr(527, 1000), Fr(3, 5), eL}):
            if e >= Fr(1, 2):
                check(q(L, e) <= (2 * (1 - e)) ** L,
                      f'L={L}: q_L({e}) = {float(q(L, e)):.6f} <= (2(1-eps))^L = {float((2 * (1 - e)) ** L):.6f}')
    import sympy as sp
    epsilon = sp.symbols('epsilon')
    for L in (1, 2, 3):
        poly = sp.sympify(data[str(L)]['q_poly'], locals={'epsilon': epsilon})
        for e in (Fr(1, 2), Fr(3, 5), Fr(2, 3), Fr(587, 1000), Fr(7, 10)):
            pv = Fr(str(sp.Rational(poly.subs(epsilon, sp.Rational(e.numerator, e.denominator)))))
            check(pv == q(L, e), f'L={L}: json polynomial at {e} = {pv} equals the transfer')
    for L, val in ((1, '1/3'), (2, '13/81'), (3, '1529/19683')):
        check(q(L, Fr(2, 3)) == Fr(val) and f'L={L}: q_L(2/3) by row transfer = {val}' in output, f'L={L}: q_L(2/3) = {val}')
    noise = ', '.join(f'${float(rows[L][0]):.3f}$' for L in range(1, 7)) + f' and ${float(rows[7][0]):.3f}$, at depths one to seven'
    check(noise in src, f'Theorem D lists the certified noises: {noise}')
    print('FAILURES:', bad)
    sys.exit(1 if bad else 0)


main()
