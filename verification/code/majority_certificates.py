#!/usr/bin/env python3
"""
Exact finite certificates for the noisy majority PCA (radius 1, binary alphabet).

Local-map representation.  At every space-time update, independently,
  with probability 1 - eps : apply the deterministic majority rule  f(u) = maj(u_{-1}, u_0, u_1);
  with probability eps/2   : output the constant 0;
  with probability eps/2   : output the constant 1.
The induced transition kernel is p(a | u) = (1 - eps) 1[a = f(u)] + eps/2 > 0.

For depth L, the output at the origin after L updates is a map
    g : {0,1}^{[-L, L]} -> {0,1}
determined by the innovation configuration c in {D, N0, N1}^{L^2}.

  q_L(eps)     = P( g is nonconstant )                                   (common-map coefficient)
  theta_L(eps) = max over input pairs (i, j) of P( g(i) != g(j) )        (fixed-pair coefficient)

Both are polynomials in eps with rational coefficients.  We enumerate every innovation
configuration, accumulate integer counts grouped by the number k of noisy updates, and
evaluate exactly in rational arithmetic.  The certificates used in the paper are

  common map : (2L+1) q_L(eps)     < 1
  fixed pair : (2L+1) theta_L(eps) < 1

Run:  python3 majority_certificates.py            (writes majority_certificates_output.txt)
"""
import itertools, sys, time, json, hashlib
from fractions import Fraction as Fr
import numpy as np
import sympy as sp

def maj(a, b, c):
    return ((a.astype(np.int16) + b + c) >= 2).astype(np.uint8)

def cone_sites(L):
    """Update sites of the depth-L backward cone of the origin, as (time, x), time = 1..L."""
    sites = []
    for s in range(1, L + 1):
        half = L - s
        for x in range(-half, half + 1):
            sites.append((s, x))
    return sites

def enumerate_counts(L):
    """Return (n_in, count_nc[k], count_pair[k][i,j]) with k = number of noisy updates."""
    sites = cone_sites(L)
    nsites = len(sites)
    assert nsites == L * L
    width = 2 * L + 1
    n_in = 2 ** width
    # all input rows, shape (n_in, width); column index x+L holds site x
    inputs = np.array([[(m >> (width - 1 - j)) & 1 for j in range(width)] for m in range(n_in)], dtype=np.uint8)
    count_nc = np.zeros(nsites + 1, dtype=np.int64)
    count_pair = np.zeros((nsites + 1, n_in, n_in), dtype=np.int64)
    zeros = np.zeros(n_in, dtype=np.uint8)
    ones = np.ones(n_in, dtype=np.uint8)
    t0 = time.time()
    for conf in itertools.product((0, 1, 2), repeat=nsites):   # 0 = D, 1 = N0, 2 = N1
        k = sum(1 for c in conf if c != 0)
        # row values are stored in dict keyed by (time, x); time 0 = inputs
        prev = {x: inputs[:, x + L] for x in range(-L, L + 1)}
        idx = 0
        for s in range(1, L + 1):
            half = L - s
            cur = {}
            for x in range(-half, half + 1):
                c = conf[idx]; idx += 1
                if c == 0:
                    cur[x] = maj(prev[x - 1], prev[x], prev[x + 1])
                elif c == 1:
                    cur[x] = zeros
                else:
                    cur[x] = ones
            prev = cur
        out = prev[0]
        if out.min() != out.max():
            count_nc[k] += 1
        count_pair[k] += (out[:, None] != out[None, :])
    sys.stderr.write(f"L={L}: enumerated {3**nsites} configurations in {time.time()-t0:.1f}s\n")
    return n_in, count_nc, count_pair

def poly_from_counts(counts, nsites):
    """Exact sympy polynomial  sum_k counts[k] (eps/2)^k (1-eps)^(nsites-k)."""
    e = sp.symbols('epsilon')
    expr = sum(sp.Integer(int(counts[k])) * (e / 2) ** k * (1 - e) ** (nsites - k) for k in range(nsites + 1))
    return sp.expand(expr)

def eval_exact(counts, nsites, eps):
    eps = Fr(eps)
    return sum(Fr(int(counts[k])) * (eps / 2) ** k * (1 - eps) ** (nsites - k) for k in range(nsites + 1))

def theta_exact(count_pair, nsites, eps):
    """Exact max over all pairs at a rational eps."""
    eps = Fr(eps)
    w = [(eps / 2) ** k * (1 - eps) ** (nsites - k) for k in range(nsites + 1)]
    n_in = count_pair.shape[1]
    best = Fr(0); best_pair = None
    for i in range(n_in):
        for j in range(i + 1, n_in):
            v = sum(Fr(int(count_pair[k, i, j])) * w[k] for k in range(nsites + 1))
            if v > best:
                best, best_pair = v, (i, j)
    return best, best_pair

def theta_float_curve(count_pair, nsites, grid):
    """Float evaluation of theta_L on a grid, for locating thresholds (certified separately)."""
    n_in = count_pair.shape[1]
    vals = []
    for eps in grid:
        w = np.array([(eps / 2) ** k * (1 - eps) ** (nsites - k) for k in range(nsites + 1)])
        M = np.tensordot(w, count_pair, axes=(0, 0))
        vals.append(M.max())
    return np.array(vals)

def main():
    e = sp.symbols('epsilon')
    lines = []
    def out(s=""):
        lines.append(s); print(s)
    out("Exact finite certificates for the noisy majority PCA (radius 1, alphabet {0,1}).")
    out("Kernel: p(a|u) = (1-eps) 1[a = maj(u)] + eps/2.  Certificates: (2L+1) q_L < 1 and (2L+1) theta_L < 1.")
    out()
    results = {}
    for L in (1, 2, 3):
        nsites = L * L
        n_in, count_nc, count_pair = enumerate_counts(L)
        qpoly = poly_from_counts(count_nc, nsites)
        out(f"== depth L = {L}: {nsites} update sites, {n_in} input rows, {3**nsites} innovation configurations ==")
        out(f"q_{L}(eps) = {sp.factor(qpoly)}")
        out(f"         = {qpoly}")
        # exact evaluation at a few rational points
        pts = [Fr(1, 2), Fr(3, 5), Fr(5, 8), Fr(2, 3), Fr(7, 10), Fr(3, 4)]
        out("  eps      q_L(eps)            (2L+1)q_L        theta_L(eps)         (2L+1)theta_L   max pair   trivial [3(1-eps)]^L")
        for eps in pts:
            q = eval_exact(count_nc, nsites, eps)
            th, pair = theta_exact(count_pair, nsites, eps)
            triv = (3 * (1 - eps)) ** L
            out(f"  {str(eps):6s}   {str(q):18s}  {float((2*L+1)*q):.6f}         {str(th):18s}  {float((2*L+1)*th):.6f}        {pair}   {float(triv):.6f}")
        # thresholds on a fine grid (float), then exact certification at the first rational grid point that passes
        grid = np.arange(0.300, 0.700, 0.001)
        qf = np.array([float(eval_exact(count_nc, nsites, Fr(round(x*1000), 1000))) for x in grid])
        thf = theta_float_curve(count_pair, nsites, grid)
        def first_pass(vals):
            ok = np.where((2 * L + 1) * vals < 1)[0]
            return None if len(ok) == 0 else ok[0]
        iq = first_pass(qf); ith = first_pass(thf)
        rec = {"L": L, "q_poly": str(qpoly)}
        if iq is not None:
            eq = Fr(round(grid[iq] * 1000), 1000)
            # move up by one grid step and certify exactly; also check monotone decrease on the grid above it
            qex = eval_exact(count_nc, nsites, eq)
            assert (2 * L + 1) * qex < 1
            mono = bool(np.all(np.diff(qf[iq:]) <= 1e-12))
            out(f"  common-map certificate: (2L+1) q_L < 1 holds exactly at eps = {eq} (value {float((2*L+1)*qex):.6f}); q_L decreasing on the grid above it: {mono}")
            rec["q_cert_eps"] = str(eq)
        if ith is not None:
            et = Fr(round(grid[ith] * 1000), 1000)
            thex, pair = theta_exact(count_pair, nsites, et)
            assert (2 * L + 1) * thex < 1
            mono = bool(np.all(np.diff(thf[ith:]) <= 1e-12))
            out(f"  fixed-pair certificate: (2L+1) theta_L < 1 holds exactly at eps = {et} (value {float((2*L+1)*thex):.6f}, pair {pair}); theta_L decreasing on the grid above it: {mono}")
            rec["theta_cert_eps"] = str(et)
        out()
        results[L] = rec
    out("Note. Each certificate is verified in exact rational arithmetic at the stated point. Monotonicity of q_L and theta_L in eps")
    out("is reported on the grid only and is not used by the paper.")
    text = "\n".join(lines) + "\n"
    with open("majority_certificates_output.txt", "w") as f:
        f.write(text)
    with open("majority_certificates_output.json", "w") as f:
        json.dump(results, f, indent=1)
    print("SHA-256(majority_certificates_output.txt) =", hashlib.sha256(text.encode()).hexdigest())

if __name__ == "__main__":
    main()


# ---------------------------------------------------------------------------------------
# Part 2.  Monotone shortcut and deeper certificates.
#
# Every local map (majority or a constant) is monotone, so the composite map g is monotone.
# A monotone Boolean map is nonconstant iff g(0...0) = 0 and g(1...1) = 1, and any pair of
# inputs that disagree under g forces nonconstancy.  Hence, for monotone rules,
#       theta_L = q_L = P(g(1...1) = 1) - P(g(0...0) = 1).
# The two probabilities are computed by an exact row-by-row transfer over the cone, which
# reaches depths beyond brute-force enumeration.  We cross-check against Part 1 for L <= 3.
# ---------------------------------------------------------------------------------------

def bits_matrix(w):
    return np.array([[(z >> (w - 1 - j)) & 1 for j in range(w)] for z in range(2 ** w)], dtype=np.int64)

def row_dp_float(L, eps):
    """(P(g(1..1)=1), P(g(0..0)=1)) in float64 via vectorized row transfer."""
    def P1(m):  # P(output 1 | majority value m), m array
        return (1 - eps) * m + eps * 0.5
    res = []
    for const in (1, 0):
        w1 = 2 * L - 1
        Z = bits_matrix(w1)
        p1 = P1(const)
        dist = np.prod(np.where(Z == 1, p1, 1 - p1), axis=1)          # row 1: iid sites
        for s in range(2, L + 1):
            wprev = 2 * (L - s + 1) + 1; w = 2 * (L - s) + 1
            X = bits_matrix(wprev)                                        # (nx, wprev)
            M = ((X[:, :-2] + X[:, 1:-1] + X[:, 2:]) >= 2).astype(float)  # (nx, w) majority of triples
            pj = P1(M)                                                    # (nx, w)
            Zb = bits_matrix(w)                                           # (nz, w)
            T = np.ones((X.shape[0], Zb.shape[0]))
            for j in range(w):
                T *= np.where(Zb[None, :, j] == 1, pj[:, j, None], 1 - pj[:, j, None])
            dist = dist @ T
        res.append(dist[1])
    return res[0], res[1]

def row_dp_exact(L, eps):
    """Same transfer in exact rational arithmetic (Fractions)."""
    eps = Fr(eps)
    def P1(m): return (1 - eps) * m + eps / 2
    res = []
    for const in (1, 0):
        w1 = 2 * L - 1
        p1 = P1(const)
        dist = {}
        for y in range(2 ** w1):
            pr = Fr(1)
            for i in range(w1):
                pr *= p1 if (y >> (w1 - 1 - i)) & 1 else (1 - p1)
            dist[y] = pr
        for s in range(2, L + 1):
            wprev = 2 * (L - s + 1) + 1; w = 2 * (L - s) + 1
            Zbits = [[(z >> (w - 1 - j)) & 1 for j in range(w)] for z in range(2 ** w)]
            new = [Fr(0)] * (2 ** w)
            for x, px in dist.items():
                if px == 0: continue
                xb = [(x >> (wprev - 1 - i)) & 1 for i in range(wprev)]
                pj = [P1(1 if (xb[j] + xb[j + 1] + xb[j + 2]) >= 2 else 0) for j in range(w)]
                qj = [1 - p for p in pj]
                for z in range(2 ** w):
                    pr = px
                    zb = Zbits[z]
                    for j in range(w):
                        pr *= pj[j] if zb[j] else qj[j]
                    new[z] += pr
            dist = {z: v for z, v in enumerate(new)}
        res.append(dist[1])
    return res[0], res[1]

def deeper(Lmax=7):
    lines = []
    def out(s=""):
        lines.append(s); print(s, flush=True)
    out()
    out("== Part 2: monotone shortcut, theta_L = q_L = P(g(1..1)=1) - P(g(0..0)=1), by exact row transfer ==")
    out("Cross-check against brute force (Part 1) at eps = 2/3:")
    for L in (1, 2, 3):
        a, b = row_dp_exact(L, Fr(2, 3))
        out(f"  L={L}: q_L(2/3) by row transfer = {a - b}")
    out("Certified fixed-pair thresholds on the grid eps = k/1000 (first grid point passing in floating point, then verified exactly):")
    certs = {}
    for L in range(1, Lmax + 1):
        t0 = time.time()
        first = None
        for k in range(300, 700):
            x = k / 1000
            a, b = row_dp_float(L, x)
            if (2 * L + 1) * (a - b) < 1 - 1e-9:
                first = k; break
        if first is None:
            out(f"  L={L}: no grid point below 0.700 passes"); continue
        ef = Fr(first, 1000)
        a, b = row_dp_exact(L, ef)
        q = a - b
        assert (2 * L + 1) * q < 1, (L, ef, q)
        certs[L] = str(ef)
        out(f"  L={L}: (2L+1) q_L < 1 holds exactly at eps = {ef}:  (2L+1) q_L = {float((2*L+1)*q):.6f}   [{time.time()-t0:.1f}s]")
    out("Exact values of (2L+1) q_L at eps = 3/5 and eps = 1/2:")
    for L in range(1, Lmax + 1):
        vals = []
        for e_ in (Fr(3, 5), Fr(1, 2)):
            a, b = row_dp_exact(L, e_)
            vals.append((e_, a - b))
        out("  L=%d: " % L + "   ".join(f"eps={e_}: {float((2*L+1)*v):.6f}" for e_, v in vals))
    text = "\n".join(lines) + "\n"
    with open("majority_certificates_output.txt", "a") as f:
        f.write(text)
    with open("majority_certificates_output.json") as f:
        res = json.load(f)
    res["deeper_certs"] = certs
    with open("majority_certificates_output.json", "w") as f:
        json.dump(res, f, indent=1)
    full = open("majority_certificates_output.txt").read()
    print("SHA-256(majority_certificates_output.txt) =", hashlib.sha256(full.encode()).hexdigest())

if __name__ == "__main__":
    deeper(Lmax=int(sys.argv[1]) if len(sys.argv) > 1 else 7)
