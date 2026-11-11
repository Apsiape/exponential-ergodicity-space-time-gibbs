"""Exact checks of two algebraic ingredients of the prefix-cap theorem
(Section 13 of the paper, Theorem 13.3).

(1) Cap identity with ALL exterior-child factors of the specification:
    gamma^tau_V(sigma_v=a) = p(a|tau_par(v)) Z_a / sum_b p(b|tau_par(v)) Z_b,
    Z_a = E h(X^a_T), h = product of the cap factors.
    The left side is computed by brute force from the raw specification
    product over every site u with S_u meeting V; the right side by forward
    transfer.  Also checks that every exterior child lies in the cap row.
(2) The message bound (Lemma 13.4 of the paper): for every intermediate row j
    and every pair of rows z, z' at Hamming distance 1,
    g_j(z) <= K^{|N|} g_j(z').  (Hamming distance k follows by chaining.)
    Reports the largest observed log-ratio in units of |N| log K.

Toy PCAs: binary alphabet, random strictly positive tables, random exterior.
Run all four cases:  python -I -B check_cap_identity_and_messages.py
Run one case:        python -I -B check_cap_identity_and_messages.py 0   (0..3)
The recorded output check_cap_identity_and_messages.out.txt concatenates the
four single-case runs 0, 1, 2, 3.
"""
import itertools, math, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import resource_guard
resource_guard.guard()
import numpy as np


def children(R, N):
    return {tuple(x[i] - y[i] for i in range(len(x))) for x in R for y in N}


def build(d, N, T):
    root = tuple([0] * d)
    rows = {0: [root]}
    prev = {root}
    for t in range(1, T + 1):
        cur = sorted(children(prev, N))
        rows[t] = cur
        prev = set(cur)
    return rows


def run_case(name, d, N, T, seed, brute=True):
    rng = np.random.default_rng(seed)
    nN = len(N)
    p1 = rng.uniform(0.05, 0.95, size=2 ** nN)       # p(1 | u), u encoded as bits in N order
    def p(a, u_idx):
        return p1[u_idx] if a == 1 else 1.0 - p1[u_idx]
    K = max(max(p1) / min(p1), max(1 - p1) / min(1 - p1))
    rows = build(d, N, T)
    V = [((x), t) for t in range(T + 1) for x in rows[t]]
    Vset = set(V)
    # exterior configuration on a bounding window
    rad = T + 3
    coords = list(itertools.product(range(-rad - 2, rad + 3), repeat=d))
    tau = {}
    for t in range(-1, T + 3):
        for x in coords:
            if (x, t) not in Vset:
                tau[(x, t)] = int(rng.integers(0, 2))
    def parents(site):
        x, t = site
        return [(tuple(x[i] + y[i] for i in range(d)), t - 1) for y in N]
    # sites u whose S_u meets V
    U = set(V)
    for s in V:
        x, t = s
        for y in N:
            U.add((tuple(x[i] - y[i] for i in range(d)), t + 1))
    ext_children = [u for u in U if u not in Vset]
    noncap = [u for u in ext_children if u[1] != T + 1]
    out = {"case": name, "d": d, "|N|": nN, "T": T, "|V|": len(V),
           "exterior_children_not_in_cap": len(noncap)}

    root = V[0]
    def uidx(vals):
        r = 0
        for b in vals:
            r = 2 * r + b
        return r

    # ---- right side: forward transfer -------------------------------------
    # state of row t = tuple of bits over rows[t]
    def row_states(t):
        n = len(rows[t])
        return ((np.arange(2 ** n, dtype=np.int64)[:, None] >> np.arange(n - 1, -1, -1)) & 1).astype(np.int8)

    def value_getter(t, states, a):
        pos = {x: i for i, x in enumerate(rows[t])}
        def get(site):
            x, tt = site
            if tt == 0 and (x, 0) == root:
                return np.full(states.shape[0], a, dtype=np.int8)
            if tt == t and x in pos:
                return states[:, pos[x]]
            return np.full(states.shape[0], tau[site], dtype=np.int8)
        return get

    def factor_vec(site, getter, value):
        idx = np.zeros(getter(parents(site)[0]).shape[0], dtype=np.int64)
        for par in parents(site):
            idx = 2 * idx + getter(par).astype(np.int64)
        pv = p1[idx]
        if isinstance(value, np.ndarray):
            return np.where(value == 1, pv, 1 - pv)
        return pv if value == 1 else 1 - pv

    cap_sites = [u for u in ext_children if u[1] == T + 1]

    def messages(a):
        """g[t] over row_states(t) for t=1..T (forward expectation of h)."""
        sT = row_states(T)
        get = value_getter(T, sT, a)
        h = np.ones(sT.shape[0])
        for u in cap_sites:
            h *= factor_vec(u, get, tau[u])
        g = {T: h}
        for t in range(T, 1, -1):
            sp = row_states(t - 1)
            sc = row_states(t)
            getp = value_getter(t - 1, sp, a)
            M = np.ones((sp.shape[0], sc.shape[0]))
            for i, x in enumerate(rows[t]):
                pv = factor_vec((x, t), getp, 1)              # p(1|parents) over prev states
                M *= np.where(sc[None, :, i] == 1, pv[:, None], 1 - pv[:, None])
            g[t - 1] = M @ g[t]
        # Z_a = sum over row-1 states of P(row1 | root a, tau) g[1]
        s1 = row_states(1)
        get0 = value_getter(0, np.zeros((1, 1), dtype=np.int8), a)
        w = np.ones(s1.shape[0])
        for i, x in enumerate(rows[1]):
            pv = factor_vec((x, 1), get0, 1)[0]
            w *= np.where(s1[:, i] == 1, pv, 1 - pv)
        return g, float(w @ g[1])

    gs, Z = {}, {}
    for a in (0, 1):
        gs[a], Z[a] = messages(a)
    pv_root = p1[uidx([tau[par] for par in parents(root)])]
    prior = np.array([1 - pv_root, pv_root])
    rhs = prior * np.array([Z[0], Z[1]])
    rhs = rhs / rhs.sum()
    out["root_marginal_forward"] = rhs.tolist()
    out["root_prior"] = prior.tolist()

    # ---- left side: brute force specification ------------------------------
    if brute:
        nV = len(V)
        ar = np.arange(2 ** nV, dtype=np.int32)
        conf = np.empty((2 ** nV, nV), dtype=np.int8)
        for i in range(nV):
            conf[:, i] = (ar >> (nV - 1 - i)) & 1
        del ar
        posV = {s: i for i, s in enumerate(V)}
        def getv(site):
            if site in posV:
                return conf[:, posV[site]]
            return np.full(conf.shape[0], tau[site], dtype=np.int8)
        wts = np.ones(conf.shape[0])
        for u in U:
            idx = np.zeros(conf.shape[0], dtype=np.int64)
            for par in parents(u):
                idx = 2 * idx + getv(par).astype(np.int64)
            pv = p1[idx]
            val = getv(u)
            wts *= np.where(val == 1, pv, 1 - pv)
        lhs1 = wts[conf[:, posV[root]] == 1].sum() / wts.sum()
        out["root_marginal_bruteforce"] = [1 - lhs1, lhs1]
        out["cap_identity_abs_err"] = abs(lhs1 - rhs[1])

    # ---- Lemma 1 on Hamming-1 pairs ----------------------------------------
    unit = nN * math.log(K)
    worst = 0.0
    for a in (0, 1):
        for t in range(1, T + 1):
            g = gs[a][t]
            n = len(rows[t])
            idx = np.arange(2 ** n)
            for i in range(n):
                flip = idx ^ (1 << (n - 1 - i))
                r = np.abs(np.log(g[idx]) - np.log(g[flip]))
                worst = max(worst, float(r.max()) / unit)
    out["lemma1_max_logratio_over_|N|logK"] = worst
    out["lemma1_holds"] = worst <= 1 + 1e-12
    return out


if __name__ == "__main__":
    cases = [
        ("1D radius 1, T=3", 1, [(-1,), (0,), (1,)], 3, 1, True),
        ("1D radius 1, T=5 (longer future, no brute force)", 1, [(-1,), (0,), (1,)], 5, 2, False),
        ("2D N={0,e1,e2}, T=3", 2, [(0, 0), (1, 0), (0, 1)], 3, 3, True),
        ("2D N={0,+-e1,+-e2}, T=2", 2, [(0, 0), (1, 0), (-1, 0), (0, 1), (0, -1)], 2, 4, True),
    ]
    all_ok = True
    if len(sys.argv) > 1:
        cases = [cases[int(sys.argv[1])]]
    for c in cases:
        for s in range(3):
            o = run_case(c[0], c[1], c[2], c[3], 100 * c[4] + s, c[5])
            ok = o["lemma1_holds"] and o["exterior_children_not_in_cap"] == 0 and \
                 (o.get("cap_identity_abs_err", 0.0) < 1e-12)
            all_ok &= ok
            print(o, flush=True)
    print("ALL CHECKS PASS" if all_ok else "SOME CHECK FAILED", flush=True)
