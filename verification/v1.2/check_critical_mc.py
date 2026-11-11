"""Refutation attempt for the critical-example lower bound (Theorem 11.8 of the paper).
Synchronous heat bath on an L x L torus, J_h=(log 3)/2, J_v=(log 2)/2, from all-plus.
Registered expectation: m_t stays ABOVE (1/2)(2t+1)^(-1/8) (the claimed bound) and decays
slowly (power-law-like, not exponentially). Also a control slightly below criticality
(couplings scaled by 0.9) should decay fast, showing the instrument can see exponential decay.
Run: python -I -B check_critical_mc.py
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from resource_guard import guard
guard()
import numpy as np

def run(L, T, scale, seed):
    rng = np.random.default_rng(seed)
    s = np.ones((L, L), dtype=np.int8)
    e2h = 3.0 ** scale   # exp(2 J_h)
    e2v = 2.0 ** scale   # exp(2 J_v)
    out = [1.0]
    for t in range(1, T + 1):
        H = (np.roll(s, 1, 1) + np.roll(s, -1, 1)).astype(np.float64)   # horizontal neighbours
        V = (np.roll(s, 1, 0) + np.roll(s, -1, 0)).astype(np.float64)   # vertical neighbours
        x = e2h ** H * e2v ** V
        p = x / (1.0 + x)
        s = np.where(rng.random((L, L)) < p, 1, -1).astype(np.int8)
        out.append(float(s.mean()))
    return np.array(out)

L, T = 1024, 400
for scale, label in ((1.0, 'critical'), (0.9, 'high-T control (couplings x0.9)')):
    ms = np.mean([run(L, T, scale, seed) for seed in (1, 2)], axis=0)
    print('==', label, flush=True)
    for t in (1, 2, 5, 10, 20, 50, 100, 200, 400):
        bound = 0.5 * (2 * t + 1) ** (-1 / 8)
        print('t=%4d  m_t=%.4f  claimed lower bound=%.4f  ok=%s' % (t, ms[t], bound, ms[t] >= bound), flush=True)
    if scale == 1.0:
        tt = np.array([50, 100, 200, 400])
        slope = np.polyfit(np.log(tt), np.log(ms[tt]), 1)[0]
        print('log-log slope over t=50..400: %.4f (theory ~ -1/(8z), z~2.17 -> -0.058)' % slope, flush=True)
