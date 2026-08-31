#!/usr/bin/env python3
"""
Point-level random-sieve experiment after a fixed wheel through 997.

For every prime r with 997 < r <= Q, choose one forbidden residue class
a_r mod r uniformly at random. This produces the truncated residual
Hardy-Littlewood pair factor

    G(h) = prod_r g_r(h),

where
    g_r(h) = r/(r-1)                  if r | h,
             r(r-2)/(r-1)^2          otherwise.

The controlled model then independently thins surviving points, using a
per-interval probability chosen to preserve the smooth prime-count mean.
This isolates the effect of residual pair structure from mean-model changes.

Input:
    first_1000_prime_square_interval_index.csv

The full numerical outputs used in the report are stored in:
    hl_point_model_key_results.csv
    hl_residual_singular_series_pair_check.csv
"""
import math
import numpy as np
import pandas as pd
from scipy.special import expi
from pathlib import Path

ROOT = Path(__file__).resolve().parent
df = pd.read_csv(ROOT / "first_1000_prime_square_interval_index.csv")
s = df[df["p_n"] >= 4001].copy().reset_index(drop=True)

A = s["p_n^2"].to_numpy(dtype=np.int64)
B = s["p_التالي^2"].to_numpy(dtype=np.int64)
obs = s["عدد_الأوليات_بين_p2_ومربع_التالي"].to_numpy(dtype=np.int64)

start, end = int(A[0]), int(B[-1])
N = end - start
starts, ends = A-start, B-start

def primes_upto(n):
    a = np.ones(n+1, dtype=bool)
    a[:2] = False
    for k in range(2, math.isqrt(n)+1):
        if a[k]:
            a[k*k:n+1:k] = False
    return np.flatnonzero(a)

primes = primes_upto(7919)

# Fixed exact wheel through 997.
base = np.ones(N, dtype=bool)
for q in primes[primes <= 997]:
    q = int(q)
    base[(-start) % q::q] = False

# Square endpoints are excluded because the experiment uses open intervals.
for x in np.unique(np.concatenate([A, B])):
    j = int(x)-start
    if 0 <= j < N:
        base[j] = False

baseC = np.array([
    np.count_nonzero(base[a:b]) for a,b in zip(starts,ends)
], dtype=np.int64)

# Smooth prime-count means, globally calibrated once as in prior experiments.
shape = expi(np.log(B.astype(float))) - expi(np.log(A.astype(float)))
alpha = obs.sum()/shape.sum()
mu = alpha*shape

# Main Q levels used in the stored experiment.
stages = [997,1259,1999,3163,3989,5003,6311,7919]
rng = np.random.default_rng(20260831)

# The residual pair factor after conditioning on the fixed wheel.
def residual_pair_factor(h, Q):
    tail = primes[(primes > 997) & (primes <= Q)]
    g = 1.0
    for r0 in tail:
        r = int(r0)
        if h % r == 0:
            g *= r/(r-1)
        else:
            g *= r*(r-2)/(r-1)**2
    return g

for Q in stages:
    tail = primes[(primes > 997) & (primes <= Q)]
    product = float(np.prod(1-1/tail.astype(float))) if len(tail) else 1.0

    # Expected final retention probability after the random tail sieve plus
    # interval-wise thinning. In each interval the smooth mean stays mu_i.
    theta = mu/(baseC*product)

    print(
        f"Q={Q:4d} tail_primes={len(tail):4d} "
        f"tail_survival={product:.12f} "
        f"theta_range=({theta.min():.6f},{theta.max():.6f})"
    )

# Example residual singular-series weights.
for h in [738, 2018, 3082, 7919*2]:
    print("h", h, "G_997_7919", residual_pair_factor(h, 7919))
