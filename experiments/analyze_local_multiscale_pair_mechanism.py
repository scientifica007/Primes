#!/usr/bin/env python3
"""
Local multiscale pair mechanism.

Purpose
-------
Replace the interval-level Gaussian pair-kernel correction by a genuinely
point-level, finite-support multiscale mechanism.

The mechanism:
  1. fixed wheel q0=7;
  2. random forbidden residue classes for 7 < p <= 7919;
  3. local zero-sum Haar redistributions at dyadic scales;
  4. theory-only NNLS fit of scale variances to the short-interval pair kernel;
  5. one scalar amplitude calibrated on training data using only R;
  6. held-out validation of rho(k), runs, and triple sign patterns.

No explicit k=3 information is used.

Repository input
----------------
first_1000_prime_square_interval_index.csv
"""

import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import nnls
from scipy.sparse import coo_matrix
from scipy.special import expi

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "first_1000_prime_square_interval_index.csv"

BASE_Q = 7
Q = 7919
CELL = 1024
MAX_LAG = 20

RAW_REALIZATIONS = 40
TRAIN_RAW = 12
TRAIN_FIELDS_PER_RAW = 25
TEST_FIELDS_PER_RAW = 24

RAW_SEEDS = [20261010, 20261011]
TRAIN_FIELD_SEED = 20261600
TEST_FIELD_SEED = 20261700

df = pd.read_csv(DATA)
s = df[df["p_n"] >= 4001].copy().reset_index(drop=True)

p = s["p_n"].to_numpy(np.int64)
p_next = s["p_التالي"].to_numpy(np.int64)
A = s["p_n^2"].to_numpy(np.int64)
B = s["p_التالي^2"].to_numpy(np.int64)
obs = s["عدد_الأوليات_بين_p2_ومربع_التالي"].to_numpy(float)

start, end = int(A[0]), int(B[-1])
N = end - start
starts = A - start
ends = B - start
M = len(A)


def primes_upto(limit):
    z = np.ones(limit + 1, dtype=bool)
    z[:2] = False
    for k in range(2, math.isqrt(limit) + 1):
        if z[k]:
            z[k*k:limit+1:k] = False
    return np.flatnonzero(z)


def lag_corrs(Z, max_lag):
    Z = np.asarray(Z, float)
    if Z.ndim == 1:
        Z = Z[None, :]
    out = np.empty((Z.shape[0], max_lag))
    for k in range(1, max_lag + 1):
        x = Z[:, :-k]
        y = Z[:, k:]
        xc = x - x.mean(axis=1, keepdims=True)
        yc = y - y.mean(axis=1, keepdims=True)
        out[:, k-1] = np.sum(xc * yc, axis=1) / np.sqrt(
            np.sum(xc * xc, axis=1) * np.sum(yc * yc, axis=1)
        )
    return out


def Delta(t):
    t = float(t)
    return 0.5 * (
        (t + 1) * math.log(t + 1)
        - 2 * t * math.log(t)
        + (0.0 if t == 1 else (t - 1) * math.log(t - 1))
    )


primes = primes_upto(Q)

shape = expi(np.log(B.astype(float))) - expi(np.log(A.astype(float)))
alpha = obs.sum() / shape.sum()
mu = alpha * shape

# ------------------------------------------------------------
# Fixed wheel and tail.
# ------------------------------------------------------------

base = np.ones(N, dtype=bool)

for q0 in primes[primes <= BASE_Q]:
    q = int(q0)
    base[(-start) % q::q] = False

for x in np.unique(np.concatenate([A, B])):
    j = int(x) - start
    if 0 <= j < N:
        base[j] = False

baseC = np.add.reduceat(base, starts).astype(np.int64)
base_var = mu * (1 - mu / baseC)
sd_base = np.sqrt(base_var)
z_obs = (obs - mu) / sd_base

R_obs = float(np.sum(z_obs * z_obs) / (M - 1))
obs_rho = lag_corrs(z_obs, MAX_LAG)[0]

tail = primes[(primes > BASE_Q) & (primes <= Q)]
tail_product = float(np.prod(1 - 1 / tail.astype(float)))
theta = mu / (baseC * tail_product)

# ------------------------------------------------------------
# Fine local segments: regular 1024-integer cells, split also at
# prime-square interval boundaries.
# ------------------------------------------------------------

regular = np.arange(0, N + CELL, CELL, dtype=np.int64)
regular = regular[regular <= N]
if regular[-1] != N:
    regular = np.append(regular, N)

boundaries = np.unique(
    np.concatenate([regular, starts, ends, [0, N]])
).astype(np.int64)

seg_starts = boundaries[:-1]
seg_ends = boundaries[1:]
seg_lengths = seg_ends - seg_starts

seg_interval = np.searchsorted(
    ends, seg_starts, side="right"
).astype(np.int16)

seg_cell = (seg_starts // CELL).astype(np.int32)
ncell = int(math.ceil(N / CELL))
cell_index = np.arange(ncell, dtype=np.int64)

seg_theta = theta[seg_interval]


def interval_sum_from_segments(values):
    return np.bincount(
        seg_interval, weights=values, minlength=M
    )


# ------------------------------------------------------------
# Generate 40 point-level random-sieve realizations and retain
# survivor counts per local segment.
# ------------------------------------------------------------

S_parts = []

for seed in RAW_SEEDS:
    rng = np.random.default_rng(seed)
    part = []
    for _ in range(RAW_REALIZATIONS // len(RAW_SEEDS)):
        mask = base.copy()
        for rv in tail:
            r = int(rv)
            forbidden = int(rng.integers(0, r))
            offset = (forbidden - start) % r
            mask[offset:N:r] = False
        part.append(
            np.add.reduceat(mask, seg_starts).astype(np.uint16)
        )
    S_parts.append(np.asarray(part))

Sseg = np.vstack(S_parts).astype(np.float64)

TRAIN_IDX = np.arange(0, TRAIN_RAW)
TEST_IDX = np.arange(TRAIN_RAW, RAW_REALIZATIONS)

# ------------------------------------------------------------
# Short-interval pair theory for the actual variable intervals.
# ------------------------------------------------------------

centers = 0.5 * (A + B)


def V(h, X):
    if h <= 0:
        return 0.0
    return h * math.log(X / h)


def interval_pair_rho(i, j):
    a, b = float(A[i]), float(B[i])
    c, d = float(A[j]), float(B[j])
    Xloc = math.sqrt(float(centers[i]) * float(centers[j]))

    cov = 0.5 * (
        V(d - a, Xloc)
        + V(c - b, Xloc)
        - V(c - a, Xloc)
        - V(d - b, Xloc)
    )

    return cov / math.sqrt(
        V(b - a, Xloc) * V(d - c, Xloc)
    )


theory_rho = np.array([
    np.mean([
        interval_pair_rho(i, i + k)
        for i in range(M - k)
    ])
    for k in range(1, MAX_LAG + 1)
])

delta = np.array([
    Delta(k) for k in range(1, MAX_LAG + 1)
])
delta_norm = delta / delta[0]

# ------------------------------------------------------------
# Exact covariance basis of local zero-sum Haar transfers.
#
# For each scale, average over four grid translations.
# The basis is computed from geometry only; no prime counts enter.
# ------------------------------------------------------------

scales = np.array([2**j for j in range(1, 16)], dtype=int)
interval_len = (ends - starts).astype(float)

basis_rows = []

for L in scales:
    shifts = sorted(set([
        0, L // 4, L // 2, (3 * L) // 4
    ]))

    cov_acc = np.zeros(MAX_LAG + 1)

    for shift in shifts:
        u = seg_cell.astype(np.int64) + shift
        block = u // L
        pos = u % L
        sign = np.where(pos < L / 2, 1.0, -1.0)

        data = (
            seg_lengths.astype(float)
            / interval_len[seg_interval]
        ) * sign

        nblocks = int(block.max()) + 1

        W = coo_matrix(
            (data, (seg_interval, block)),
            shape=(M, nblocks)
        ).tocsr()

        cov_acc[0] += float(np.mean(
            np.asarray(W.multiply(W).sum(axis=1)).ravel()
        ))

        for k in range(1, MAX_LAG + 1):
            vals = np.asarray(
                W[:-k].multiply(W[k:]).sum(axis=1)
            ).ravel()
            cov_acc[k] += float(np.mean(vals))

    cov_acc /= len(shifts)
    basis_rows.append(cov_acc)

basis = np.asarray(basis_rows).T

# Target covariance: unit variance plus short-interval correlations.
target = np.concatenate([[1.0], theory_rho])

lag_weights = np.ones(MAX_LAG + 1)
lag_weights[0] = 3.0
lag_weights[1:6] = 2.0

variance_weights, nnls_residual = nnls(
    basis * lag_weights[:, None],
    target * lag_weights
)

fit_cov = basis @ variance_weights
fit_corr = fit_cov[1:] / fit_cov[0]

scale_table = pd.DataFrame({
    "scale_cells": scales,
    "scale_integers": scales * CELL,
    "variance_weight_v": variance_weights,
})

scale_table.to_csv(
    ROOT / "local_multiscale_scale_weights_q7.csv",
    index=False, encoding="utf-8-sig"
)

fit_table = pd.DataFrame({
    "lag": np.arange(0, MAX_LAG + 1),
    "target_covariance": target,
    "fitted_local_multiscale_covariance": fit_cov,
})
fit_table["absolute_error"] = np.abs(
    fit_table["target_covariance"]
    - fit_table["fitted_local_multiscale_covariance"]
)

fit_table.to_csv(
    ROOT / "local_multiscale_theory_kernel_fit_q7.csv",
    index=False, encoding="utf-8-sig"
)

active = np.where(variance_weights > 1e-10)[0]
active_scales = scales[active]
active_v = variance_weights[active]


def local_multiscale_field(rng):
    field = np.zeros(ncell, dtype=float)

    for L, var_weight in zip(active_scales, active_v):
        shifts = sorted(set([
            0, L // 4, L // 2, (3 * L) // 4
        ]))

        shift = int(rng.choice(shifts))

        u = cell_index + shift
        block = u // L
        pos = u % L

        sign = np.where(pos < L / 2, 1.0, -1.0)

        coeff = rng.choice(
            np.array([-1.0, 1.0]),
            size=int(block.max()) + 1
        )

        field += (
            math.sqrt(var_weight)
            * coeff[block]
            * sign
        )

    return field


# ------------------------------------------------------------
# Calibrate ONE amplitude using only R on the 12 training raw sieves.
# ------------------------------------------------------------

rng = np.random.default_rng(TRAIN_FIELD_SEED)

c0_terms = []
c1_terms = []
c2_terms = []

for ridx in TRAIN_IDX:
    ss = Sseg[ridx]

    S_interval = interval_sum_from_segments(ss)
    Y0 = theta * S_interval
    Z0 = (Y0 - mu) / sd_base

    for _ in range(TRAIN_FIELDS_PER_RAW):
        field = local_multiscale_field(rng)

        delta_count = interval_sum_from_segments(
            field[seg_cell] * ss
        )

        Zd = delta_count / sd_base

        c0_terms.append(
            np.sum(Z0 * Z0) / (M - 1)
        )
        c1_terms.append(
            np.sum(Z0 * Zd) / (M - 1)
        )
        c2_terms.append(
            np.sum(Zd * Zd) / (M - 1)
        )

c0 = float(np.mean(c0_terms))
c1 = float(np.mean(c1_terms))
c2 = float(np.mean(c2_terms))

disc = (
    (2 * c1) ** 2
    - 4 * c2 * (c0 - R_obs)
)

roots = [
    (-(2 * c1) + math.sqrt(disc)) / (2 * c2),
    (-(2 * c1) - math.sqrt(disc)) / (2 * c2),
]

amplitude = min(r for r in roots if r >= 0)

# ------------------------------------------------------------
# Held-out validation: 28 raw sieves x 24 local fields = 672 runs.
# ------------------------------------------------------------

rng = np.random.default_rng(TEST_FIELD_SEED)

Rvals = []
runs = []
rhos = []
monotone = []
pattern_counts = []
meanY = np.zeros(M)
clip_count = 0
total_probabilities = 0

for ridx in TEST_IDX:
    ss = Sseg[ridx]

    for _ in range(TEST_FIELDS_PER_RAW):
        field = local_multiscale_field(rng)

        pseg = (
            seg_theta
            + amplitude * field[seg_cell]
        )

        clip_count += int(np.sum(
            (pseg < 0) | (pseg > 1)
        ))
        total_probabilities += len(pseg)

        pseg = np.clip(pseg, 0, 1)

        target_counts = pseg * ss
        lo = np.floor(target_counts)
        frac = target_counts - lo

        K = lo + rng.binomial(1, frac)

        Y = np.bincount(
            seg_interval,
            weights=K,
            minlength=M
        )

        meanY += Y

        Z = (Y - mu) / sd_base

        Rvals.append(
            np.sum(Z * Z) / (M - 1)
        )

        sg = Z >= 0

        runs.append(
            1 + np.sum(
                sg[1:] != sg[:-1]
            )
        )

        rhos.append(
            lag_corrs(Z, MAX_LAG)[0]
        )

        monotone.append(np.sum(
            ((Z[:-2] < Z[1:-1]) & (Z[1:-1] < Z[2:]))
            |
            ((Z[:-2] > Z[1:-1]) & (Z[1:-1] > Z[2:]))
        ))

        codes = (
            sg[:-2].astype(np.int8) * 4
            + sg[1:-1].astype(np.int8) * 2
            + sg[2:].astype(np.int8)
        )

        pattern_counts.append(
            np.bincount(codes, minlength=8)
        )

Rvals = np.asarray(Rvals)
runs = np.asarray(runs)
rhos = np.asarray(rhos)
monotone = np.asarray(monotone)
pattern_counts = np.asarray(pattern_counts)

meanrho = rhos.mean(axis=0)
sdrho = rhos.std(axis=0, ddof=1)


def two_sided_p(obs_value, simulations):
    simulations = np.asarray(simulations)
    center = simulations.mean()

    return float(
        (
            1
            + np.sum(
                np.abs(simulations - center)
                >= abs(obs_value - center)
            )
        )
        / (len(simulations) + 1)
    )


amp_delta = float(
    np.dot(-delta, meanrho)
    / np.dot(delta, delta)
)

shape_rmse = float(np.sqrt(np.mean(
    (
        meanrho / meanrho[0]
        - delta_norm
    ) ** 2
)))

kernel_z_rmse = float(np.sqrt(np.mean(
    (
        (meanrho - theory_rho)
        / sdrho
    ) ** 2
)))

obs_runs = int(
    1 + np.sum(
        (z_obs[1:] >= 0)
        != (z_obs[:-1] >= 0)
    )
)

trip = np.column_stack([
    z_obs[:-2],
    z_obs[1:-1],
    z_obs[2:]
])

obs_monotone = int(np.sum(
    (
        (trip[:, 0] < trip[:, 1])
        & (trip[:, 1] < trip[:, 2])
    )
    |
    (
        (trip[:, 0] > trip[:, 1])
        & (trip[:, 1] > trip[:, 2])
    )
))

summary = pd.DataFrame([
    ["simulations", len(Rvals)],
    ["amplitude_calibrated_from_train_R", amplitude],
    ["R_observed", R_obs],
    ["R_mc_mean", float(Rvals.mean())],
    ["R_mc_sd", float(Rvals.std(ddof=1))],
    ["R_p", two_sided_p(R_obs, Rvals)],
    ["rho1_observed", float(obs_rho[0])],
    ["rho1_theory", float(theory_rho[0])],
    ["rho1_mc_mean", float(meanrho[0])],
    ["rho2_theory", float(theory_rho[1])],
    ["rho2_mc_mean", float(meanrho[1])],
    ["delta_amplitude", amp_delta],
    ["normalized_Delta_shape_RMSE", shape_rmse],
    ["kernel_z_rmse_vs_theory", kernel_z_rmse],
    ["runs_observed", obs_runs],
    ["runs_mc_mean", float(runs.mean())],
    ["runs_mc_sd", float(runs.std(ddof=1))],
    ["runs_p", two_sided_p(obs_runs, runs)],
    ["monotone_triples_observed", obs_monotone],
    ["monotone_triples_mc_mean", float(monotone.mean())],
    ["monotone_triples_mc_sd", float(monotone.std(ddof=1))],
    ["monotone_triples_p", two_sided_p(
        obs_monotone, monotone
    )],
    ["clip_fraction", clip_count / total_probabilities],
    ["mean_bias_z_rms", float(np.sqrt(np.mean(
        (
            (meanY / len(Rvals) - mu)
            / sd_base
        ) ** 2
    )))],
    ["nnls_theory_kernel_residual", nnls_residual],
], columns=["metric", "value"])

summary.to_csv(
    ROOT / "local_multiscale_point_model_summary_q7.csv",
    index=False, encoding="utf-8-sig"
)

# ------------------------------------------------------------
# Full rho(k) tests.
# ------------------------------------------------------------

obs_z = (obs_rho - meanrho) / sdrho

T_obs = float(np.sqrt(np.mean(obs_z ** 2)))
M_obs = float(np.max(np.abs(obs_z)))

sim_z = (
    rhos - meanrho[None, :]
) / sdrho[None, :]

T_sim = np.sqrt(np.mean(sim_z ** 2, axis=1))
M_sim = np.max(np.abs(sim_z), axis=1)

global_test = pd.DataFrame([
    [
        "rho_kernel_RMS_standardized",
        T_obs,
        float(T_sim.mean()),
        float(T_sim.std(ddof=1)),
        float(
            (1 + np.sum(T_sim >= T_obs))
            / (len(T_sim) + 1)
        )
    ],
    [
        "rho_kernel_max_abs_standardized",
        M_obs,
        float(M_sim.mean()),
        float(M_sim.std(ddof=1)),
        float(
            (1 + np.sum(M_sim >= M_obs))
            / (len(M_sim) + 1)
        )
    ],
], columns=[
    "metric", "observed", "mc_mean",
    "mc_sd", "upper_tail_p"
])

global_test.to_csv(
    ROOT / "local_multiscale_point_model_global_rho_test_q7.csv",
    index=False, encoding="utf-8-sig"
)

rho_rows = []

for k in range(MAX_LAG):
    sims = rhos[:, k]

    rho_rows.append({
        "lag": k + 1,
        "Delta_normalized": float(delta_norm[k]),
        "observed_rho": float(obs_rho[k]),
        "short_interval_theory_rho": float(theory_rho[k]),
        "local_multiscale_rho_mean": float(meanrho[k]),
        "local_multiscale_rho_sd": float(sdrho[k]),
        "observed_z_vs_model": float(obs_z[k]),
        "two_sided_p_observed_rho": two_sided_p(
            obs_rho[k], sims
        ),
    })

rho_table = pd.DataFrame(rho_rows)

rho_table.to_csv(
    ROOT / "local_multiscale_point_model_rho_q7.csv",
    index=False, encoding="utf-8-sig"
)

# ------------------------------------------------------------
# Triple sign patterns, still only as diagnostics.
# ------------------------------------------------------------

sg_obs = z_obs >= 0

obs_codes = (
    sg_obs[:-2].astype(np.int8) * 4
    + sg_obs[1:-1].astype(np.int8) * 2
    + sg_obs[2:].astype(np.int8)
)

obs_patterns = np.bincount(
    obs_codes, minlength=8
)

pattern_rows = []

for k in range(8):
    sims = pattern_counts[:, k]

    pattern_rows.append({
        "sign_pattern": f"{k:03b}",
        "observed_count": int(obs_patterns[k]),
        "mc_mean": float(sims.mean()),
        "mc_sd": float(sims.std(ddof=1)),
        "two_sided_p": two_sided_p(
            obs_patterns[k], sims
        ),
    })

patterns = pd.DataFrame(pattern_rows)

patterns.to_csv(
    ROOT / "local_multiscale_point_model_sign_triples_q7.csv",
    index=False, encoding="utf-8-sig"
)

print(summary.to_string(index=False))
print()
print(global_test.to_string(index=False))
print()
print(rho_table.to_string(index=False))
print()
print(patterns.to_string(index=False))
