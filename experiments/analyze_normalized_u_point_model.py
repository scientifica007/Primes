#!/usr/bin/env python3
"""
Frozen normalized-u spectrum inside the actual random-sieve point model.
Primary experiment
------------------
1. Use the p>=4001 prime-square sample (450 intervals).
2. Fixed wheel q0=7.
3. Random forbidden residue class for every prime 7<q<=7919.
4. Split each physical prime-square interval into exactly 64 local-u cells.
5. Load the already-frozen normalized-u reference spectrum from
   normalized_u_reference_spectrum.csv.
6. Apply the spectrum to random-sieve survivors with one scalar amplitude.
7. Calibrate that amplitude using ONLY R on 12 training raw-sieve realizations.
   Calibration includes the expected variance of randomized integer rounding.
8. Test on 28 held-out raw-sieve realizations x 24 fields = 672 simulations.
9. Report R, rho(k), runs, and triple-sign diagnostics.
No target lag is used to fit the spectrum or amplitude.
No k=3 information is used.
"""
import math
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.special import expi
ROOT = Path(__file__).resolve().parent
INTERVAL_FILE = ROOT / "first_1000_prime_square_interval_index.csv"
SPECTRUM_FILE = ROOT / "normalized_u_reference_spectrum.csv"
BASE_Q = 7
Q = 7919
CPI = 64
MAX_LAG = 20
RAW_REALIZATIONS = 40
TRAIN_RAW = 12
TRAIN_FIELDS = 25
TEST_FIELDS = 24
RAW_SEEDS = [20261010, 20261011]
CALIBRATION_SEED = 20260990
TEST_SEED = 20260970
def primes_upto(n):
    z = np.ones(n + 1, dtype=bool)
    z[:2] = False
    for q in range(2, math.isqrt(n) + 1):
        if z[q]:
            z[q*q:n+1:q] = False
    return np.flatnonzero(z)
def lag_corrs(Z, max_lag=MAX_LAG):
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
            np.sum(xc * xc, axis=1)
            * np.sum(yc * yc, axis=1)
        )
    return out
def V(h, X):
    if h <= 0:
        return 0.0
    return h * math.log(X / h)
def interval_pair_rho(A, B, i, j):
    a, b = float(A[i]), float(B[i])
    c, d = float(A[j]), float(B[j])
    centers = 0.5 * (A + B)
    Xloc = math.sqrt(float(centers[i]) * float(centers[j]))
    cov = 0.5 * (
        V(d - a, Xloc)
        + V(c - b, Xloc)
        - V(c - a, Xloc)
        - V(d - b, Xloc)
    )
    return cov / math.sqrt(
        V(b - a, Xloc)
        * V(d - c, Xloc)
    )
def Delta(t):
    t = float(t)
    return 0.5 * (
        (t + 1) * math.log(t + 1)
        - 2 * t * math.log(t)
        + (
            0.0
            if t == 1
            else (t - 1) * math.log(t - 1)
        )
    )
def two_sided_p(observed, simulations):
    simulations = np.asarray(simulations)
    center = simulations.mean()
    return float(
        (
            1
            + np.sum(
                np.abs(simulations - center)
                >= abs(observed - center)
            )
        )
        / (len(simulations) + 1)
    )
df = pd.read_csv(INTERVAL_FILE)
s = df[df["p_n"] >= 4001].copy().reset_index(drop=True)
A = s["p_n^2"].to_numpy(np.int64)
B = s["p_التالي^2"].to_numpy(np.int64)
obs = s["عدد_الأوليات_بين_p2_ومربع_التالي"].to_numpy(float)
M = len(A)
start = int(A[0])
end = int(B[-1])
N = end - start
shape = (
    expi(np.log(B.astype(float)))
    - expi(np.log(A.astype(float)))
)
alpha = obs.sum() / shape.sum()
mu = alpha * shape
primes = primes_upto(Q)
base = np.ones(N, dtype=bool)
for q0 in primes[primes <= BASE_Q]:
    q = int(q0)
    base[(-start) % q::q] = False
for x in np.unique(np.concatenate([A, B])):
    j = int(x) - start
    if 0 <= j < N:
        base[j] = False
interval_starts = A - start
baseC = np.add.reduceat(
    base,
    interval_starts
).astype(np.int64)
base_var = (
    mu
    * (1 - mu / baseC)
)
sd_base = np.sqrt(base_var)
z_obs = (
    obs - mu
) / sd_base
R_obs = float(
    np.sum(z_obs * z_obs)
    / (M - 1)
)
obs_runs = int(
    1 + np.sum(
        (z_obs[1:] >= 0)
        != (z_obs[:-1] >= 0)
    )
)
obs_rho = lag_corrs(
    z_obs
)[0]
tail = primes[
    (primes > BASE_Q)
    & (primes <= Q)
]
tail_product = float(
    np.prod(
        1
        - 1 / tail.astype(float)
    )
)
theta = (
    mu
    / (baseC * tail_product)
)
cell_starts = []
cell_intervals = []
for i, (a, b) in enumerate(zip(A, B)):
    H = int(b - a)
    bounds = (
        a
        + (
            np.arange(
                CPI + 1,
                dtype=np.int64
            )
            * H
        )
        // CPI
    )
    if not np.all(
        np.diff(bounds) > 0
    ):
        raise RuntimeError(
            "An interval is too short for the requested u resolution."
        )
    cell_starts.extend(
        (bounds[:-1] - start).tolist()
    )
    cell_intervals.extend(
        [i] * CPI
    )
cell_starts = np.asarray(
    cell_starts,
    dtype=np.int64
)
cell_intervals = np.asarray(
    cell_intervals,
    dtype=np.int16
)
NCELL = len(cell_starts)
u_index = np.arange(
    NCELL,
    dtype=np.int64
)
def interval_sum_cells(values):
    return np.add.reduceat(
        values,
        np.arange(
            0,
            NCELL,
            CPI
        )
    )
parts = []
per_seed = (
    RAW_REALIZATIONS
    // len(RAW_SEEDS)
)
for seed in RAW_SEEDS:
    rng = np.random.default_rng(seed)
    part = np.empty(
        (per_seed, NCELL),
        dtype=np.uint16
    )
    for j in range(per_seed):
        mask = base.copy()
        for rv in tail:
            r = int(rv)
            forbidden = int(
                rng.integers(0, r)
            )
            offset = (
                forbidden - start
            ) % r
            mask[
                offset:N:r
            ] = False
        part[j] = np.add.reduceat(
            mask,
            cell_starts
        ).astype(np.uint16)
    parts.append(part)
S_cells = np.vstack(
    parts
).astype(np.float64)
TRAIN_IDX = np.arange(
    0,
    TRAIN_RAW
)
TEST_IDX = np.arange(
    TRAIN_RAW,
    RAW_REALIZATIONS
)
spectrum = pd.read_csv(
    SPECTRUM_FILE
)
active = spectrum[
    spectrum["coefficient_v"]
    > 1e-12
].copy()
active_scales = active[
    "support_cells"
].to_numpy(np.int64)
active_v = active[
    "coefficient_v"
].to_numpy(float)
def normalized_u_field(rng):
    field = np.zeros(
        NCELL,
        dtype=float
    )
    for L, var_weight in zip(
        active_scales,
        active_v
    ):
        shifts = np.array(
            sorted(set([
                0,
                int(round(L / 4)),
                int(round(L / 2)),
                int(round(3 * L / 4)),
            ])),
            dtype=np.int64
        )
        shift = int(
            rng.choice(shifts)
        )
        u = u_index + shift
        block = u // L
        pos = u % L
        sign = np.where(
            pos < L / 2,
            1.0,
            -1.0
        )
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
def expected_training_R(amplitude):
    rng = np.random.default_rng(
        CALIBRATION_SEED
    )
    values = []
    for ridx in TRAIN_IDX:
        S = S_cells[ridx]
        for _ in range(TRAIN_FIELDS):
            field = normalized_u_field(
                rng
            )
            pcell = (
                theta[cell_intervals]
                + amplitude * field
            )
            target = pcell * S
            frac = (
                target
                - np.floor(target)
            )
            Y_cont = interval_sum_cells(
                target
            )
            Z_cont = (
                Y_cont - mu
            ) / sd_base
            rounding_variance = interval_sum_cells(
                frac * (1 - frac)
            )
            expected_square = (
                Z_cont * Z_cont
                + rounding_variance
                / base_var
            )
            values.append(
                np.sum(
                    expected_square
                )
                / (M - 1)
            )
    return float(
        np.mean(values)
    )
lo = 0.0
hi = 0.004
while expected_training_R(hi) < R_obs:
    hi *= 1.5
for _ in range(18):
    mid = (lo + hi) / 2
    if (
        expected_training_R(mid)
        < R_obs
    ):
        lo = mid
    else:
        hi = mid
amplitude = (
    lo + hi
) / 2
rng = np.random.default_rng(
    TEST_SEED
)
Rvals = []
runs = []
rhos = []
monotone = []
patterns = []
clip_count = 0
total_probabilities = 0
for ridx in TEST_IDX:
    S = S_cells[ridx]
    for _ in range(TEST_FIELDS):
        field = normalized_u_field(
            rng
        )
        pcell = (
            theta[cell_intervals]
            + amplitude * field
        )
        clip_count += int(
            np.sum(
                (pcell < 0)
                | (pcell > 1)
            )
        )
        total_probabilities += (
            len(pcell)
        )
        pcell = np.clip(
            pcell,
            0,
            1
        )
        target = pcell * S
        lo_int = np.floor(target)
        K = (
            lo_int
            + rng.binomial(
                1,
                target - lo_int
            )
        )
        Y = interval_sum_cells(K)
        Z = (
            Y - mu
        ) / sd_base
        Rvals.append(
            np.sum(Z * Z)
            / (M - 1)
        )
        sg = Z >= 0
        runs.append(
            1 + np.sum(
                sg[1:] != sg[:-1]
            )
        )
        rhos.append(
            lag_corrs(Z)[0]
        )
        monotone.append(
            np.sum(
                (
                    (Z[:-2] < Z[1:-1])
                    & (Z[1:-1] < Z[2:])
                )
                |
                (
                    (Z[:-2] > Z[1:-1])
                    & (Z[1:-1] > Z[2:])
                )
            )
        )
        codes = (
            sg[:-2].astype(np.int8) * 4
            + sg[1:-1].astype(np.int8) * 2
            + sg[2:].astype(np.int8)
        )
        patterns.append(
            np.bincount(
                codes,
                minlength=8
            )
        )
Rvals = np.asarray(Rvals)
runs = np.asarray(runs)
rhos = np.asarray(rhos)
monotone = np.asarray(monotone)
patterns = np.asarray(patterns)
meanrho = rhos.mean(axis=0)
sdrho = rhos.std(
    axis=0,
    ddof=1
)
theory_rho = np.array([
    np.mean([
        interval_pair_rho(
            A, B,
            i, i + k
        )
        for i
        in range(M - k)
    ])
    for k
    in range(
        1,
        MAX_LAG + 1
    )
])
delta = np.array([
    Delta(k)
    for k
    in range(
        1,
        MAX_LAG + 1
    )
])
delta_norm = (
    delta / delta[0]
)
kernel_rmse = float(
    np.sqrt(np.mean(
        (
            meanrho
            - theory_rho
        ) ** 2
    ))
)
shape_rmse = float(
    np.sqrt(np.mean(
        (
            meanrho
            / meanrho[0]
            - delta_norm
        ) ** 2
    ))
)
observed_z = (
    obs_rho - meanrho
) / sdrho
T_obs = float(
    np.sqrt(
        np.mean(
            observed_z ** 2
        )
    )
)
M_obs = float(
    np.max(
        np.abs(
            observed_z
        )
    )
)
simulation_z = (
    rhos
    - meanrho[None, :]
) / sdrho[None, :]
T_sim = np.sqrt(
    np.mean(
        simulation_z ** 2,
        axis=1
    )
)
M_sim = np.max(
    np.abs(simulation_z),
    axis=1
)
global_RMS_p = float(
    (
        1
        + np.sum(
            T_sim >= T_obs
        )
    )
    / (len(T_sim) + 1)
)
global_max_p = float(
    (
        1
        + np.sum(
            M_sim >= M_obs
        )
    )
    / (len(M_sim) + 1)
)
trip = np.column_stack([
    z_obs[:-2],
    z_obs[1:-1],
    z_obs[2:]
])
obs_monotone = int(
    np.sum(
        (
            (trip[:, 0] < trip[:, 1])
            & (trip[:, 1] < trip[:, 2])
        )
        |
        (
            (trip[:, 0] > trip[:, 1])
            & (trip[:, 1] > trip[:, 2])
        )
    )
)
summary = pd.DataFrame([
    [
        "reference_spectrum",
        "uniform normalized-u spectrum frozen from X~2e7"
    ],
    [
        "training_raw_sieves",
        TRAIN_RAW
    ],
    [
        "heldout_raw_sieves",
        len(TEST_IDX)
    ],
    [
        "heldout_simulations",
        len(Rvals)
    ],
    [
        "amplitude_from_training_R_after_rounding",
        amplitude
    ],
    [
        "R_observed",
        R_obs
    ],
    [
        "R_mc_mean",
        float(Rvals.mean())
    ],
    [
        "R_mc_sd",
        float(Rvals.std(ddof=1))
    ],
    [
        "R_two_sided_p",
        two_sided_p(
            R_obs, Rvals
        )
    ],
    [
        "rho1_observed",
        float(obs_rho[0])
    ],
    [
        "rho1_mc_mean",
        float(meanrho[0])
    ],
    [
        "rho1_mc_sd",
        float(sdrho[0])
    ],
    [
        "rho1_two_sided_p",
        two_sided_p(
            obs_rho[0],
            rhos[:, 0]
        )
    ],
    [
        "rho2_observed",
        float(obs_rho[1])
    ],
    [
        "rho2_mc_mean",
        float(meanrho[1])
    ],
    [
        "kernel_rmse_vs_short_interval_theory",
        kernel_rmse
    ],
    [
        "normalized_Delta_shape_RMSE",
        shape_rmse
    ],
    [
        "rho_kernel_global_RMS_p",
        global_RMS_p
    ],
    [
        "rho_kernel_global_max_p",
        global_max_p
    ],
    [
        "runs_observed",
        obs_runs
    ],
    [
        "runs_mc_mean",
        float(runs.mean())
    ],
    [
        "runs_mc_sd",
        float(runs.std(ddof=1))
    ],
    [
        "runs_two_sided_p",
        two_sided_p(
            obs_runs, runs
        )
    ],
    [
        "monotone_triples_observed",
        obs_monotone
    ],
    [
        "monotone_triples_mc_mean",
        float(monotone.mean())
    ],
    [
        "monotone_triples_mc_sd",
        float(monotone.std(ddof=1))
    ],
    [
        "monotone_triples_two_sided_p",
        two_sided_p(
            obs_monotone,
            monotone
        )
    ],
    [
        "clip_fraction",
        clip_count
        / total_probabilities
    ],
], columns=[
    "metric", "value"
])
summary.to_csv(
    ROOT /
    "normalized_u_point_model_summary_q7.csv",
    index=False,
    encoding="utf-8-sig"
)
rho_rows = []
for k in range(MAX_LAG):
    rho_rows.append({
        "lag":
            k + 1,
        "observed_rho":
            float(obs_rho[k]),
        "short_interval_theory_rho":
            float(theory_rho[k]),
        "normalized_u_point_model_rho_mean":
            float(meanrho[k]),
        "normalized_u_point_model_rho_sd":
            float(sdrho[k]),
        "observed_z_vs_model":
            float(observed_z[k]),
        "two_sided_p_observed_rho":
            two_sided_p(
                obs_rho[k],
                rhos[:, k]
            ),
    })
pd.DataFrame(
    rho_rows
).to_csv(
    ROOT /
    "normalized_u_point_model_rho_q7.csv",
    index=False,
    encoding="utf-8-sig"
)
pd.DataFrame([
    [
        "rho_kernel_RMS_standardized",
        T_obs,
        global_RMS_p
    ],
    [
        "rho_kernel_max_abs_standardized",
        M_obs,
        global_max_p
    ],
], columns=[
    "metric",
    "observed_statistic",
    "upper_tail_p"
]).to_csv(
    ROOT /
    "normalized_u_point_model_global_rho_test_q7.csv",
    index=False,
    encoding="utf-8-sig"
)
sg_obs = z_obs >= 0
obs_codes = (
    sg_obs[:-2].astype(np.int8) * 4
    + sg_obs[1:-1].astype(np.int8) * 2
    + sg_obs[2:].astype(np.int8)
)
obs_patterns = np.bincount(
    obs_codes,
    minlength=8
)
pattern_rows = []
for k in range(8):
    sims = patterns[:, k]
    pattern_rows.append({
        "sign_pattern":
            f"{k:03b}",
        "observed_count":
            int(obs_patterns[k]),
        "mc_mean":
            float(sims.mean()),
        "mc_sd":
            float(sims.std(ddof=1)),
        "two_sided_p":
            two_sided_p(
                obs_patterns[k],
                sims
            ),
    })
pd.DataFrame(
    pattern_rows
).to_csv(
    ROOT /
    "normalized_u_point_model_sign_triples_q7.csv",
    index=False,
    encoding="utf-8-sig"
)
print(summary.to_string(index=False))
