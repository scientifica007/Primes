#!/usr/bin/env python3
"""
Point-level full rho(k) kernel experiment.

Goal
----
Before adding explicit Hardy--Littlewood k=3 information, test whether a
point-level model using pair structure only can reproduce the full interval
correlation kernel rho(k).

Stages
------
1. Fixed wheel through q0=7.
2. Random forbidden residue class for each prime q0 < p <= Q=7919.
3. Diagnose how independent Bernoulli thinning attenuates covariance.
4. Sweep block-quota thinning:
      block size 1 = Bernoulli thinning;
      one huge block = almost deterministic quota thinning.
5. Build a "pair-kernel quota" bridge model:
      theta_i*S_i + sqrt(lambda*V_i)*G_i,
   where G is Gaussian with a covariance matrix fixed by the short-interval
   pair theory.  lambda is chosen ONLY from the missing marginal variance,
   not by fitting individual lags.
6. Test rho(k), R, runs, and triple sign patterns.

This is deliberately a pair-only experiment.  It does not add k=3 terms.

Repository input
----------------
first_1000_prime_square_interval_index.csv
"""

import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.special import expi

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "first_1000_prime_square_interval_index.csv"

BASE_Q = 7
Q = 7919
OUTER = 80
MAX_LAG = 20
MAX_COV_LAG = 50
BLOCK_INNER = 20
PAIR_INNER = 50
SEED = 20260901

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
z_obs = (obs - mu) / np.sqrt(base_var)

tail = primes[(primes > BASE_Q) & (primes <= Q)]
tail_product = float(np.prod(1 - 1 / tail.astype(float)))
theta = mu / (baseC * tail_product)

rng = np.random.default_rng(SEED)
S_mat = np.empty((OUTER, M), dtype=np.int32)
for outer in range(OUTER):
    mask = base.copy()
    for rv in tail:
        r = int(rv)
        forbidden = int(rng.integers(0, r))
        offset = (forbidden - start) % r
        mask[offset:N:r] = False
    S_mat[outer] = np.add.reduceat(mask, starts).astype(np.int32)

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
    return cov / math.sqrt(V(b - a, Xloc) * V(d - c, Xloc))


theory_rho = np.array([
    np.mean([interval_pair_rho(i, i + k) for i in range(M - k)])
    for k in range(1, MAX_LAG + 1)
])

obs_rho = lag_corrs(z_obs, MAX_LAG)[0]
delta_shape = np.array([Delta(k) for k in range(1, MAX_LAG + 1)])
delta_norm = delta_shape / delta_shape[0]


def fit_delta(rho):
    x = -delta_shape
    amplitude = float(np.dot(x, rho) / np.dot(x, x))
    pred = amplitude * x
    rmse = float(np.sqrt(np.mean((rho - pred) ** 2)))
    normalized = rho / rho[0]
    shape_rmse = float(np.sqrt(np.mean((normalized - delta_norm) ** 2)))
    return amplitude, rmse, shape_rmse


S_mean = S_mat.mean(axis=0)
S_sd = S_mat.std(axis=0, ddof=1)
Z_raw = (S_mat - S_mean) / S_sd
raw_rho = lag_corrs(Z_raw, MAX_LAG)

blocks = [1,2,4,8,16,32,64,128,256,512,1024,2048,4096,8192,32768]
sweep_rows = []
kernel_rows = []
rng = np.random.default_rng(SEED + 1)

for g in blocks:
    q = S_mat // g
    rem = S_mat % g
    eg = theta * g
    ag = np.floor(eg).astype(np.int64)
    fg = eg - ag
    K = q[None, :, :] * ag[None, None, :]
    K = K + rng.binomial(q[None, :, :], fg[None, None, :], size=(BLOCK_INNER,) + q.shape)
    er = theta[None, :] * rem
    ar = np.floor(er).astype(np.int64)
    fr = er - ar
    K = K + ar[None, :, :] + rng.binomial(1, fr[None, :, :], size=(BLOCK_INNER,) + rem.shape)
    K = K.astype(float).reshape(-1, M)
    Z = (K - mu) / np.sqrt(base_var)
    rhos = lag_corrs(Z, MAX_LAG)
    rho_mean = rhos.mean(axis=0)
    rho_sd = rhos.std(axis=0, ddof=1)
    Rvals = np.sum(Z * Z, axis=1) / (M - 1)
    sg = Z >= 0
    runs = 1 + np.sum(sg[:, 1:] != sg[:, :-1], axis=1)
    amp, delta_rmse, delta_shape_rmse = fit_delta(rho_mean)
    R_obs = float(np.sum(z_obs * z_obs) / (M - 1))
    zR = (Rvals.mean() - R_obs) / Rvals.std(ddof=1)
    kernel_z_rmse = float(np.sqrt(np.mean(((rho_mean - theory_rho) / rho_sd) ** 2)))
    sweep_rows.append({
        "base_q": BASE_Q,
        "block_size_g": g,
        "simulations": len(Rvals),
        "R_observed": R_obs,
        "R_mc_mean": float(Rvals.mean()),
        "R_mc_sd": float(Rvals.std(ddof=1)),
        "R_error_z": float(zR),
        "runs_observed": int(1 + np.sum((z_obs[1:] >= 0) != (z_obs[:-1] >= 0))),
        "runs_mc_mean": float(runs.mean()),
        "runs_mc_sd": float(runs.std(ddof=1)),
        "rho1_mc_mean": float(rho_mean[0]),
        "rho2_mc_mean": float(rho_mean[1]),
        "rho5_mc_mean": float(rho_mean[4]),
        "delta_amplitude": amp,
        "delta_fit_rmse": delta_rmse,
        "delta_normalized_shape_rmse": delta_shape_rmse,
        "kernel_z_rmse_vs_theory": kernel_z_rmse,
        "joint_score_R_kernel": float(math.sqrt(zR * zR + kernel_z_rmse * kernel_z_rmse)),
    })
    for k in range(1, MAX_LAG + 1):
        kernel_rows.append({
            "base_q": BASE_Q,
            "block_size_g": g,
            "lag": k,
            "rho_mc_mean": float(rho_mean[k-1]),
            "rho_mc_sd": float(rho_sd[k-1]),
            "short_interval_theory_rho": float(theory_rho[k-1]),
            "Delta_normalized": float(delta_norm[k-1]),
        })

pd.DataFrame(sweep_rows).to_csv(ROOT / "block_quota_thinning_sweep_q7.csv", index=False, encoding="utf-8-sig")
pd.DataFrame(kernel_rows).to_csv(ROOT / "block_quota_full_rho_kernels_q7.csv", index=False, encoding="utf-8-sig")

Z_struct = (theta[None, :] * S_mat - mu) / np.sqrt(base_var)
R_struct = float(np.mean(np.sum(Z_struct * Z_struct, axis=1) / (M - 1)))
rho_struct = lag_corrs(Z_struct, MAX_LAG).mean(axis=0)
R_obs = float(np.sum(z_obs * z_obs) / (M - 1))
lambda_missing = max(0.0, R_obs - R_struct)

R_pair = np.eye(M)
for lag in range(1, MAX_COV_LAG + 1):
    for i in range(M - lag):
        rr = interval_pair_rho(i, i + lag)
        R_pair[i, i + lag] = rr
        R_pair[i + lag, i] = rr

min_eigenvalue = float(np.linalg.eigvalsh(R_pair).min())
if min_eigenvalue <= 0:
    raise RuntimeError(f"Pair covariance matrix is not positive definite: {min_eigenvalue}")
chol = np.linalg.cholesky(R_pair)

rng = np.random.default_rng(SEED + 2)
all_rho = []
all_R = []
all_runs = []
all_mono = []
pattern_counts = []

for S in S_mat:
    G = rng.standard_normal((PAIR_INNER, M)) @ chol.T
    target = theta[None, :] * S[None, :] + np.sqrt(base_var[None, :] * lambda_missing) * G
    lo = np.floor(target)
    frac = target - lo
    K = lo + rng.binomial(1, np.clip(frac, 0, 1), size=target.shape)
    K = np.clip(K, 0, S[None, :]).astype(float)
    Z = (K - mu) / np.sqrt(base_var)
    all_rho.append(lag_corrs(Z, MAX_LAG))
    all_R.extend((np.sum(Z * Z, axis=1) / (M - 1)).tolist())
    sg = Z >= 0
    all_runs.extend((1 + np.sum(sg[:, 1:] != sg[:, :-1], axis=1)).tolist())
    mono = np.sum(((Z[:, :-2] < Z[:, 1:-1]) & (Z[:, 1:-1] < Z[:, 2:])) | ((Z[:, :-2] > Z[:, 1:-1]) & (Z[:, 1:-1] > Z[:, 2:])), axis=1)
    all_mono.extend(mono.tolist())
    codes = sg[:, :-2].astype(np.int8) * 4 + sg[:, 1:-1].astype(np.int8) * 2 + sg[:, 2:].astype(np.int8)
    for row in codes:
        pattern_counts.append(np.bincount(row, minlength=8))

all_rho = np.vstack(all_rho)
all_R = np.asarray(all_R)
all_runs = np.asarray(all_runs)
all_mono = np.asarray(all_mono)
pattern_counts = np.asarray(pattern_counts)
rho_pairkernel = all_rho.mean(axis=0)
rho_pairkernel_sd = all_rho.std(axis=0, ddof=1)
amp, _, delta_shape_rmse = fit_delta(rho_pairkernel)

trip = np.column_stack([z_obs[:-2], z_obs[1:-1], z_obs[2:]])
obs_mono = int(np.sum(((trip[:, 0] < trip[:, 1]) & (trip[:, 1] < trip[:, 2])) | ((trip[:, 0] > trip[:, 1]) & (trip[:, 1] > trip[:, 2]))))
sg_obs = z_obs >= 0
obs_codes = sg_obs[:-2].astype(np.int8) * 4 + sg_obs[1:-1].astype(np.int8) * 2 + sg_obs[2:].astype(np.int8)
obs_patterns = np.bincount(obs_codes, minlength=8)
obs_runs = int(1 + np.sum(sg_obs[1:] != sg_obs[:-1]))


def two_sided_p(obs_value, simulations):
    simulations = np.asarray(simulations)
    center = simulations.mean()
    return float((1 + np.sum(np.abs(simulations - center) >= abs(obs_value - center))) / (len(simulations) + 1))

summary = pd.DataFrame([
    ["base_q", BASE_Q],
    ["lambda_missing_variance", lambda_missing],
    ["R_struct_before_latent", R_struct],
    ["R_observed", R_obs],
    ["R_pairkernel_mc_mean", float(all_R.mean())],
    ["R_pairkernel_mc_sd", float(all_R.std(ddof=1))],
    ["R_pairkernel_p", two_sided_p(R_obs, all_R)],
    ["rho1_observed", float(obs_rho[0])],
    ["rho1_struct", float(rho_struct[0])],
    ["rho1_pairkernel", float(rho_pairkernel[0])],
    ["rho1_short_interval_theory", float(theory_rho[0])],
    ["delta_amplitude", amp],
    ["delta_normalized_shape_rmse", delta_shape_rmse],
    ["runs_observed", obs_runs],
    ["runs_pairkernel_mc_mean", float(all_runs.mean())],
    ["runs_pairkernel_mc_sd", float(all_runs.std(ddof=1))],
    ["runs_pairkernel_p", two_sided_p(obs_runs, all_runs)],
    ["monotone_triples_observed", obs_mono],
    ["monotone_triples_mc_mean", float(all_mono.mean())],
    ["monotone_triples_mc_sd", float(all_mono.std(ddof=1))],
    ["monotone_triples_p", two_sided_p(obs_mono, all_mono)],
    ["pair_covariance_min_eigenvalue", min_eigenvalue],
], columns=["metric", "value"])
summary.to_csv(ROOT / "pair_kernel_point_model_summary_q7.csv", index=False, encoding="utf-8-sig")

obs_z = (obs_rho - rho_pairkernel) / rho_pairkernel_sd
sim_z = (all_rho - rho_pairkernel[None, :]) / rho_pairkernel_sd[None, :]
T_obs = float(np.sqrt(np.mean(obs_z ** 2)))
M_obs = float(np.max(np.abs(obs_z)))
T_sim = np.sqrt(np.mean(sim_z ** 2, axis=1))
M_sim = np.max(np.abs(sim_z), axis=1)

global_test = pd.DataFrame([
    ["rho_kernel_rms_standardized_deviation", T_obs, float(T_sim.mean()), float(T_sim.std(ddof=1)), float((1 + np.sum(T_sim >= T_obs)) / (len(T_sim) + 1))],
    ["rho_kernel_max_abs_standardized_deviation", M_obs, float(M_sim.mean()), float(M_sim.std(ddof=1)), float((1 + np.sum(M_sim >= M_obs)) / (len(M_sim) + 1))],
], columns=["metric", "observed", "mc_mean", "mc_sd", "upper_tail_p"])
global_test.to_csv(ROOT / "pair_kernel_point_model_global_rho_test_q7.csv", index=False, encoding="utf-8-sig")

rho_rows = []
for k in range(1, MAX_LAG + 1):
    sims = all_rho[:, k-1]
    rho_rows.append({
        "lag": k,
        "Delta_normalized": float(delta_norm[k-1]),
        "observed_rho": float(obs_rho[k-1]),
        "short_interval_theory_rho": float(theory_rho[k-1]),
        "structural_random_sieve_rho": float(rho_struct[k-1]),
        "pairkernel_point_model_rho_mean": float(rho_pairkernel[k-1]),
        "pairkernel_point_model_rho_sd": float(rho_pairkernel_sd[k-1]),
        "observed_z_vs_pairkernel": float(obs_z[k-1]),
        "two_sided_p_observed_rho": two_sided_p(obs_rho[k-1], sims),
    })
pd.DataFrame(rho_rows).to_csv(ROOT / "pair_kernel_point_model_rho_q7.csv", index=False, encoding="utf-8-sig")

pattern_rows = []
for k in range(8):
    sims = pattern_counts[:, k]
    pattern_rows.append({
        "sign_pattern": f"{k:03b}",
        "observed_count": int(obs_patterns[k]),
        "pairkernel_mean_count": float(sims.mean()),
        "pairkernel_sd": float(sims.std(ddof=1)),
        "two_sided_p": two_sided_p(obs_patterns[k], sims),
    })
pd.DataFrame(pattern_rows).to_csv(ROOT / "pair_kernel_point_model_sign_triples_q7.csv", index=False, encoding="utf-8-sig")

print(summary.to_string(index=False))
print(global_test.to_string(index=False))
