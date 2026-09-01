#!/usr/bin/env python3
"""
Theory-aligned comparison with:
  - Sun-Kai Leung (2026): joint distribution of primes in multiple short intervals
  - Hardy--Littlewood / short-interval covariance viewpoint
  - the project's prime-square interval pair-Gaussian model

Outputs:
  leung2026_weighted_interval_summary.csv
  leung2026_lag_correlation_comparison.csv
  leung2026_triple_ordering_comparison.csv
  prime_square_pair_gaussian_higher_order_diagnostics.csv
  prime_square_pair_gaussian_sign_triples.csv
  theory_aligned_comparison_summary.csv

The experiment deliberately checks whether pair covariance plus a multivariate
Gaussian approximation already explains local "triple" behavior before adding
an explicit k=3 Hardy--Littlewood layer.

Requirements:
  scipy, numpy, pandas

Repository inputs:
  first_1000_prime_square_interval_index.csv
  primes_up_to_1001st_prime_square.csv
  negative_correlation_details.csv
"""

import math
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent

PRIME_FILE = ROOT / "primes_up_to_1001st_prime_square.csv"
primes = pd.read_csv(PRIME_FILE, usecols=["prime"])["prime"].to_numpy(np.int64)

MAX_Y = 62_800_000
powers = []
weights = []
for p in primes:
    p = int(p)
    if p > MAX_Y:
        break
    v = p
    lp = math.log(p)
    while v <= MAX_Y:
        powers.append(v)
        weights.append(lp)
        if v > MAX_Y // p:
            break
        v *= p

order = np.argsort(powers)
pp = np.asarray(powers, dtype=np.int64)[order]
ww = np.asarray(weights, dtype=np.float64)[order]
cw = np.concatenate([[0.0], np.cumsum(ww)])

def psi(y):
    y = np.asarray(y, dtype=np.float64)
    idx = np.searchsorted(pp, y, side="right")
    return cw[idx]

gamma = 0.5772156649015328606
c0 = 1 - gamma - math.log(2 * math.pi)

def Delta(t):
    t = float(t)
    a = (t + 1) * math.log(t + 1)
    b = 2 * t * math.log(t)
    c = 0.0 if abs(t - 1) < 1e-15 else (t - 1) * math.log(t - 1)
    return 0.5 * (a - b + c)

deltas = [0.001, 0.002, 0.004, 0.008]
r = 4
samples = 8000
x_lo = 20_000_000.0
x_hi = 55_000_000.0
summary_rows = []
lag_rows = []
ordering_rows = []
rng = np.random.default_rng(20260901)

for delta in deltas:
    u = np.linspace(math.log(x_lo), math.log(x_hi), samples)
    x = np.exp(u)
    Z = np.empty((samples, r), dtype=float)
    Vdelta = delta * math.log(1 / delta) + c0 * delta

    for j in range(r):
        t = j
        lo = (1 + (t - 0.5) * delta) * x
        hi = (1 + (t + 0.5) * delta) * x
        weighted = psi(hi) - psi(lo)
        E = (weighted - delta * x) / np.sqrt(x)
        Z[:, j] = E / math.sqrt(Vdelta)

    emp_corr = np.corrcoef(Z, rowvar=False)
    emp_var = Z.var(axis=0, ddof=1)
    emp_mean = Z.mean(axis=0)
    centered = (Z - Z.mean(axis=0)) / Z.std(axis=0, ddof=0)
    emp_skew = np.mean(centered ** 3, axis=0)
    emp_kurt = np.mean(centered ** 4, axis=0)

    C_lead = np.eye(r)
    C_sec = np.eye(r)
    Dlead = math.log(1 / delta)
    Dsec = Dlead + c0
    for j in range(r):
        for k in range(j + 1, r):
            d = abs(j - k)
            C_lead[j, k] = C_lead[k, j] = -Delta(d) / Dlead
            C_sec[j, k] = C_sec[k, j] = -Delta(d) / Dsec

    sims = 400_000
    G = rng.multivariate_normal(np.zeros(r), C_sec, size=sims)
    emp_runs = 1 + np.sum((Z[:, 1:] >= 0) != (Z[:, :-1] >= 0), axis=1)
    gauss_runs = 1 + np.sum((G[:, 1:] >= 0) != (G[:, :-1] >= 0), axis=1)

    perms = [(0,1,2),(0,2,1),(1,0,2),(1,2,0),(2,0,1),(2,1,0)]
    oe = np.argsort(Z[:, :3], axis=1)
    og = np.argsort(G[:, :3], axis=1)
    for perm in perms:
        pvec = np.asarray(perm)
        fe = np.mean(np.all(oe == pvec, axis=1))
        fg = np.mean(np.all(og == pvec, axis=1))
        ordering_rows.append({
            "delta": delta,
            "ordering": "".join(map(str, perm)),
            "empirical_frequency": float(fe),
            "gaussian_theory_frequency": float(fg),
            "difference": float(fe - fg),
        })

    summary_rows.append({
        "delta": delta,
        "samples_log_weighted": samples,
        "variance_prediction_Vdelta": Vdelta,
        "empirical_mean_z_avg": float(emp_mean.mean()),
        "empirical_variance_z_avg": float(emp_var.mean()),
        "empirical_skewness_avg": float(emp_skew.mean()),
        "empirical_kurtosis_avg": float(emp_kurt.mean()),
        "empirical_runs4_mean": float(emp_runs.mean()),
        "gaussian_theory_runs4_mean": float(gauss_runs.mean()),
        "runs4_difference": float(emp_runs.mean() - gauss_runs.mean()),
        "corr_matrix_rmse_leading": float(np.sqrt(np.mean((emp_corr - C_lead) ** 2))),
        "corr_matrix_rmse_secondary": float(np.sqrt(np.mean((emp_corr - C_sec) ** 2))),
    })

    for d in range(1, r):
        vals = [emp_corr[j, j + d] for j in range(r - d)]
        lag_rows.append({
            "delta": delta,
            "lag_in_interval_units": d,
            "empirical_corr_mean": float(np.mean(vals)),
            "empirical_corr_min": float(np.min(vals)),
            "empirical_corr_max": float(np.max(vals)),
            "leung_leading_prediction": float(-Delta(d) / Dlead),
            "leung_secondary_prediction": float(-Delta(d) / Dsec),
            "Delta_d": float(Delta(d)),
        })

summary = pd.DataFrame(summary_rows)
lags = pd.DataFrame(lag_rows)
orders = pd.DataFrame(ordering_rows)
summary.to_csv(ROOT / "leung2026_weighted_interval_summary.csv", index=False, encoding="utf-8-sig")
lags.to_csv(ROOT / "leung2026_lag_correlation_comparison.csv", index=False, encoding="utf-8-sig")
orders.to_csv(ROOT / "leung2026_triple_ordering_comparison.csv", index=False, encoding="utf-8-sig")

detail = pd.read_csv(ROOT / "negative_correlation_details.csv")
z_obs = detail["z_A"].to_numpy(float)
p = detail["p_n"].to_numpy(np.int64)
p_next = detail["p_next"].to_numpy(np.int64)
A = p * p
B = p_next * p_next
centers = 0.5 * (A + B)

def V_short(h, X):
    if h <= 0:
        return 0.0
    return h * math.log(X / h)

def pair_rho(i, j):
    a, b = float(A[i]), float(B[i])
    c, d = float(A[j]), float(B[j])
    Xloc = math.sqrt(float(centers[i]) * float(centers[j]))
    cov = 0.5 * (V_short(d-a, Xloc) + V_short(c-b, Xloc) - V_short(c-a, Xloc) - V_short(d-b, Xloc))
    return cov / math.sqrt(V_short(b-a, Xloc) * V_short(d-c, Xloc))

N = len(z_obs)
max_lag = 50
R = np.eye(N)
for lag in range(1, max_lag + 1):
    for i in range(N - lag):
        rr = pair_rho(i, i + lag)
        R[i, i + lag] = rr
        R[i + lag, i] = rr

eigmin = float(np.linalg.eigvalsh(R).min())
if eigmin <= 0:
    raise RuntimeError(f"Pair covariance matrix is not positive definite: {eigmin}")
chol = np.linalg.cholesky(R)

trip = np.column_stack([z_obs[:-2], z_obs[1:-1], z_obs[2:]])
obs_monotone = int(np.sum(((trip[:,0] < trip[:,1]) & (trip[:,1] < trip[:,2])) | ((trip[:,0] > trip[:,1]) & (trip[:,1] > trip[:,2]))))
obs_turning = len(trip) - obs_monotone

def sign_codes(X):
    return ((X[:,0] >= 0).astype(int) * 4 + (X[:,1] >= 0).astype(int) * 2 + (X[:,2] >= 0).astype(int))

obs_codes = sign_codes(trip)
obs_pattern_counts = np.bincount(obs_codes, minlength=8)
rng = np.random.default_rng(20260901)
REPS = 50_000
BATCH = 250
mono = np.empty(REPS, dtype=np.int16)
pattern_counts = np.empty((REPS, 8), dtype=np.int16)
runs = np.empty(REPS, dtype=np.int16)

for st in range(0, REPS, BATCH):
    m = min(BATCH, REPS - st)
    G = rng.standard_normal((m, N)) @ chol.T
    T0, T1, T2 = G[:, :-2], G[:, 1:-1], G[:, 2:]
    mono[st:st+m] = np.sum(((T0 < T1) & (T1 < T2)) | ((T0 > T1) & (T1 > T2)), axis=1)
    sg = G >= 0
    runs[st:st+m] = 1 + np.sum(sg[:,1:] != sg[:,:-1], axis=1)
    codes = sg[:,:-2].astype(np.int8) * 4 + sg[:,1:-1].astype(np.int8) * 2 + sg[:,2:].astype(np.int8)
    for k in range(8):
        pattern_counts[st:st+m, k] = np.sum(codes == k, axis=1)

def two_sided_p(obs, sims):
    center = sims.mean()
    return (1 + np.sum(np.abs(sims - center) >= abs(obs - center))) / (len(sims) + 1)

higher = pd.DataFrame([
    ["min_eigenvalue_covariance", eigmin, np.nan, np.nan, np.nan],
    ["runs", 256.0, float(runs.mean()), float(runs.std(ddof=1)), two_sided_p(256, runs)],
    ["monotone_triples", float(obs_monotone), float(mono.mean()), float(mono.std(ddof=1)), two_sided_p(obs_monotone, mono)],
    ["turning_triples", float(obs_turning), float((448-mono).mean()), float((448-mono).std(ddof=1)), two_sided_p(obs_turning, 448-mono)],
], columns=["metric","observed","pair_gaussian_mean","pair_gaussian_sd","two_sided_p"])

pattern_rows = []
for k in range(8):
    sims = pattern_counts[:, k]
    pattern_rows.append([f"{k:03b}", int(obs_pattern_counts[k]), float(sims.mean()), float(sims.std(ddof=1)), two_sided_p(obs_pattern_counts[k], sims)])
patterns = pd.DataFrame(pattern_rows, columns=["sign_pattern","observed_count","pair_gaussian_mean","pair_gaussian_sd","two_sided_p"])

higher.to_csv(ROOT / "prime_square_pair_gaussian_higher_order_diagnostics.csv", index=False, encoding="utf-8-sig")
patterns.to_csv(ROOT / "prime_square_pair_gaussian_sign_triples.csv", index=False, encoding="utf-8-sig")

lag001 = lags[(lags["delta"] == 0.001) & (lags["lag_in_interval_units"] == 1)].iloc[0]
lag002 = lags[(lags["delta"] == 0.002) & (lags["lag_in_interval_units"] == 1)].iloc[0]

decision = pd.DataFrame([
    ["Leung 2026 adjacent negative correlation", "weighted psi counts in equal-relative short intervals", "prediction", "rho_1 approximately -log(2)/log(1/delta), with secondary correction from the variance expansion"],
    ["weighted moving-window test delta=0.001", "8000 log-spaced anchors", "observed", f"rho1={lag001['empirical_corr_mean']:.6f}; secondary prediction={lag001['leung_secondary_prediction']:.6f}"],
    ["weighted moving-window test delta=0.002", "8000 log-spaced anchors", "observed", f"rho1={lag002['empirical_corr_mean']:.6f}; secondary prediction={lag002['leung_secondary_prediction']:.6f}"],
    ["prime-square monotone triples", "448 consecutive triples", "observed vs pair-Gaussian", f"observed={obs_monotone}; model mean={mono.mean():.3f}; p={two_sided_p(obs_monotone, mono):.3f}"],
    ["prime-square runs", "450 consecutive intervals", "observed vs pair-Gaussian", f"observed=256; model mean={runs.mean():.3f}; p={two_sided_p(256, runs):.3f}"],
    ["decision on explicit k=3 model", "current evidence", "research decision", "Do not promote triples yet; first make the point-level model reproduce the full pair covariance kernel and test scaling."],
], columns=["item","scope","type","result"])
decision.to_csv(ROOT / "theory_aligned_comparison_summary.csv", index=False, encoding="utf-8-sig")

print(summary.to_string(index=False))
print()
print(lags.to_string(index=False))
print()
print(higher.to_string(index=False))
print()
print(patterns.to_string(index=False))
