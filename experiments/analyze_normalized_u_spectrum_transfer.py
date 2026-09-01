#!/usr/bin/env python3
"""
Normalized-u fine-spectrum transfer experiment.

Each prime-square interval (p_n^2, p_{n+1}^2) is mapped to one unit:

    u = n + (x - p_n^2) / (p_{n+1}^2 - p_n^2).

A dense Haar library is built in u-space. One fine spectrum is fitted
on a single reference panel near X ~ 2e7. The spectrum is then frozen
and transferred to prime-square panels through X ~ 2e12.

No target-panel NNLS refitting is allowed. Only the leading logarithmic
amplitude factor

    alpha = log((X/H)_ref) / log((X/H)_target)

is applied.

Primary model: uniform interval weight in u.
Robustness models: q_i proportional to sqrt(H_i), and to
sqrt(H_i / log(X_i/H_i)).

This is a geometry/theory transfer experiment. It does not count actual
primes up to 1e12; target kernels come from the short-interval covariance law.
"""

import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import nnls
from scipy.sparse import coo_matrix

ROOT = Path(__file__).resolve().parent
X_TARGETS = [2e7, 2e8, 2e9, 2e10, 2e11, 2e12]
PANEL_M = 128
OFFSETS = [-128, 0, 128]
MAX_LAG = 20
CELLS_PER_INTERVAL = 64
STEPS_PER_OCTAVE = 16
REFERENCE_PANEL = "X2e+07_off+0"


def primes_upto(n):
    z = np.ones(n + 1, dtype=bool)
    z[:2] = False
    for p in range(2, math.isqrt(n) + 1):
        if z[p]:
            z[p*p:n+1:p] = False
    return np.flatnonzero(z)


def V(h, X):
    if h <= 0:
        return 0.0
    return h * math.log(X / h)


def pair_rho(A, B, i, j):
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
    return cov / math.sqrt(V(b - a, Xloc) * V(d - c, Xloc))


def log_scales(ncells, steps):
    values = []
    x = 1.0
    xmax = math.log2(ncells)
    while x <= xmax + 1e-9:
        values.append(int(round(2**x)))
        x += 1 / steps
    return np.array(sorted(set(v for v in values if 2 <= v <= ncells)), dtype=int)


P_MAX = int(math.sqrt(max(X_TARGETS)) * 1.05) + 5000
primes = primes_upto(P_MAX)
panel_data = {}

for Xtarget in X_TARGETS:
    center = int(np.argmin(np.abs(primes - math.sqrt(Xtarget))))
    for offset in OFFSETS:
        start_index = center + offset - PANEL_M // 2
        ps = primes[start_index:start_index + PANEL_M + 1].astype(np.int64)
        A = ps[:-1] ** 2
        B = ps[1:] ** 2
        theory = np.array([
            np.mean([pair_rho(A, B, i, i + k) for i in range(PANEL_M - k)])
            for k in range(1, MAX_LAG + 1)
        ])
        panel_id = f"X{Xtarget:.0e}_off{offset:+d}"
        panel_data[panel_id] = (A, B, theory)


NCELL = PANEL_M * CELLS_PER_INTERVAL
SCALES = log_scales(NCELL, STEPS_PER_OCTAVE)
u_cells = np.arange(NCELL, dtype=np.int64)
row_index = np.repeat(np.arange(PANEL_M, dtype=np.int32), CELLS_PER_INTERVAL)
cell_weight = np.full(NCELL, 1 / CELLS_PER_INTERVAL, dtype=float)

gram_scales = []
for L in SCALES:
    shifts = sorted(set([0, int(round(L / 4)), int(round(L / 2)), int(round(3 * L / 4))]))
    Gavg = np.zeros((PANEL_M, PANEL_M), dtype=float)
    for shift in shifts:
        u = u_cells + shift
        block = u // L
        pos = u % L
        sign = np.where(pos < L / 2, 1.0, -1.0)
        W = coo_matrix(
            (cell_weight * sign, (row_index, block)),
            shape=(PANEL_M, int(block.max()) + 1),
        ).tocsr()
        Gavg += (W @ W.T).toarray()
    Gavg /= len(shifts)
    gram_scales.append(Gavg)
gram_scales = np.asarray(gram_scales)


def basis_from_q(q):
    q = np.asarray(q, float)
    basis = np.empty((MAX_LAG + 1, len(SCALES)))
    q2 = q * q
    for j, G in enumerate(gram_scales):
        basis[0, j] = np.mean(q2 * np.diag(G))
        for k in range(1, MAX_LAG + 1):
            basis[k, j] = np.mean(q[:-k] * q[k:] * np.diag(G, k=k))
    return basis


def panel_q(A, B, mode):
    H = (B - A).astype(float)
    centers = 0.5 * (A + B)
    if mode == "uniform":
        q = np.ones(PANEL_M)
    elif mode == "sqrtH":
        q = np.sqrt(H)
        q /= math.sqrt(np.mean(q * q))
    elif mode == "count_over_shortvar":
        q = np.sqrt(H / np.log(centers / H))
        q /= math.sqrt(np.mean(q * q))
    else:
        raise ValueError(mode)
    return q


Aref, Bref, theory_ref = panel_data[REFERENCE_PANEL]
Href = (Bref - Aref).astype(float)
Xref = 0.5 * (Aref + Bref)
Rref = float(np.median(Xref) / np.median(Href))

modes = ["uniform", "sqrtH", "count_over_shortvar"]
reference_rows = []
transfer_rows = []
lag_weights = np.ones(MAX_LAG + 1)
lag_weights[0] = 3.0
lag_weights[1:6] = 2.0
primary_v = None
primary_basis = None

for mode in modes:
    qref = panel_q(Aref, Bref, mode)
    basis_ref = basis_from_q(qref)
    target_ref = np.concatenate([[1.0], theory_ref])
    v, residual = nnls(basis_ref * lag_weights[:, None], target_ref * lag_weights)
    fit_ref = basis_ref @ v
    rho_ref = fit_ref[1:] / fit_ref[0]
    reference_rows.append({
        "mode": mode,
        "nnls_residual": float(residual),
        "kernel_rmse": float(np.sqrt(np.mean((rho_ref - theory_ref) ** 2))),
        "shape_rmse": float(np.sqrt(np.mean((rho_ref / rho_ref[0] - theory_ref / theory_ref[0]) ** 2))),
        "rho1_fit": float(rho_ref[0]),
        "rho1_target": float(theory_ref[0]),
        "nonzero_scales": int(np.sum(v > 1e-12)),
    })
    if mode == "uniform":
        primary_v = v.copy()
        primary_basis = basis_ref.copy()

    for panel_id, (A, B, theory) in panel_data.items():
        q = panel_q(A, B, mode)
        basis_target = basis_from_q(q)
        cov = basis_target @ v
        rho_shape = cov[1:] / cov[0]
        H = (B - A).astype(float)
        X = 0.5 * (A + B)
        Hmed = float(np.median(H))
        Xmed = float(np.median(X))
        R = Xmed / Hmed
        alpha_log = math.log(Rref) / math.log(R)
        prediction = alpha_log * rho_shape
        alpha_best = float(np.dot(rho_shape, theory) / np.dot(rho_shape, rho_shape))
        best_prediction = alpha_best * rho_shape
        transfer_rows.append({
            "mode": mode,
            "panel_id": panel_id,
            "median_X": Xmed,
            "median_H": Hmed,
            "X_over_H": R,
            "alpha_log": alpha_log,
            "alpha_best": alpha_best,
            "alpha_relative_error": (alpha_best - alpha_log) / alpha_log,
            "rho1_pred": float(prediction[0]),
            "rho1_target": float(theory[0]),
            "kernel_rmse_log": float(np.sqrt(np.mean((prediction - theory) ** 2))),
            "shape_rmse": float(np.sqrt(np.mean((prediction / prediction[0] - theory / theory[0]) ** 2))),
            "kernel_rmse_best_scalar": float(np.sqrt(np.mean((best_prediction - theory) ** 2))),
        })

reference = pd.DataFrame(reference_rows)
transfer = pd.DataFrame(transfer_rows)
reference.to_csv(ROOT / "normalized_u_reference_fit.csv", index=False, encoding="utf-8-sig")
transfer.to_csv(ROOT / "normalized_u_fine_spectrum_transfer.csv", index=False, encoding="utf-8-sig")

contribution = primary_v * primary_basis[0]
spectrum = pd.DataFrame({
    "support_cells": SCALES,
    "support_intervals_u": SCALES / CELLS_PER_INTERVAL,
    "coefficient_v": primary_v,
    "variance_contribution": contribution,
})
spectrum["variance_fraction"] = spectrum["variance_contribution"] / spectrum["variance_contribution"].sum()
spectrum.to_csv(ROOT / "normalized_u_reference_spectrum.csv", index=False, encoding="utf-8-sig")

uniform_basis = basis_from_q(np.ones(PANEL_M))
robust_rows = []
reference_candidates = [panel_id for panel_id in panel_data if panel_id.startswith("X2e+07")]

for reference_id in reference_candidates:
    A0, B0, theory0 = panel_data[reference_id]
    v0, _ = nnls(
        uniform_basis * lag_weights[:, None],
        np.concatenate([[1.0], theory0]) * lag_weights,
    )
    cov0 = uniform_basis @ v0
    rho0 = cov0[1:] / cov0[0]
    R0 = float(np.median(0.5 * (A0 + B0)) / np.median(B0 - A0))
    for target_id, (A, B, theory) in panel_data.items():
        R = float(np.median(0.5 * (A + B)) / np.median(B - A))
        alpha = math.log(R0) / math.log(R)
        prediction = alpha * rho0
        robust_rows.append({
            "reference_panel": reference_id,
            "target_panel": target_id,
            "reference_X_over_H": R0,
            "target_X_over_H": R,
            "kernel_rmse": float(np.sqrt(np.mean((prediction - theory) ** 2))),
            "shape_rmse": float(np.sqrt(np.mean((prediction / prediction[0] - theory / theory[0]) ** 2))),
            "rho1_pred": float(prediction[0]),
            "rho1_target": float(theory[0]),
        })

robust = pd.DataFrame(robust_rows)
robust.to_csv(ROOT / "normalized_u_reference_robustness.csv", index=False, encoding="utf-8-sig")
summary_rows = []
for reference_id, group in robust.groupby("reference_panel"):
    held_out = group[group["target_panel"] != reference_id]
    high = held_out[held_out["target_panel"].str.startswith("X2e+12")]
    summary_rows.append({
        "reference_panel": reference_id,
        "held_out_targets": len(held_out),
        "mean_kernel_rmse": float(held_out["kernel_rmse"].mean()),
        "max_kernel_rmse": float(held_out["kernel_rmse"].max()),
        "mean_shape_rmse": float(held_out["shape_rmse"].mean()),
        "max_shape_rmse": float(held_out["shape_rmse"].max()),
        "mean_kernel_rmse_at_2e12": float(high["kernel_rmse"].mean()),
    })
robust_summary = pd.DataFrame(summary_rows)
robust_summary.to_csv(ROOT / "normalized_u_reference_robustness_summary.csv", index=False, encoding="utf-8-sig")

old_path = ROOT / "wide_X_universal_spectrum_transfer.csv"
if old_path.exists():
    old = pd.read_csv(old_path)
    old_mean = float(old["kernel_rmse_theory_alpha"].mean())
    old_max = float(old["kernel_rmse_theory_alpha"].max())
else:
    old_mean = np.nan
    old_max = np.nan

final_rows = []
for mode in modes:
    d = transfer[(transfer["mode"] == mode) & (transfer["panel_id"] != REFERENCE_PANEL)]
    final_rows.append({
        "model": mode,
        "reference_panel": REFERENCE_PANEL,
        "held_out_targets": len(d),
        "mean_kernel_rmse": float(d["kernel_rmse_log"].mean()),
        "max_kernel_rmse": float(d["kernel_rmse_log"].max()),
        "mean_shape_rmse": float(d["shape_rmse"].mean()),
        "max_shape_rmse": float(d["shape_rmse"].max()),
        "mean_abs_alpha_relative_error": float(d["alpha_relative_error"].abs().mean()),
        "max_abs_alpha_relative_error": float(d["alpha_relative_error"].abs().max()),
        "improvement_factor_vs_old_mean": float(old_mean / d["kernel_rmse_log"].mean()) if np.isfinite(old_mean) else np.nan,
    })

if np.isfinite(old_mean):
    final_rows.append({
        "model": "old_physical_median_H_transfer",
        "reference_panel": "previous experiment",
        "held_out_targets": len(old),
        "mean_kernel_rmse": old_mean,
        "max_kernel_rmse": old_max,
        "mean_shape_rmse": np.nan,
        "max_shape_rmse": np.nan,
        "mean_abs_alpha_relative_error": np.nan,
        "max_abs_alpha_relative_error": np.nan,
        "improvement_factor_vs_old_mean": 1.0,
    })

final_summary = pd.DataFrame(final_rows)
final_summary.to_csv(ROOT / "normalized_u_transfer_summary.csv", index=False, encoding="utf-8-sig")

print(reference.to_string(index=False))
print()
print(final_summary.to_string(index=False))
print()
print(robust_summary.to_string(index=False))
