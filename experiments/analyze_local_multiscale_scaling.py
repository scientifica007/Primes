#!/usr/bin/env python3
"""
Dimensionless scaling experiment for the local multiscale Haar representation.

Question
--------
Do the theory-trained local scale contributions behave approximately as a
function of L/H (and X/H), rather than of the absolute scale L alone?

Important methodological point
------------------------------
Individual NNLS coefficients v_s are not identifiable when several Haar
scales produce nearly collinear covariance shapes.  The script therefore
tracks the actual variance contribution

    c_s = v_s * C_s(0)

and aggregates c_s into broad dimensionless L/H bands.

Experiments
-----------
1. Synthetic equal-H controls at X = 20m, 40m, 60m and X/H = 250,500,1000.
2. Spectrum-transfer test: fit at X=20m, scale supports by H_target/H_source,
   and compare to the target kernel without refitting.
3. Four contiguous prime-square panels.
4. Cross-panel transfer from panel 1 to panels 2--4.

Input
-----
first_1000_prime_square_interval_index.csv
"""

import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import nnls
from scipy.sparse import coo_matrix

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "first_1000_prime_square_interval_index.csv"

CELL = 1024
MAX_LAG = 20

df = pd.read_csv(DATA)
s = df[df["p_n"] >= 4001].copy().reset_index(drop=True)

A_all = s["p_n^2"].to_numpy(np.int64)
B_all = s["p_التالي^2"].to_numpy(np.int64)


def Delta(t):
    t = float(t)
    return 0.5 * (
        (t + 1) * math.log(t + 1)
        - 2 * t * math.log(t)
        + (0.0 if t == 1 else (t - 1) * math.log(t - 1))
    )


def V(h, X):
    if h <= 0:
        return 0.0
    return h * math.log(X / h)


def pair_rho_interval(A, B, i, j):
    a, b = float(A[i]), float(B[i])
    c, d = float(A[j]), float(B[j])

    centers = 0.5 * (A + B)
    Xloc = math.sqrt(
        float(centers[i]) * float(centers[j])
    )

    cov = 0.5 * (
        V(d - a, Xloc)
        + V(c - b, Xloc)
        - V(c - a, Xloc)
        - V(d - b, Xloc)
    )

    return cov / math.sqrt(
        V(b - a, Xloc) * V(d - c, Xloc)
    )


def make_log_scales(ncell, steps_per_octave):
    vals = []
    x = 1.0
    maxlog = math.log2(max(2, ncell))

    while x <= maxlog + 1e-9:
        vals.append(int(round(2**x)))
        x += 1 / steps_per_octave

    return np.array(
        sorted(set(
            v for v in vals
            if 2 <= v <= ncell
        )),
        dtype=int
    )


def make_haar_basis(
    A,
    B,
    cell=CELL,
    max_lag=MAX_LAG,
    steps_per_octave=4,
):
    origin = int(A[0])
    starts = (A - origin).astype(np.int64)
    ends = (B - origin).astype(np.int64)
    N = int(ends[-1])
    M = len(A)

    regular = np.arange(
        0, N + cell, cell, dtype=np.int64
    )
    regular = regular[regular <= N]

    if len(regular) == 0 or regular[-1] != N:
        regular = np.append(regular, N)

    boundaries = np.unique(np.concatenate([
        regular, starts, ends, [0, N]
    ])).astype(np.int64)

    seg_starts = boundaries[:-1]
    seg_ends = boundaries[1:]
    seg_lengths = seg_ends - seg_starts

    seg_interval = np.searchsorted(
        ends, seg_starts, side="right"
    ).astype(np.int32)

    seg_cell = (
        seg_starts // cell
    ).astype(np.int64)

    interval_len = (
        ends - starts
    ).astype(float)

    ncell = int(math.ceil(N / cell))

    scales = make_log_scales(
        ncell, steps_per_octave
    )

    rows = []

    for L in scales:
        shifts = sorted(set([
            0,
            int(round(L / 4)),
            int(round(L / 2)),
            int(round(3 * L / 4)),
        ]))

        cov_acc = np.zeros(max_lag + 1)

        for shift in shifts:
            u = seg_cell + shift
            block = u // L
            pos = u % L

            sign = np.where(
                pos < L / 2, 1.0, -1.0
            )

            data = (
                seg_lengths.astype(float)
                / interval_len[seg_interval]
            ) * sign

            W = coo_matrix(
                (
                    data,
                    (seg_interval, block)
                ),
                shape=(
                    M,
                    int(block.max()) + 1
                )
            ).tocsr()

            cov_acc[0] += float(np.mean(
                np.asarray(
                    W.multiply(W).sum(axis=1)
                ).ravel()
            ))

            for k in range(1, max_lag + 1):
                vals = np.asarray(
                    W[:-k]
                    .multiply(W[k:])
                    .sum(axis=1)
                ).ravel()

                cov_acc[k] += float(
                    np.mean(vals)
                )

        cov_acc /= len(shifts)
        rows.append(cov_acc)

    return scales, np.asarray(rows).T


def fit_case(
    A,
    B,
    steps_per_octave=4,
):
    max_lag = min(
        MAX_LAG,
        len(A) - 2
    )

    scales, basis = make_haar_basis(
        A, B,
        max_lag=max_lag,
        steps_per_octave=steps_per_octave
    )

    theory = np.array([
        np.mean([
            pair_rho_interval(
                A, B, i, i + k
            )
            for i in range(len(A) - k)
        ])
        for k in range(1, max_lag + 1)
    ])

    target = np.concatenate([
        [1.0], theory
    ])

    lag_weights = np.ones(
        max_lag + 1
    )
    lag_weights[0] = 3.0
    lag_weights[
        1:min(6, max_lag + 1)
    ] = 2.0

    v, residual = nnls(
        basis * lag_weights[:, None],
        target * lag_weights
    )

    variance_contribution = (
        v * basis[0]
    )

    return {
        "scales": scales,
        "basis": basis,
        "theory": theory,
        "target": target,
        "v": v,
        "contribution": variance_contribution,
        "residual": residual,
    }


def band_fractions(
    scales,
    contribution,
    H,
):
    ratios = scales * CELL / H
    total = contribution.sum()

    definitions = [
        ("L/H<2", 0, 2),
        ("2<=L/H<4", 2, 4),
        ("4<=L/H<8", 4, 8),
        ("L/H>=8", 8, np.inf),
    ]

    result = {}

    for name, lo, hi in definitions:
        mask = (
            (ratios >= lo)
            & (ratios < hi)
        )

        result[name] = float(
            contribution[mask].sum()
            / total
        )

    return result


def synthetic_case(
    ratio,
    X0,
    steps_per_octave=4,
):
    H = int(round(X0 / ratio))
    M = 128

    a0 = int(round(
        X0 - M * H / 2
    ))

    if a0 <= H:
        a0 = int(X0)

    A = np.array([
        a0 + i * H
        for i in range(M)
    ], dtype=np.int64)

    B = A + H

    fitted = fit_case(
        A, B,
        steps_per_octave
    )

    fitted.update({
        "ratio": ratio,
        "X": X0,
        "H": H,
        "A": A,
        "B": B,
    })

    return fitted


# Synthetic dimensionless band collapse.
synthetic_rows = []

for ratio in [250, 500, 1000]:
    for X0 in [20_000_000, 40_000_000, 60_000_000]:
        case = synthetic_case(ratio, X0, 4)
        row = {
            "ratio": ratio,
            "X": X0,
            "H": case["H"],
            "residual": case["residual"],
        }
        row.update(band_fractions(
            case["scales"],
            case["contribution"],
            case["H"]
        ))
        synthetic_rows.append(row)

synthetic = pd.DataFrame(synthetic_rows)
synthetic.to_csv(
    ROOT / "local_multiscale_scaling_synthetic_bands.csv",
    index=False,
    encoding="utf-8-sig"
)


def transfer_spectrum(source, target):
    scale_factor = target["H"] / source["H"]
    target_contribution = np.zeros(len(target["scales"]))

    for L, contribution in zip(source["scales"], source["contribution"]):
        if contribution <= 1e-12:
            continue
        desired = L * scale_factor
        j = int(np.argmin(np.abs(np.log(target["scales"] / desired))))
        target_contribution[j] += contribution

    target_v = np.divide(
        target_contribution,
        target["basis"][0],
        out=np.zeros_like(target_contribution),
        where=target["basis"][0] > 0,
    )

    fit = target["basis"] @ target_v
    correlation = fit[1:] / fit[0]
    kernel_rmse = float(np.sqrt(np.mean((correlation - target["theory"]) ** 2)))

    delta_norm = np.array([
        Delta(k) for k in range(1, len(correlation) + 1)
    ])
    delta_norm /= delta_norm[0]

    shape_rmse = float(np.sqrt(np.mean((correlation / correlation[0] - delta_norm) ** 2)))

    return scale_factor, kernel_rmse, shape_rmse, float(correlation[0]), float(target["theory"][0])


# Spectrum transfer.
transfer_rows = []

for steps in [8]:
    for ratio in [250, 500, 1000]:
        source = synthetic_case(ratio, 20_000_000, steps)
        for Xtarget in [40_000_000, 60_000_000]:
            target = synthetic_case(ratio, Xtarget, steps)
            scale_factor, kernel_rmse, shape_rmse, rho1_pred, rho1_target = transfer_spectrum(source, target)
            transfer_rows.append({
                "X_over_H": ratio,
                "target_X": Xtarget,
                "H_scale_factor": scale_factor,
                "kernel_rmse": kernel_rmse,
                "shape_rmse": shape_rmse,
                "rho1_pred": rho1_pred,
                "rho1_target": rho1_target,
                "scale_steps_per_octave": steps,
            })

source = synthetic_case(500, 20_000_000, 16)
for Xtarget in [40_000_000, 60_000_000]:
    target = synthetic_case(500, Xtarget, 16)
    scale_factor, kernel_rmse, shape_rmse, rho1_pred, rho1_target = transfer_spectrum(source, target)
    transfer_rows.append({
        "X_over_H": 500,
        "target_X": Xtarget,
        "H_scale_factor": scale_factor,
        "kernel_rmse": kernel_rmse,
        "shape_rmse": shape_rmse,
        "rho1_pred": rho1_pred,
        "rho1_target": rho1_target,
        "scale_steps_per_octave": 16,
    })

transfer = pd.DataFrame(transfer_rows)
transfer.to_csv(
    ROOT / "local_multiscale_scaling_synthetic_transfer.csv",
    index=False,
    encoding="utf-8-sig"
)


# Four prime-square panels.
panel_indices = np.array_split(np.arange(len(A_all)), 4)
prime_rows = []
prime_cases = {}

for panel_id, idx in enumerate(panel_indices, 1):
    A = A_all[idx]
    B = B_all[idx]
    case = fit_case(A, B, 4)
    Hmed = float(np.median(B - A))
    Xmed = float(np.median((A + B) / 2))
    row = {
        "panel": panel_id,
        "agg": 1,
        "M": len(A),
        "median_X": Xmed,
        "median_H": Hmed,
        "X_over_H": Xmed / Hmed,
        "residual": case["residual"],
    }
    row.update(band_fractions(case["scales"], case["contribution"], Hmed))
    prime_rows.append(row)
    prime_cases[panel_id] = {
        **fit_case(A, B, 8),
        "A": A,
        "B": B,
        "Hmed": Hmed,
        "Xmed": Xmed,
        "Rmed": Xmed / Hmed,
    }

prime = pd.DataFrame(prime_rows)
prime.to_csv(
    ROOT / "local_multiscale_scaling_prime_panels.csv",
    index=False,
    encoding="utf-8-sig"
)


# Prime-panel transfer from panel 1.
source = prime_cases[1]
source_contribution = source["contribution"]
prime_transfer_rows = []

for panel_id in [2, 3, 4]:
    target = prime_cases[panel_id]
    H_scale = target["Hmed"] / source["Hmed"]
    target_contribution = np.zeros(len(target["scales"]))

    for L, contribution in zip(source["scales"], source_contribution):
        if contribution <= 1e-12:
            continue
        desired = L * H_scale
        j = int(np.argmin(np.abs(np.log(target["scales"] / desired))))
        target_contribution[j] += contribution

    target_v = np.divide(
        target_contribution,
        target["basis"][0],
        out=np.zeros_like(target_contribution),
        where=target["basis"][0] > 0,
    )

    fit = target["basis"] @ target_v
    transferred_rho = fit[1:] / fit[0]
    amplitude_factor = math.log(source["Rmed"]) / math.log(target["Rmed"])
    prediction = amplitude_factor * transferred_rho

    prime_transfer_rows.append({
        "source_panel": 1,
        "target_panel": panel_id,
        "source_X_over_H": source["Rmed"],
        "target_X_over_H": target["Rmed"],
        "H_scale_factor": H_scale,
        "log_amplitude_factor": amplitude_factor,
        "kernel_rmse_log_factor": float(np.sqrt(np.mean((prediction - target["theory"]) ** 2))),
        "shape_rmse_log_factor": float(np.sqrt(np.mean((prediction / prediction[0] - target["theory"] / target["theory"][0]) ** 2))),
        "rho1_pred_log_factor": float(prediction[0]),
        "rho1_target": float(target["theory"][0]),
    })

prime_transfer = pd.DataFrame(prime_transfer_rows)
prime_transfer.to_csv(
    ROOT / "local_multiscale_scaling_prime_transfer.csv",
    index=False,
    encoding="utf-8-sig"
)

print(synthetic.to_string(index=False))
print()
print(transfer.to_string(index=False))
print()
print(prime.to_string(index=False))
print()
print(prime_transfer.to_string(index=False))
