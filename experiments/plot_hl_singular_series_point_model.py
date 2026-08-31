#!/usr/bin/env python3
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent
PLOTS = ROOT / "plots"
PLOTS.mkdir(exist_ok=True)

stage = pd.read_csv(ROOT / "hl_point_model_key_results.csv")
pair = pd.read_csv(ROOT / "hl_residual_singular_series_pair_check.csv")

fig, ax = plt.subplots(figsize=(10, 6))
ax.plot(stage["Q"], stage["lag1_mc_mean"], marker="o", label="Point model Monte Carlo mean")
ax.axhline(stage["lag1_observed"].iloc[0], linewidth=1.5, label="Observed")
ax.set_xscale("log")
ax.set_xlabel("Random-sieve cutoff Q after fixed wheel 997")
ax.set_ylabel("Lag-1 autocorrelation")
ax.set_title("Residual singular-series pair structure moves correlation in the right direction")
ax.grid(True, alpha=0.25)
ax.legend()
fig.tight_layout()
fig.savefig(PLOTS / "hl_point_model_lag1_progression.svg")
plt.close(fig)

fig, ax = plt.subplots(figsize=(10, 6))
ax.plot(stage["Q"], stage["R_mc_mean"], marker="o", label="Model mean R")
ax.axhline(stage["R_observed"].iloc[0], linewidth=1.5, label="Observed R")
ax.set_xscale("log")
ax.set_xlabel("Random-sieve cutoff Q after fixed wheel 997")
ax.set_ylabel("R = dispersion relative to q=997 baseline")
ax.set_title("Deeper pair-producing sieve over-regularizes marginal dispersion")
ax.grid(True, alpha=0.25)
ax.legend()
fig.tight_layout()
fig.savefig(PLOTS / "hl_point_model_dispersion_progression.svg")
plt.close(fig)

fig, ax = plt.subplots(figsize=(9, 6))
ax.scatter(pair["theory"], pair["empirical"], s=18, alpha=0.65)
lo = min(pair["theory"].min(), pair["empirical"].min())
hi = max(pair["theory"].max(), pair["empirical"].max())
ax.plot([lo, hi], [lo, hi], linewidth=1)
ax.set_xlabel("Residual singular-series weight after wheel 997")
ax.set_ylabel("Empirical prime-pair ratio among wheel-997 survivor pairs")
ax.set_title("Pair-level signal after conditioning on wheel 997")
ax.grid(True, alpha=0.25)
fig.tight_layout()
fig.savefig(PLOTS / "hl_residual_pair_weight_check.svg")
plt.close(fig)
